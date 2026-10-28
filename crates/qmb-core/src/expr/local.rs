use alloc::{sync::Arc, vec::Vec};

use super::parts::{ActExpr, AdjointExpr, ProductExpr, ScaleExpr, SumExpr};
use crate::{ExpressionError, Identity, LinearCombination};

/// A shared local state over state primitive S, operator primitive O and coefficient C.
///
/// Cloning shares immutable nodes without cloning primitives or coefficients.
/// Construction requests no algebra or numeric capabilities.
///
/// Different state primitives cannot be added implicitly:
///
/// ```compile_fail
/// use qmb_core::LocalState;
/// let a = LocalState::<u8, (), ()>::primitive(0);
/// let b = LocalState::<u16, (), ()>::primitive(0);
/// let _ = a + b;
/// ```
///
/// Coefficient precision also stays fixed:
///
/// ```compile_fail
/// use qmb_core::LocalState;
/// let a = LocalState::<(), (), f32>::zero();
/// let b = LocalState::<(), (), f64>::zero();
/// let _ = a + b;
/// ```
///
/// Shared ownership does not make a non-Send primitive transferable:
///
/// ```compile_fail
/// use std::rc::Rc;
/// use qmb_core::LocalState;
/// fn needs_send<T: Send>() {}
/// needs_send::<LocalState<Rc<()>, (), ()>>();
/// ```
pub struct LocalState<S, O, C> {
    pub(crate) root: Option<Arc<StateNode<S, O, C>>>,
    // A structural property, not an evaluation cache. It permits O(1) reuse of
    // already evaluated subtrees without rescanning long expression chains.
    pub(crate) has_action: bool,
}

/// A shared local operator over primitive O and coefficient C.
///
/// Its primitive type must match the state's delayed action:
///
/// ```compile_fail
/// use qmb_core::{LocalOperator, LocalState};
/// let state = LocalState::<(), u8, ()>::zero();
/// let other = LocalOperator::<u16, ()>::primitive(0);
/// let _ = state.act(&other);
/// ```
pub struct LocalOperator<O, C> {
    pub(crate) root: Option<Arc<OperatorNode<O, C>>>,
}

pub(crate) enum StateNode<S, O, C> {
    Primitive(S),
    Scale(ScaleExpr<LocalState<S, O, C>, C>),
    Sum(SumExpr<LocalState<S, O, C>>),
    Act(ActExpr<LocalState<S, O, C>, LocalOperator<O, C>>),
}

pub(crate) enum OperatorNode<O, C> {
    Primitive(O),
    Scale(ScaleExpr<LocalOperator<O, C>, C>),
    Sum(SumExpr<LocalOperator<O, C>>),
    Product(ProductExpr<LocalOperator<O, C>>),
    Adjoint(AdjointExpr<LocalOperator<O, C>>),
}

impl<S, O, C> Clone for LocalState<S, O, C> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
            has_action: self.has_action,
        }
    }
}
impl<O, C> Clone for LocalOperator<O, C> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
        }
    }
}

impl<S, O, C> LocalState<S, O, C> {
    pub(crate) fn has_pending_action(&self) -> bool {
        self.has_action
    }
    pub(crate) fn new(node: StateNode<S, O, C>) -> Self {
        let has_action = match &node {
            StateNode::Primitive(_) => false,
            StateNode::Scale(s) => s.operand.has_action,
            StateNode::Sum(s) => s.iter().any(|s| s.has_action),
            StateNode::Act(_) => true,
        };
        Self {
            root: Some(Arc::new(node)),
            has_action,
        }
    }
    pub(crate) fn node(&self) -> &StateNode<S, O, C> {
        self.root.as_deref().expect("live state handle")
    }

    /// Stores a primitive without generating a coefficient or running a rule.
    pub fn primitive(primitive: S) -> Self {
        Self::new(StateNode::Primitive(primitive))
    }
    /// Constructs the empty sum, without numeric bounds.
    pub fn zero() -> Self {
        Self::new(StateNode::Sum(SumExpr::Empty))
    }
    /// Constructs an ordered sum; one operand is returned unchanged.
    pub fn sum(operands: impl IntoIterator<Item = Self>) -> Self {
        Self::sum_collected(operands.into_iter().collect())
    }
    fn sum_collected(operands: SumExpr<Self>) -> Self {
        match operands {
            SumExpr::One(value) => value,
            values => Self::new(StateNode::Sum(values)),
        }
    }
    /// Collects fallible operands in order, stopping at the first input error.
    /// No primitive rules are evaluated and errors are returned unchanged.
    pub fn try_sum(
        operands: impl IntoIterator<Item = Result<Self, ExpressionError>>,
    ) -> Result<Self, ExpressionError> {
        Ok(Self::sum_collected(
            operands.into_iter().collect::<Result<_, _>>()?,
        ))
    }
    /// Incorporates mapped terms as ordinary scale and sum nodes.
    pub fn from_linear_combination(terms: LinearCombination<Self, C>) -> Self {
        Self::sum(
            terms
                .into_iter()
                .map(|term| term.operand.into_scaled(term.coefficient)),
        )
    }
    pub(crate) fn into_scaled(self, coefficient: C) -> Self {
        Self::new(StateNode::Scale(ScaleExpr {
            coefficient,
            operand: self,
        }))
    }
    /// Constructs a scale node sharing this expression, even for zero or one.
    /// Converts the input to C exactly once without changing coefficient precision.
    pub fn scaled(&self, coefficient: impl Into<C>) -> Self {
        self.clone().into_scaled(coefficient.into())
    }
    pub(crate) fn into_act(self, operator: LocalOperator<O, C>) -> Self {
        Self::new(StateNode::Act(ActExpr {
            state: self,
            operator,
        }))
    }
    /// Constructs a delayed left action. `psi.act(b).act(a)` means `A(B psi)`.
    pub fn act(&self, operator: &LocalOperator<O, C>) -> Self {
        self.clone().into_act(operator.clone())
    }
}

impl<O, C> LocalOperator<O, C> {
    pub(crate) fn new(node: OperatorNode<O, C>) -> Self {
        Self {
            root: Some(Arc::new(node)),
        }
    }
    pub(crate) fn node(&self) -> &OperatorNode<O, C> {
        self.root.as_deref().expect("live operator handle")
    }
    /// Stores a primitive without generating a coefficient or running a rule.
    pub fn primitive(primitive: O) -> Self {
        Self::new(OperatorNode::Primitive(primitive))
    }
    /// Constructs the empty sum, without numeric bounds.
    pub fn zero() -> Self {
        Self::new(OperatorNode::Sum(SumExpr::Empty))
    }
    /// Constructs an ordered sum; one operand is returned unchanged.
    pub fn sum(operands: impl IntoIterator<Item = Self>) -> Self {
        Self::sum_collected(operands.into_iter().collect())
    }
    fn sum_collected(operands: SumExpr<Self>) -> Self {
        match operands {
            SumExpr::One(value) => value,
            values => Self::new(OperatorNode::Sum(values)),
        }
    }
    /// Collects fallible operands in order, stopping at the first input error.
    /// No primitive rules are evaluated and errors are returned unchanged.
    pub fn try_sum(
        operands: impl IntoIterator<Item = Result<Self, ExpressionError>>,
    ) -> Result<Self, ExpressionError> {
        Ok(Self::sum_collected(
            operands.into_iter().collect::<Result<_, _>>()?,
        ))
    }
    /// Incorporates mapped terms as ordinary scale and sum nodes.
    pub fn from_linear_combination(terms: LinearCombination<Self, C>) -> Self {
        Self::sum(
            terms
                .into_iter()
                .map(|term| term.operand.into_scaled(term.coefficient)),
        )
    }
    pub(crate) fn into_scaled(self, coefficient: C) -> Self {
        Self::new(OperatorNode::Scale(ScaleExpr {
            coefficient,
            operand: self,
        }))
    }
    /// Constructs a scale node sharing this expression, even for zero or one.
    /// Converts the input to C exactly once without changing coefficient precision.
    pub fn scaled(&self, coefficient: impl Into<C>) -> Self {
        self.clone().into_scaled(coefficient.into())
    }
    /// Constructs the identity primitive without generating numeric coefficients.
    pub fn identity() -> Self
    where
        O: Identity,
    {
        <Self as Identity>::identity()
    }
    pub(crate) fn into_product(self, rhs: Self) -> Self {
        Self::new(OperatorNode::Product(ProductExpr { lhs: self, rhs }))
    }
    /// Constructs the ordered product `self * rhs`; `rhs` acts first.
    pub fn product(&self, rhs: &Self) -> Self {
        self.clone().into_product(rhs.clone())
    }
    /// Constructs a product in input order. The empty product is identity.
    pub fn product_all(operands: impl IntoIterator<Item = Self>) -> Self
    where
        O: Identity,
    {
        let mut iter = operands.into_iter();
        match iter.next() {
            None => Self::identity(),
            Some(first) => iter.fold(first, Self::into_product),
        }
    }
    pub(crate) fn into_adjoint(self) -> Self {
        Self::new(OperatorNode::Adjoint(AdjointExpr { operand: self }))
    }
    /// Constructs a delayed adjoint, without distribution or primitive rules.
    pub fn adjoint(&self) -> Self {
        self.clone().into_adjoint()
    }
    /// Constructs `self * rhs + rhs * self`.
    pub fn anticommutator(&self, rhs: &Self) -> Self {
        Self::sum([self.product(rhs), rhs.product(self)])
    }
}

impl<S, O, C> Default for LocalState<S, O, C> {
    fn default() -> Self {
        Self::zero()
    }
}
impl<O, C> Default for LocalOperator<O, C> {
    fn default() -> Self {
        Self::zero()
    }
}

enum Pending<S, O, C> {
    State(Arc<StateNode<S, O, C>>),
    Operator(Arc<OperatorNode<O, C>>),
}

fn state_child<S, O, C>(mut child: LocalState<S, O, C>, stack: &mut Vec<Pending<S, O, C>>) {
    if let Some(root) = child.root.take() {
        stack.push(Pending::State(root));
    }
}
fn operator_child<S, O, C>(mut child: LocalOperator<O, C>, stack: &mut Vec<Pending<S, O, C>>) {
    if let Some(root) = child.root.take() {
        stack.push(Pending::Operator(root));
    }
}

// Arc::into_inner also avoids recursive final-release races between threads.
// Taking children out before their handles drop keeps the whole mixed DAG on
// this work stack. Empty/shared/primitive roots need no work-stack allocation.
fn release<S, O, C>(root: Pending<S, O, C>) {
    let mut stack = Vec::new();
    let mut first = Some(root);
    while let Some(root) = first.take().or_else(|| stack.pop()) {
        match root {
            Pending::State(root) => {
                if let Some(node) = Arc::into_inner(root) {
                    match node {
                        StateNode::Primitive(_) => {}
                        StateNode::Scale(s) => state_child(s.operand, &mut stack),
                        StateNode::Sum(s) => {
                            for child in s {
                                state_child(child, &mut stack)
                            }
                        }
                        StateNode::Act(a) => {
                            state_child(a.state, &mut stack);
                            operator_child(a.operator, &mut stack);
                        }
                    }
                }
            }
            Pending::Operator(root) => {
                if let Some(node) = Arc::into_inner(root) {
                    match node {
                        OperatorNode::Primitive(_) => {}
                        OperatorNode::Scale(s) => operator_child(s.operand, &mut stack),
                        OperatorNode::Sum(s) => {
                            for child in s {
                                operator_child(child, &mut stack)
                            }
                        }
                        OperatorNode::Product(p) => {
                            operator_child(p.lhs, &mut stack);
                            operator_child(p.rhs, &mut stack);
                        }
                        OperatorNode::Adjoint(a) => operator_child(a.operand, &mut stack),
                    }
                }
            }
        }
    }
}

impl<S, O, C> Drop for LocalState<S, O, C> {
    fn drop(&mut self) {
        if let Some(root) = self.root.take() {
            release(Pending::State(root))
        }
    }
}
impl<O, C> Drop for LocalOperator<O, C> {
    fn drop(&mut self) {
        if let Some(root) = self.root.take() {
            // Operators have no state descendants. The mixed work stack needs
            // no concrete state primitive when releasing an operator root.
            release::<core::convert::Infallible, O, C>(Pending::Operator(root))
        }
    }
}

impl<O: Identity, C> Identity for OperatorNode<O, C> {
    fn identity() -> Self {
        Self::Primitive(O::identity())
    }
}

impl<O: Identity, C> Identity for LocalOperator<O, C> {
    fn identity() -> Self {
        Self::new(OperatorNode::identity())
    }
}
