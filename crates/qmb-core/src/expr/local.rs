use alloc::{sync::Arc, vec::Vec};

use super::parts::{ActExpr, AdjointExpr, ProductExpr, ScaleExpr, SumExpr};
use crate::{LinearCombination, LocalStructure};

/// An immutable local state expression with shared subexpressions.
///
/// Cloning shares nodes and requires no cloning or numeric capabilities of the
/// structure, primitives, or coefficients. Operations preserve their inputs.
///
/// Coefficient precisions cannot be mixed implicitly:
///
/// ```compile_fail
/// use core::num::NonZeroUsize;
/// use num_complex::Complex;
/// use qmb_core::{LocalState, LocalStructure};
/// struct L;
/// impl LocalStructure for L {
///     type State = u8;
///     type Operator = u8;
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(2).unwrap();
/// }
/// let a = LocalState::<L, Complex<f32>>::primitive(0);
/// let b = LocalState::<L, Complex<f64>>::primitive(0);
/// let _ = a + b;
/// ```
///
/// Shared ownership does not make non-Send primitives transferable:
///
/// ```compile_fail
/// use core::num::NonZeroUsize;
/// use std::rc::Rc;
/// use num_complex::Complex;
/// use qmb_core::{LocalState, LocalStructure};
/// struct L;
/// impl LocalStructure for L {
///     type State = Rc<()>;
///     type Operator = ();
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
/// }
/// fn needs_send<C: Send>() {}
/// needs_send::<LocalState<L, Complex<f64>>>();
/// ```
pub struct LocalState<L: LocalStructure, C> {
    pub(super) root: Option<Arc<StateNode<L, C>>>,
    // A structural property, not an evaluation cache. It permits O(1) reuse of
    // already evaluated subtrees without rescanning long expression chains.
    pub(super) has_action: bool,
}

/// An immutable local operator expression, independent of a chosen basis.
///
/// Distinct local structures stay distinct even if they use the same primitives:
///
/// ```compile_fail
/// use core::num::NonZeroUsize;
/// use num_complex::Complex;
/// use qmb_core::{LocalState, LocalOperator, LocalStructure};
/// struct A;
/// struct B;
/// impl LocalStructure for A {
///     type State = u8;
///     type Operator = u8;
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(2).unwrap();
/// }
/// impl LocalStructure for B {
///     type State = u8;
///     type Operator = u8;
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(2).unwrap();
/// }
/// let state = LocalState::<A, Complex<f64>>::primitive(0);
/// let op = LocalOperator::<B, Complex<f64>>::primitive(0);
/// let _ = state.act(&op);
/// ```
pub struct LocalOperator<L: LocalStructure, C> {
    pub(super) root: Option<Arc<OperatorNode<L, C>>>,
}

pub(super) enum StateNode<L: LocalStructure, C> {
    Primitive(L::State),
    Scale(ScaleExpr<LocalState<L, C>, C>),
    Sum(SumExpr<LocalState<L, C>>),
    Act(ActExpr<LocalState<L, C>, LocalOperator<L, C>>),
}

pub(super) enum OperatorNode<L: LocalStructure, C> {
    Primitive(L::Operator),
    Identity,
    Scale(ScaleExpr<LocalOperator<L, C>, C>),
    Sum(SumExpr<LocalOperator<L, C>>),
    Product(ProductExpr<LocalOperator<L, C>>),
    Adjoint(AdjointExpr<LocalOperator<L, C>>),
}

impl<L: LocalStructure, C> Clone for LocalState<L, C> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
            has_action: self.has_action,
        }
    }
}
impl<L: LocalStructure, C> Clone for LocalOperator<L, C> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
        }
    }
}

impl<L: LocalStructure, C> LocalState<L, C> {
    pub(crate) fn has_pending_action(&self) -> bool {
        self.has_action
    }
    pub(super) fn new(node: StateNode<L, C>) -> Self {
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
    pub(super) fn node(&self) -> &StateNode<L, C> {
        self.root.as_deref().expect("live state handle")
    }

    /// Stores a primitive without generating a coefficient or running a rule.
    pub fn primitive(primitive: L::State) -> Self {
        Self::new(StateNode::Primitive(primitive))
    }
    /// Constructs the empty sum, without numeric bounds.
    pub fn zero() -> Self {
        Self::new(StateNode::Sum(SumExpr::Empty))
    }
    /// Constructs an ordered sum; one operand is returned unchanged.
    pub fn sum(operands: impl IntoIterator<Item = Self>) -> Self {
        match operands.into_iter().collect() {
            SumExpr::One(value) => value,
            values => Self::new(StateNode::Sum(values)),
        }
    }
    /// Incorporates mapped terms as ordinary scale and sum nodes.
    pub fn from_linear_combination(terms: LinearCombination<Self, C>) -> Self {
        Self::sum(
            terms
                .into_iter()
                .map(|term| term.operand.into_scaled(term.coefficient)),
        )
    }
    pub(super) fn into_scaled(self, coefficient: C) -> Self {
        Self::new(StateNode::Scale(ScaleExpr {
            coefficient,
            operand: self,
        }))
    }
    /// Constructs a scale node sharing this expression, even for zero or one.
    pub fn scaled(&self, coefficient: C) -> Self {
        self.clone().into_scaled(coefficient)
    }
    pub(super) fn into_act(self, operator: LocalOperator<L, C>) -> Self {
        Self::new(StateNode::Act(ActExpr {
            state: self,
            operator,
        }))
    }
    /// Constructs a delayed left action. `psi.act(b).act(a)` means `A(B psi)`.
    pub fn act(&self, operator: &LocalOperator<L, C>) -> Self {
        self.clone().into_act(operator.clone())
    }
}

impl<L: LocalStructure, C> LocalOperator<L, C> {
    pub(super) fn new(node: OperatorNode<L, C>) -> Self {
        Self {
            root: Some(Arc::new(node)),
        }
    }
    pub(super) fn node(&self) -> &OperatorNode<L, C> {
        self.root.as_deref().expect("live operator handle")
    }
    /// Stores a primitive without generating a coefficient or running a rule.
    pub fn primitive(primitive: L::Operator) -> Self {
        Self::new(OperatorNode::Primitive(primitive))
    }
    /// Constructs the empty sum, without numeric bounds.
    pub fn zero() -> Self {
        Self::new(OperatorNode::Sum(SumExpr::Empty))
    }
    /// Constructs an ordered sum; one operand is returned unchanged.
    pub fn sum(operands: impl IntoIterator<Item = Self>) -> Self {
        match operands.into_iter().collect() {
            SumExpr::One(value) => value,
            values => Self::new(OperatorNode::Sum(values)),
        }
    }
    /// Incorporates mapped terms as ordinary scale and sum nodes.
    pub fn from_linear_combination(terms: LinearCombination<Self, C>) -> Self {
        Self::sum(
            terms
                .into_iter()
                .map(|term| term.operand.into_scaled(term.coefficient)),
        )
    }
    pub(super) fn into_scaled(self, coefficient: C) -> Self {
        Self::new(OperatorNode::Scale(ScaleExpr {
            coefficient,
            operand: self,
        }))
    }
    /// Constructs a scale node sharing this expression, even for zero or one.
    pub fn scaled(&self, coefficient: C) -> Self {
        self.clone().into_scaled(coefficient)
    }
    /// Constructs a symbolic identity, without invoking the identity rule.
    pub fn identity() -> Self {
        Self::new(OperatorNode::Identity)
    }
    pub(super) fn into_product(self, rhs: Self) -> Self {
        Self::new(OperatorNode::Product(ProductExpr { lhs: self, rhs }))
    }
    /// Constructs the ordered product `self * rhs`; `rhs` acts first.
    pub fn product(&self, rhs: &Self) -> Self {
        self.clone().into_product(rhs.clone())
    }
    /// Constructs a product in input order. The empty product is identity.
    pub fn product_all(operands: impl IntoIterator<Item = Self>) -> Self {
        let mut iter = operands.into_iter();
        match iter.next() {
            None => Self::identity(),
            Some(first) => iter.fold(first, Self::into_product),
        }
    }
    pub(super) fn into_adjoint(self) -> Self {
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

impl<L: LocalStructure, C> Default for LocalState<L, C> {
    fn default() -> Self {
        Self::zero()
    }
}
impl<L: LocalStructure, C> Default for LocalOperator<L, C> {
    fn default() -> Self {
        Self::zero()
    }
}

enum Pending<L: LocalStructure, C> {
    State(Arc<StateNode<L, C>>),
    Operator(Arc<OperatorNode<L, C>>),
}

fn state_child<L: LocalStructure, C>(mut child: LocalState<L, C>, stack: &mut Vec<Pending<L, C>>) {
    if let Some(root) = child.root.take() {
        stack.push(Pending::State(root));
    }
}
fn operator_child<L: LocalStructure, C>(
    mut child: LocalOperator<L, C>,
    stack: &mut Vec<Pending<L, C>>,
) {
    if let Some(root) = child.root.take() {
        stack.push(Pending::Operator(root));
    }
}

// Arc::into_inner also avoids recursive final-release races between threads.
// Taking children out before their handles drop keeps the whole mixed DAG on
// this work stack. Empty/shared/primitive roots need no work-stack allocation.
fn release<L: LocalStructure, C>(root: Pending<L, C>) {
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
                        OperatorNode::Primitive(_) | OperatorNode::Identity => {}
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

impl<L: LocalStructure, C> Drop for LocalState<L, C> {
    fn drop(&mut self) {
        if let Some(root) = self.root.take() {
            release(Pending::State(root))
        }
    }
}
impl<L: LocalStructure, C> Drop for LocalOperator<L, C> {
    fn drop(&mut self) {
        if let Some(root) = self.root.take() {
            release(Pending::Operator(root))
        }
    }
}
