use alloc::{sync::Arc, vec::Vec};

use super::parts::{ActExpr, AdjointExpr, ProductExpr, ScaleExpr, SumExpr};
use crate::{Complex, LinearCombination, LocalStructure};

/// An immutable local state expression with shared subexpressions.
///
/// Cloning shares nodes and requires no cloning or numeric capabilities of the
/// structure, primitives, or coefficients. Operations preserve their inputs.
///
/// Coefficient precisions cannot be mixed implicitly:
///
/// ```compile_fail
/// use core::num::NonZeroUsize;
/// use qmb_core::{LocalState, LocalStructure};
/// struct L;
/// impl LocalStructure for L {
///     type State = u8;
///     type Operator = u8;
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(2).unwrap();
/// }
/// let a = LocalState::<L, f32>::primitive(0);
/// let b = LocalState::<L, f64>::primitive(0);
/// let _ = a + b;
/// ```
///
/// Shared ownership does not make non-Send primitives transferable:
///
/// ```compile_fail
/// use core::num::NonZeroUsize;
/// use std::rc::Rc;
/// use qmb_core::{LocalState, LocalStructure};
/// struct L;
/// impl LocalStructure for L {
///     type State = Rc<()>;
///     type Operator = ();
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
/// }
/// fn needs_send<T: Send>() {}
/// needs_send::<LocalState<L, f64>>();
/// ```
pub struct LocalState<L: LocalStructure, T> {
    pub(super) root: Option<Arc<StateNode<L, T>>>,
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
/// let state = LocalState::<A, f64>::primitive(0);
/// let op = LocalOperator::<B, f64>::primitive(0);
/// let _ = state.act(&op);
/// ```
pub struct LocalOperator<L: LocalStructure, T> {
    pub(super) root: Option<Arc<OperatorNode<L, T>>>,
}

pub(super) enum StateNode<L: LocalStructure, T> {
    Primitive(L::State),
    Scale(ScaleExpr<LocalState<L, T>, T>),
    Sum(SumExpr<LocalState<L, T>>),
    Act(ActExpr<LocalState<L, T>, LocalOperator<L, T>>),
}

pub(super) enum OperatorNode<L: LocalStructure, T> {
    Primitive(L::Operator),
    Identity,
    Scale(ScaleExpr<LocalOperator<L, T>, T>),
    Sum(SumExpr<LocalOperator<L, T>>),
    Product(ProductExpr<LocalOperator<L, T>>),
    Adjoint(AdjointExpr<LocalOperator<L, T>>),
}

impl<L: LocalStructure, T> Clone for LocalState<L, T> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
            has_action: self.has_action,
        }
    }
}
impl<L: LocalStructure, T> Clone for LocalOperator<L, T> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
        }
    }
}

impl<L: LocalStructure, T> LocalState<L, T> {
    pub(super) fn new(node: StateNode<L, T>) -> Self {
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
    pub(super) fn node(&self) -> &StateNode<L, T> {
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
    pub fn from_linear_combination(terms: LinearCombination<Self, T>) -> Self {
        Self::sum(
            terms
                .into_iter()
                .map(|term| term.operand.into_scaled(term.coefficient)),
        )
    }
    pub(super) fn into_scaled(self, coefficient: Complex<T>) -> Self {
        Self::new(StateNode::Scale(ScaleExpr {
            coefficient,
            operand: self,
        }))
    }
    /// Constructs a scale node sharing this expression, even for zero or one.
    pub fn scaled(&self, coefficient: Complex<T>) -> Self {
        self.clone().into_scaled(coefficient)
    }
    pub(super) fn into_act(self, operator: LocalOperator<L, T>) -> Self {
        Self::new(StateNode::Act(ActExpr {
            state: self,
            operator,
        }))
    }
    /// Constructs a delayed left action. `psi.act(b).act(a)` means `A(B psi)`.
    pub fn act(&self, operator: &LocalOperator<L, T>) -> Self {
        self.clone().into_act(operator.clone())
    }
}

impl<L: LocalStructure, T> LocalOperator<L, T> {
    pub(super) fn new(node: OperatorNode<L, T>) -> Self {
        Self {
            root: Some(Arc::new(node)),
        }
    }
    pub(super) fn node(&self) -> &OperatorNode<L, T> {
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
    pub fn from_linear_combination(terms: LinearCombination<Self, T>) -> Self {
        Self::sum(
            terms
                .into_iter()
                .map(|term| term.operand.into_scaled(term.coefficient)),
        )
    }
    pub(super) fn into_scaled(self, coefficient: Complex<T>) -> Self {
        Self::new(OperatorNode::Scale(ScaleExpr {
            coefficient,
            operand: self,
        }))
    }
    /// Constructs a scale node sharing this expression, even for zero or one.
    pub fn scaled(&self, coefficient: Complex<T>) -> Self {
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

impl<L: LocalStructure, T> Default for LocalState<L, T> {
    fn default() -> Self {
        Self::zero()
    }
}
impl<L: LocalStructure, T> Default for LocalOperator<L, T> {
    fn default() -> Self {
        Self::zero()
    }
}

enum Pending<L: LocalStructure, T> {
    State(Arc<StateNode<L, T>>),
    Operator(Arc<OperatorNode<L, T>>),
}

fn state_child<L: LocalStructure, T>(mut child: LocalState<L, T>, stack: &mut Vec<Pending<L, T>>) {
    if let Some(root) = child.root.take() {
        stack.push(Pending::State(root));
    }
}
fn operator_child<L: LocalStructure, T>(
    mut child: LocalOperator<L, T>,
    stack: &mut Vec<Pending<L, T>>,
) {
    if let Some(root) = child.root.take() {
        stack.push(Pending::Operator(root));
    }
}

// Arc::into_inner also avoids recursive final-release races between threads.
// Taking children out before their handles drop keeps the whole mixed DAG on
// this work stack. Empty/shared/primitive roots need no work-stack allocation.
fn release<L: LocalStructure, T>(root: Pending<L, T>) {
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

impl<L: LocalStructure, T> Drop for LocalState<L, T> {
    fn drop(&mut self) {
        if let Some(root) = self.root.take() {
            release(Pending::State(root))
        }
    }
}
impl<L: LocalStructure, T> Drop for LocalOperator<L, T> {
    fn drop(&mut self) {
        if let Some(root) = self.root.take() {
            release(Pending::Operator(root))
        }
    }
}
