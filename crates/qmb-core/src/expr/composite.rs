//! Placed states and operators on ordinary tensor spaces.

use alloc::{sync::Arc, vec::Vec};

use super::parts::{ActExpr, AdjointExpr, ProductExpr, ScaleExpr, SumExpr};
use crate::{HilbertSpace, LocalOperator, LocalState, LocalStructure, SpaceError};

/// Errors connecting, placing or evaluating composite expressions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum ExpressionError {
    /// A local expression must be placed on a leaf, not an arbitrary tensor.
    #[error("local placement requires a factor space")]
    NotAFactor,
    /// Operands have incompatible ordered factor spaces.
    #[error("incompatible Hilbert spaces")]
    IncompatibleSpaces,
    /// An operator refers to an index absent from its selected ambient space.
    #[error("referenced index is absent from the ambient space")]
    MissingIndex,
    /// No binding or explicit ambient space is available.
    #[error("operator ambient space is not specified")]
    UnboundOperator,
    /// An invalid tensor space was requested.
    #[error(transparent)]
    Space(#[from] SpaceError),
}

/// A composite state with an immutable, always definite Hilbert space.
///
/// `I` supplies index equality; `L` fixes the local algebra and `C` is the whole
/// coefficient. Construction and sharing require no arithmetic or primitive Clone.
///
/// Different coefficient precisions cannot be added:
///
/// ```compile_fail
/// use core::num::NonZeroUsize;
/// use num_complex::Complex;
/// use qmb_core::{HilbertSpace, LocalStructure, State};
/// struct L;
/// impl LocalStructure for L {
///     type State = ();
///     type Operator = ();
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
/// }
/// let factor = HilbertSpace::<u8, L>::factor(0);
/// let single = State::<_, _, Complex<f32>>::zero(&factor);
/// let double = State::<_, _, Complex<f64>>::zero(&factor);
/// let _ = single + double;
/// ```
///
/// Non-Send indices are supported locally, without making the state Send:
///
/// ```compile_fail
/// use core::num::NonZeroUsize;
/// use std::rc::Rc;
/// use num_complex::Complex;
/// use qmb_core::{LocalStructure, State};
/// struct L;
/// impl LocalStructure for L {
///     type State = ();
///     type Operator = ();
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
/// }
/// fn needs_send<C: Send>() {}
/// needs_send::<State<Rc<()>, L, Complex<f64>>>();
/// ```
pub struct State<I, L: LocalStructure, C> {
    pub(super) expression: StateExpr<I, L, C>,
    pub(super) space: HilbertSpace<I, L>,
}

/// A placed operator with an optional ambient-space binding.
///
/// Bindings are held only by this outer handle. Combining operators shares
/// their syntax without retaining old bindings inside child nodes.
///
/// Checked binary operators return `Result`; named methods take shared borrows.
/// Scaling and adjoints are lazy and preserve an existing binding. Evaluation
/// of an unbound operator requires an explicit space, except that action can
/// use the state's space. See [`State::act`] and [`Operator::trace_in`].
pub struct Operator<I, L: LocalStructure, C> {
    pub(super) expression: OpExpr<I, L, C>,
    pub(super) ambient: Option<HilbertSpace<I, L>>,
}

pub(super) struct StateExpr<I, L: LocalStructure, C> {
    root: Option<Arc<StateNode<I, L, C>>>,
    pub(super) pending_action: bool,
}
pub(super) struct OpExpr<I, L: LocalStructure, C> {
    root: Option<Arc<OpNode<I, L, C>>>,
}

pub(super) enum StateNode<I, L: LocalStructure, C> {
    Local(HilbertSpace<I, L>, LocalState<L, C>),
    Unit,
    Scale(ScaleExpr<StateExpr<I, L, C>, C>),
    Sum(SumExpr<StateExpr<I, L, C>>),
    Tensor(SumExpr<StateExpr<I, L, C>>),
    Act(ActExpr<StateExpr<I, L, C>, OpExpr<I, L, C>>),
}
pub(super) enum OpNode<I, L: LocalStructure, C> {
    Local(HilbertSpace<I, L>, LocalOperator<L, C>),
    Identity,
    Scale(ScaleExpr<OpExpr<I, L, C>, C>),
    Sum(SumExpr<OpExpr<I, L, C>>),
    Product(ProductExpr<OpExpr<I, L, C>>),
    Adjoint(AdjointExpr<OpExpr<I, L, C>>),
}

impl<I, L: LocalStructure, C> Clone for StateExpr<I, L, C> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
            pending_action: self.pending_action,
        }
    }
}
impl<I, L: LocalStructure, C> Clone for OpExpr<I, L, C> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
        }
    }
}
impl<I, L: LocalStructure, C> Clone for State<I, L, C> {
    fn clone(&self) -> Self {
        Self {
            expression: self.expression.clone(),
            space: self.space.clone(),
        }
    }
}
impl<I, L: LocalStructure, C> Clone for Operator<I, L, C> {
    fn clone(&self) -> Self {
        Self {
            expression: self.expression.clone(),
            ambient: self.ambient.clone(),
        }
    }
}

impl<I, L: LocalStructure, C> StateExpr<I, L, C> {
    pub(super) fn new(node: StateNode<I, L, C>) -> Self {
        let pending_action = match &node {
            StateNode::Local(_, local) => local.has_pending_action(),
            StateNode::Unit => false,
            StateNode::Scale(s) => s.operand.pending_action,
            StateNode::Sum(s) | StateNode::Tensor(s) => s.iter().any(|c| c.pending_action),
            StateNode::Act(_) => true,
        };
        Self {
            root: Some(Arc::new(node)),
            pending_action,
        }
    }
    pub(super) fn node(&self) -> &StateNode<I, L, C> {
        self.root.as_deref().expect("live state syntax")
    }
    pub(super) fn scaled(self, coefficient: C) -> Self {
        Self::new(StateNode::Scale(ScaleExpr {
            operand: self,
            coefficient,
        }))
    }
    pub(super) fn sum(children: SumExpr<Self>) -> Self {
        match children {
            SumExpr::One(child) => child,
            children => Self::new(StateNode::Sum(children)),
        }
    }
    pub(super) fn same_node(&self, rhs: &Self) -> bool {
        Arc::ptr_eq(
            self.root.as_ref().expect("live state"),
            rhs.root.as_ref().expect("live state"),
        )
    }
}

impl<I, L: LocalStructure, C> OpExpr<I, L, C> {
    pub(super) fn new(node: OpNode<I, L, C>) -> Self {
        Self {
            root: Some(Arc::new(node)),
        }
    }
    pub(super) fn node(&self) -> &OpNode<I, L, C> {
        self.root.as_deref().expect("live operator syntax")
    }
    fn scaled(self, coefficient: C) -> Self {
        Self::new(OpNode::Scale(ScaleExpr {
            operand: self,
            coefficient,
        }))
    }
}

impl<I, L: LocalStructure, C> State<I, L, C> {
    /// Places a whole local expression on a leaf, preserving its index.
    pub fn place(
        local: &LocalState<L, C>,
        factor: &HilbertSpace<I, L>,
    ) -> Result<Self, ExpressionError> {
        if factor.factor_index().is_none() {
            return Err(ExpressionError::NotAFactor);
        }
        Ok(Self {
            space: factor.clone(),
            expression: StateExpr::new(StateNode::Local(factor.clone(), local.clone())),
        })
    }
    /// Constructs a zero state with an explicit space, without numeric bounds.
    pub fn zero(space: &HilbertSpace<I, L>) -> Self {
        Self {
            space: space.clone(),
            expression: StateExpr::sum(SumExpr::Empty),
        }
    }
    /// Constructs the unit state on the empty tensor space, without generating 1.
    pub fn unit() -> Self {
        Self {
            space: HilbertSpace::unit(),
            expression: StateExpr::new(StateNode::Unit),
        }
    }
    /// Borrows the state's definite target space.
    pub fn space(&self) -> &HilbertSpace<I, L> {
        &self.space
    }
    /// Returns the scalar source space of this ket.
    pub fn domain(&self) -> HilbertSpace<I, L> {
        HilbertSpace::unit()
    }
    /// Borrows the ket's target space.
    pub fn codomain(&self) -> &HilbertSpace<I, L> {
        &self.space
    }
    /// Constructs a coefficient node sharing this state.
    pub fn scaled(&self, coefficient: C) -> Self {
        self.clone().into_scaled(coefficient)
    }
    pub(super) fn into_scaled(self, coefficient: C) -> Self {
        Self {
            space: self.space,
            expression: self.expression.scaled(coefficient),
        }
    }
}

impl<I: Eq, L: LocalStructure, C> State<I, L, C> {
    /// Constructs an ordered sum on an explicit space, including an empty sum.
    pub fn sum(
        space: &HilbertSpace<I, L>,
        states: impl IntoIterator<Item = Self>,
    ) -> Result<Self, ExpressionError> {
        let expressions = states
            .into_iter()
            .map(|state| {
                if !space.is_compatible_with(&state.space) {
                    return Err(ExpressionError::IncompatibleSpaces);
                }
                Ok(state.expression)
            })
            .collect::<Result<SumExpr<_>, _>>()?;
        Ok(Self {
            space: space.clone(),
            expression: StateExpr::sum(expressions),
        })
    }
    /// Adds states on compatible ordered spaces, without evaluating local rules.
    pub fn add(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        self.clone().into_add(rhs.clone())
    }
    pub(super) fn into_add(self, rhs: Self) -> Result<Self, ExpressionError> {
        if !self.space.is_compatible_with(&rhs.space) {
            return Err(ExpressionError::IncompatibleSpaces);
        }
        Ok(Self {
            space: self.space,
            expression: StateExpr::sum([self.expression, rhs.expression].into_iter().collect()),
        })
    }
    /// Forms a tensor product in input order. Empty input is the unit state.
    pub fn tensor(states: impl IntoIterator<Item = Self>) -> Result<Self, ExpressionError> {
        let states: SumExpr<Self> = states.into_iter().collect();
        let space = HilbertSpace::tensor(states.iter().map(|s| s.space.clone()))?;
        let parts: SumExpr<_> = states.into_iter().map(|s| s.expression).collect();
        let expression = match parts {
            SumExpr::Empty => StateExpr::new(StateNode::Unit),
            parts => StateExpr::new(StateNode::Tensor(parts)),
        };
        Ok(Self { space, expression })
    }
    /// Tensors this state with `rhs`, rejecting repeated factors.
    pub fn tensor_product(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        Self::tensor([self.clone(), rhs.clone()])
    }
    /// Constructs a delayed action using a binding or this state's space.
    pub fn act(&self, operator: &Operator<I, L, C>) -> Result<Self, ExpressionError> {
        self.act_selected(operator, None)
    }
    /// Constructs a delayed action using `ambient` in preference to an old binding.
    pub fn act_in(
        &self,
        operator: &Operator<I, L, C>,
        ambient: &HilbertSpace<I, L>,
    ) -> Result<Self, ExpressionError> {
        self.act_selected(operator, Some(ambient))
    }
    fn act_selected(
        &self,
        operator: &Operator<I, L, C>,
        ambient: Option<&HilbertSpace<I, L>>,
    ) -> Result<Self, ExpressionError> {
        let space = operator.select_space(ambient, Some(&self.space))?;
        if !space.is_compatible_with(&self.space) {
            return Err(ExpressionError::IncompatibleSpaces);
        }
        Ok(Self {
            space: space.clone(),
            expression: StateExpr::new(StateNode::Act(ActExpr {
                state: self.expression.clone(),
                operator: operator.expression.clone(),
            })),
        })
    }
    pub(super) fn into_act(self, operator: Operator<I, L, C>) -> Result<Self, ExpressionError> {
        let space = operator.select_space(None, Some(&self.space))?;
        if !space.is_compatible_with(&self.space) {
            return Err(ExpressionError::IncompatibleSpaces);
        }
        Ok(Self {
            space: space.clone(),
            expression: StateExpr::new(StateNode::Act(ActExpr {
                state: self.expression,
                operator: operator.expression,
            })),
        })
    }
}

impl<I, L: LocalStructure, C> Operator<I, L, C> {
    /// Places a whole local operator on a factor without binding an ambient space.
    pub fn place(
        local: &LocalOperator<L, C>,
        factor: &HilbertSpace<I, L>,
    ) -> Result<Self, ExpressionError> {
        if factor.factor_index().is_none() {
            return Err(ExpressionError::NotAFactor);
        }
        Ok(Self {
            ambient: None,
            expression: OpExpr::new(OpNode::Local(factor.clone(), local.clone())),
        })
    }
    /// Constructs an unbound zero (the empty sum).
    pub fn zero() -> Self {
        Self {
            ambient: None,
            expression: OpExpr::new(OpNode::Sum(SumExpr::Empty)),
        }
    }
    /// Constructs an unbound symbolic identity (the empty product).
    pub fn identity() -> Self {
        Self {
            ambient: None,
            expression: OpExpr::new(OpNode::Identity),
        }
    }
    /// Constructs a scalar multiple of identity without performing arithmetic.
    pub fn scalar(coefficient: C) -> Self {
        Self::identity().into_scaled(coefficient)
    }
    /// Borrows the ambient space, if bound.
    pub fn space(&self) -> Option<&HilbertSpace<I, L>> {
        self.ambient.as_ref()
    }
    /// Borrows the source space of a bound endomorphism.
    pub fn domain(&self) -> Option<&HilbertSpace<I, L>> {
        self.space()
    }
    /// Borrows the target space of a bound endomorphism.
    pub fn codomain(&self) -> Option<&HilbertSpace<I, L>> {
        self.space()
    }
    /// Constructs a coefficient node, preserving the optional binding.
    pub fn scaled(&self, coefficient: C) -> Self {
        self.clone().into_scaled(coefficient)
    }
    pub(super) fn into_scaled(self, coefficient: C) -> Self {
        Self {
            ambient: self.ambient,
            expression: self.expression.scaled(coefficient),
        }
    }
    /// Constructs a delayed adjoint, preserving the optional binding.
    pub fn adjoint(&self) -> Self {
        Self {
            ambient: self.ambient.clone(),
            expression: OpExpr::new(OpNode::Adjoint(AdjointExpr {
                operand: self.expression.clone(),
            })),
        }
    }
}

impl<I: Eq, L: LocalStructure, C> Operator<I, L, C> {
    fn combined_space(operands: &[Self]) -> Result<Option<HilbertSpace<I, L>>, ExpressionError> {
        // Decide whether every operand is bound before comparing any pair.
        if operands.is_empty() || operands.iter().any(|op| op.ambient.is_none()) {
            return Ok(None);
        }
        let first = operands[0].ambient.as_ref().expect("all operands bound");
        if operands
            .iter()
            .any(|op| !first.is_compatible_with(op.ambient.as_ref().expect("bound")))
        {
            return Err(ExpressionError::IncompatibleSpaces);
        }
        Ok(Some(first.clone()))
    }
    /// Constructs an ordered sum. A single unbound operand makes the result unbound.
    pub fn sum(operators: impl IntoIterator<Item = Self>) -> Result<Self, ExpressionError> {
        let operators: SumExpr<_> = operators.into_iter().collect();
        let ambient = Self::combined_space(operators.as_slice())?;
        let parts: SumExpr<_> = operators.into_iter().map(|o| o.expression).collect();
        let expression = match parts {
            SumExpr::One(op) => op,
            parts => OpExpr::new(OpNode::Sum(parts)),
        };
        Ok(Self {
            expression,
            ambient,
        })
    }
    /// Adds operators, validating compatible bindings when both are bound.
    pub fn add(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        Self::sum([self.clone(), rhs.clone()])
    }
    pub(super) fn into_add(self, rhs: Self) -> Result<Self, ExpressionError> {
        Self::sum([self, rhs])
    }
    /// Constructs a product in input order; the empty product is unbound identity.
    pub fn product_all(operators: impl IntoIterator<Item = Self>) -> Result<Self, ExpressionError> {
        let operators: SumExpr<_> = operators.into_iter().collect();
        let ambient = Self::combined_space(operators.as_slice())?;
        let mut parts = operators.into_iter().map(|o| o.expression);
        let expression = match parts.next() {
            None => OpExpr::new(OpNode::Identity),
            Some(first) => parts.fold(first, |lhs, rhs| {
                OpExpr::new(OpNode::Product(ProductExpr { lhs, rhs }))
            }),
        };
        Ok(Self {
            expression,
            ambient,
        })
    }
    /// Constructs `self * rhs`, with the right operand acting first.
    pub fn product(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        Self::product_all([self.clone(), rhs.clone()])
    }
    pub(super) fn into_product(self, rhs: Self) -> Result<Self, ExpressionError> {
        Self::product_all([self, rhs])
    }
    /// Constructs the anticommutator without evaluating local rules.
    pub fn anticommutator(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        Self::sum([self.product(rhs)?, rhs.product(self)?])
    }
    /// Binds or rebinds this syntax to an explicit space, leaving the input intact.
    pub fn bind(&self, ambient: &HilbertSpace<I, L>) -> Result<Self, ExpressionError> {
        self.validate(ambient)?;
        Ok(Self {
            expression: self.expression.clone(),
            ambient: Some(ambient.clone()),
        })
    }
    pub(super) fn select_space<'a>(
        &'a self,
        explicit: Option<&'a HilbertSpace<I, L>>,
        state: Option<&'a HilbertSpace<I, L>>,
    ) -> Result<&'a HilbertSpace<I, L>, ExpressionError> {
        let space = explicit
            .or(self.ambient.as_ref())
            .or(state)
            .ok_or(ExpressionError::UnboundOperator)?;
        self.validate(space)?;
        Ok(space)
    }
    fn validate(&self, space: &HilbertSpace<I, L>) -> Result<(), ExpressionError> {
        let mut pending = alloc::vec![&self.expression];
        while let Some(op) = pending.pop() {
            match op.node() {
                OpNode::Local(factor, _) => {
                    if !space.contains_index(factor.factor_index().expect("placed factor")) {
                        return Err(ExpressionError::MissingIndex);
                    }
                }
                OpNode::Identity => {}
                OpNode::Scale(s) => pending.push(&s.operand),
                OpNode::Sum(s) => pending.extend(s.iter()),
                OpNode::Product(p) => {
                    pending.push(&p.lhs);
                    pending.push(&p.rhs);
                }
                OpNode::Adjoint(a) => pending.push(&a.operand),
            }
        }
        Ok(())
    }
}

impl<L: LocalStructure, C> LocalState<L, C> {
    /// Places this whole expression on an existing factor, without evaluating it.
    pub fn at<I>(&self, factor: &HilbertSpace<I, L>) -> Result<State<I, L, C>, ExpressionError> {
        State::place(self, factor)
    }
}
impl<L: LocalStructure, C> LocalOperator<L, C> {
    /// Places this whole expression on an existing factor, without binding it.
    pub fn at<I>(&self, factor: &HilbertSpace<I, L>) -> Result<Operator<I, L, C>, ExpressionError> {
        Operator::place(self, factor)
    }
}

impl<I, L: LocalStructure, C> Drop for StateExpr<I, L, C> {
    fn drop(&mut self) {
        let mut first = self.root.take();
        let mut pending = Vec::new();
        while let Some(root) = first.take().or_else(|| pending.pop()) {
            if let Some(node) = Arc::into_inner(root) {
                let mut push = |mut child: Self| {
                    if let Some(root) = child.root.take() {
                        pending.push(root);
                    }
                };
                match node {
                    StateNode::Local(_, _) | StateNode::Unit => {}
                    StateNode::Scale(s) => push(s.operand),
                    StateNode::Sum(s) | StateNode::Tensor(s) => {
                        for child in s {
                            push(child);
                        }
                    }
                    StateNode::Act(a) => push(a.state),
                }
            }
        }
    }
}
impl<I, L: LocalStructure, C> Drop for OpExpr<I, L, C> {
    fn drop(&mut self) {
        let mut first = self.root.take();
        let mut pending = Vec::new();
        while let Some(root) = first.take().or_else(|| pending.pop()) {
            if let Some(node) = Arc::into_inner(root) {
                let mut push = |mut child: Self| {
                    if let Some(root) = child.root.take() {
                        pending.push(root);
                    }
                };
                match node {
                    OpNode::Local(_, _) | OpNode::Identity => {}
                    OpNode::Scale(s) => push(s.operand),
                    OpNode::Sum(s) => {
                        for child in s {
                            push(child);
                        }
                    }
                    OpNode::Product(p) => {
                        push(p.lhs);
                        push(p.rhs);
                    }
                    OpNode::Adjoint(a) => push(a.operand),
                }
            }
        }
    }
}
