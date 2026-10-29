//! Placed states and operators on ordinary and graded tensor spaces.

use alloc::{sync::Arc, vec::Vec};

use super::parts::{AdjointExpr, ProductExpr, ScaleExpr, SumExpr};
use crate::{HilbertSpace, Identity, LocalOperator, LocalState, SpaceError};

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

/// A composite state with primitive types S and O and an immutable factor space.
///
/// Tensor products require the same primitive types even though their spaces
/// carry only index structure:
///
/// ```compile_fail
/// use qmb_core::{HilbertSpace, LocalState};
/// let a = HilbertSpace::factor(0u8);
/// let b = HilbertSpace::factor(1u8);
/// let left = LocalState::<u8, (), ()>::primitive(0).at(&a).unwrap();
/// let right = LocalState::<u16, (), ()>::primitive(0).at(&b).unwrap();
/// let h = a.tensor(&b).unwrap();
/// let _ = qmb_core::State::tensor_in(&h, [left, right]);
/// ```
///
/// Composite action also requires matching operator primitives:
///
/// ```compile_fail
/// use qmb_core::{HilbertSpace, LocalState, Operator};
/// let h = HilbertSpace::factor(0u8);
/// let state = LocalState::<(), u8, ()>::zero().at(&h).unwrap();
/// let operator = Operator::<u8, u16, ()>::identity().bind(&h).unwrap();
/// let _ = state.act(&operator);
/// ```
///
/// Coefficient precisions cannot be mixed:
///
/// ```compile_fail
/// use qmb_core::{HilbertSpace, State};
/// let h = HilbertSpace::factor(0u8);
/// let a = State::<_, (), (), f32>::zero(&h);
/// let b = State::<_, (), (), f64>::zero(&h);
/// let _ = a + b;
/// ```
///
/// Shared syntax still follows the index's Send bound:
///
/// ```compile_fail
/// use std::rc::Rc;
/// use qmb_core::State;
/// fn needs_send<T: Send>() {}
/// needs_send::<State<Rc<()>, (), (), ()>>();
/// ```
pub struct State<I, S, O, C> {
    pub(crate) expression: StateExpr<I, S, O, C>,
    pub(crate) space: HilbertSpace<I>,
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
pub struct Operator<I, O, C> {
    pub(crate) expression: OpExpr<I, O, C>,
    pub(crate) ambient: Option<HilbertSpace<I>>,
}

pub(crate) struct StateExpr<I, S, O, C> {
    root: Option<Arc<StateNode<I, S, O, C>>>,
    pub(crate) pending_action: bool,
}
pub(crate) struct OpExpr<I, O, C> {
    root: Option<Arc<OpNode<I, O, C>>>,
}

pub(crate) struct StateAction<I, S, O, C> {
    pub(crate) state: StateExpr<I, S, O, C>,
    pub(crate) operator: OpExpr<I, O, C>,
    pub(crate) space: HilbertSpace<I>,
}

pub(crate) enum StateNode<I, S, O, C> {
    Local(HilbertSpace<I>, LocalState<S, O, C>),
    Unit,
    Scale(ScaleExpr<StateExpr<I, S, O, C>, C>),
    Sum(SumExpr<StateExpr<I, S, O, C>>),
    Tensor(SumExpr<StateExpr<I, S, O, C>>),
    Act(StateAction<I, S, O, C>),
}
pub(crate) enum OpNode<I, O, C> {
    Local(HilbertSpace<I>, LocalOperator<O, C>),
    Identity,
    Scale(ScaleExpr<OpExpr<I, O, C>, C>),
    Sum(SumExpr<OpExpr<I, O, C>>),
    Product(ProductExpr<OpExpr<I, O, C>>),
    Adjoint(AdjointExpr<OpExpr<I, O, C>>),
}

impl<I, S, O, C> Clone for StateExpr<I, S, O, C> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
            pending_action: self.pending_action,
        }
    }
}
impl<I, O, C> Clone for OpExpr<I, O, C> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
        }
    }
}
impl<I, S, O, C> Clone for State<I, S, O, C> {
    fn clone(&self) -> Self {
        Self {
            expression: self.expression.clone(),
            space: self.space.clone(),
        }
    }
}
impl<I, O, C> Clone for Operator<I, O, C> {
    fn clone(&self) -> Self {
        Self {
            expression: self.expression.clone(),
            ambient: self.ambient.clone(),
        }
    }
}

impl<I, S, O, C> StateExpr<I, S, O, C> {
    pub(crate) fn new(node: StateNode<I, S, O, C>) -> Self {
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
    pub(crate) fn node(&self) -> &StateNode<I, S, O, C> {
        self.root.as_deref().expect("live state syntax")
    }
    pub(crate) fn scaled(self, coefficient: C) -> Self {
        Self::new(StateNode::Scale(ScaleExpr {
            operand: self,
            coefficient,
        }))
    }
    pub(crate) fn sum(children: SumExpr<Self>) -> Self {
        match children {
            SumExpr::One(child) => child,
            children => Self::new(StateNode::Sum(children)),
        }
    }
    pub(crate) fn shares_node_with(&self, rhs: &Self) -> bool {
        Arc::ptr_eq(
            self.root.as_ref().expect("live state"),
            rhs.root.as_ref().expect("live state"),
        )
    }
}

impl<I, O, C> OpExpr<I, O, C> {
    pub(crate) fn new(node: OpNode<I, O, C>) -> Self {
        Self {
            root: Some(Arc::new(node)),
        }
    }
    pub(crate) fn node(&self) -> &OpNode<I, O, C> {
        self.root.as_deref().expect("live operator syntax")
    }
    fn scaled(self, coefficient: C) -> Self {
        Self::new(OpNode::Scale(ScaleExpr {
            operand: self,
            coefficient,
        }))
    }
}

impl<I, S, O, C> State<I, S, O, C> {
    /// Places a whole local expression on a leaf, preserving its index.
    pub fn place(
        local: &LocalState<S, O, C>,
        factor: &HilbertSpace<I>,
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
    pub fn zero(space: &HilbertSpace<I>) -> Self {
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
    pub fn space(&self) -> &HilbertSpace<I> {
        &self.space
    }
    /// Returns the scalar source space of this ket.
    pub fn domain(&self) -> HilbertSpace<I> {
        HilbertSpace::unit()
    }
    /// Borrows the ket's target space.
    pub fn codomain(&self) -> &HilbertSpace<I> {
        &self.space
    }
    /// Constructs a coefficient node sharing this state.
    /// Converts the input to C exactly once without changing coefficient precision.
    pub fn scaled(&self, coefficient: impl Into<C>) -> Self {
        self.clone().into_scaled(coefficient.into())
    }
    pub(crate) fn into_scaled(self, coefficient: C) -> Self {
        Self {
            space: self.space,
            expression: self.expression.scaled(coefficient),
        }
    }
}

impl<I: Eq, S, O, C> State<I, S, O, C> {
    /// Constructs an ordered sum on an explicit space, including an empty sum.
    pub fn sum(
        space: &HilbertSpace<I>,
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
    /// Collects fallible states before checking compatibility with `space`.
    ///
    /// The first input error stops iteration and takes precedence over space
    /// errors in earlier successful operands. Empty input is zero on `space`.
    /// Errors do not roll back work already performed by the input producer.
    pub fn try_sum(
        space: &HilbertSpace<I>,
        states: impl IntoIterator<Item = Result<Self, ExpressionError>>,
    ) -> Result<Self, ExpressionError> {
        let states: SumExpr<Self> = states.into_iter().collect::<Result<_, _>>()?;
        Self::sum(space, states)
    }
    pub(crate) fn into_add(self, rhs: Self) -> Result<Self, ExpressionError> {
        if !self.space.is_compatible_with(&rhs.space) {
            return Err(ExpressionError::IncompatibleSpaces);
        }
        Ok(Self {
            space: self.space,
            expression: StateExpr::sum([self.expression, rhs.expression].into_iter().collect()),
        })
    }
    /// Tensors placed states on an explicit space, preserving its shared handle.
    ///
    /// Ordered input indices must cover the space exactly. Multi-factor inputs
    /// must have a compatible tensor kind; grouping may differ. Empty input is
    /// accepted only on a unit space. No primitive evaluation is performed.
    pub fn tensor_in(
        space: &HilbertSpace<I>,
        states: impl IntoIterator<Item = Self>,
    ) -> Result<Self, ExpressionError> {
        Self::tensor_collected(space, states.into_iter().collect())
    }

    fn tensor_collected(
        space: &HilbertSpace<I>,
        states: SumExpr<Self>,
    ) -> Result<Self, ExpressionError> {
        if states.iter().any(|s| {
            s.space
                .effective_kind()
                .is_some_and(|kind| Some(kind) != space.effective_kind())
        }) || !states
            .iter()
            .flat_map(|s| s.space.indices())
            .eq(space.indices())
        {
            return Err(ExpressionError::IncompatibleSpaces);
        }
        let parts: SumExpr<_> = states.into_iter().map(|s| s.expression).collect();
        let expression = match parts {
            SumExpr::Empty => StateExpr::new(StateNode::Unit),
            SumExpr::One(node) => node,
            parts => StateExpr::new(StateNode::Tensor(parts)),
        };
        Ok(Self {
            space: space.clone(),
            expression,
        })
    }

    /// Collects fallible inputs before validating coverage of the explicit space.
    /// The first input error stops consumption and precedes all space errors.
    pub fn try_tensor_in(
        space: &HilbertSpace<I>,
        states: impl IntoIterator<Item = Result<Self, ExpressionError>>,
    ) -> Result<Self, ExpressionError> {
        let states: SumExpr<Self> = states.into_iter().collect::<Result<_, _>>()?;
        Self::tensor_collected(space, states)
    }
    /// Constructs a delayed action using a binding or this state's space.
    pub fn act(&self, operator: &Operator<I, O, C>) -> Result<Self, ExpressionError> {
        self.act_selected(operator, None)
    }
    /// Constructs a delayed action using `ambient` in preference to an old binding.
    pub fn act_in(
        &self,
        operator: &Operator<I, O, C>,
        ambient: &HilbertSpace<I>,
    ) -> Result<Self, ExpressionError> {
        self.act_selected(operator, Some(ambient))
    }
    fn act_selected(
        &self,
        operator: &Operator<I, O, C>,
        ambient: Option<&HilbertSpace<I>>,
    ) -> Result<Self, ExpressionError> {
        let space = operator.select_space(ambient, Some(&self.space))?;
        if !space.is_compatible_with(&self.space) {
            return Err(ExpressionError::IncompatibleSpaces);
        }
        Ok(Self {
            space: space.clone(),
            expression: StateExpr::new(StateNode::Act(StateAction {
                space: space.clone(),
                state: self.expression.clone(),
                operator: operator.expression.clone(),
            })),
        })
    }
    pub(crate) fn into_act(self, operator: Operator<I, O, C>) -> Result<Self, ExpressionError> {
        let space = operator.select_space(None, Some(&self.space))?;
        if !space.is_compatible_with(&self.space) {
            return Err(ExpressionError::IncompatibleSpaces);
        }
        Ok(Self {
            space: space.clone(),
            expression: StateExpr::new(StateNode::Act(StateAction {
                space: space.clone(),
                state: self.expression,
                operator: operator.expression,
            })),
        })
    }
}

impl<I, O, C> Operator<I, O, C> {
    /// Places a whole local operator on a factor without binding an ambient space.
    pub fn place(
        local: &LocalOperator<O, C>,
        factor: &HilbertSpace<I>,
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
    pub fn space(&self) -> Option<&HilbertSpace<I>> {
        self.ambient.as_ref()
    }
    /// Borrows the source space of a bound endomorphism.
    pub fn domain(&self) -> Option<&HilbertSpace<I>> {
        self.space()
    }
    /// Borrows the target space of a bound endomorphism.
    pub fn codomain(&self) -> Option<&HilbertSpace<I>> {
        self.space()
    }
    /// Constructs a coefficient node, preserving the optional binding.
    /// Converts the input to C exactly once without changing coefficient precision.
    pub fn scaled(&self, coefficient: impl Into<C>) -> Self {
        self.clone().into_scaled(coefficient.into())
    }
    pub(crate) fn into_scaled(self, coefficient: C) -> Self {
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

impl<I: Eq, O, C> Operator<I, O, C> {
    fn combined_space(operands: &[Self]) -> Result<Option<HilbertSpace<I>>, ExpressionError> {
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
        Self::sum_collected(operators.into_iter().collect())
    }
    /// Collects fallible operators before validating the whole sum's binding.
    ///
    /// Stops at the first input error and returns it unchanged, ahead of any
    /// binding errors in successful operands. Empty input is unbound zero.
    /// An unbound operand makes the entire sum unbound, just as in [`Self::sum`].
    pub fn try_sum(
        operators: impl IntoIterator<Item = Result<Self, ExpressionError>>,
    ) -> Result<Self, ExpressionError> {
        Self::sum_collected(operators.into_iter().collect::<Result<_, _>>()?)
    }
    fn sum_collected(operators: SumExpr<Self>) -> Result<Self, ExpressionError> {
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
    pub(crate) fn into_add(self, rhs: Self) -> Result<Self, ExpressionError> {
        Self::sum([self, rhs])
    }
    /// Constructs a product in input order; the empty product is unbound identity.
    pub fn product(operators: impl IntoIterator<Item = Self>) -> Result<Self, ExpressionError> {
        Self::product_collected(operators.into_iter().collect())
    }
    /// Collects fallible operators before validating the whole product's binding.
    ///
    /// The first input error stops iteration and takes precedence over binding
    /// errors. Order is preserved, with the right operand acting first.
    /// Empty input is unbound identity. No primitive rules are evaluated.
    pub fn try_product(
        operators: impl IntoIterator<Item = Result<Self, ExpressionError>>,
    ) -> Result<Self, ExpressionError> {
        Self::product_collected(operators.into_iter().collect::<Result<_, _>>()?)
    }
    fn product_collected(operators: SumExpr<Self>) -> Result<Self, ExpressionError> {
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
    pub(crate) fn into_product(self, rhs: Self) -> Result<Self, ExpressionError> {
        Self::product([self, rhs])
    }
    /// Constructs the anticommutator without evaluating local rules.
    pub fn anticommutator(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        Self::sum([(self * rhs)?, (rhs * self)?])
    }
    /// Binds or rebinds this syntax to an explicit space, leaving the input intact.
    pub fn bind(&self, ambient: &HilbertSpace<I>) -> Result<Self, ExpressionError> {
        self.validate(ambient)?;
        Ok(Self {
            expression: self.expression.clone(),
            ambient: Some(ambient.clone()),
        })
    }
    pub(crate) fn select_space<'a>(
        &'a self,
        explicit: Option<&'a HilbertSpace<I>>,
        state: Option<&'a HilbertSpace<I>>,
    ) -> Result<&'a HilbertSpace<I>, ExpressionError> {
        let space = explicit
            .or(self.ambient.as_ref())
            .or(state)
            .ok_or(ExpressionError::UnboundOperator)?;
        self.validate(space)?;
        Ok(space)
    }
    fn validate(&self, space: &HilbertSpace<I>) -> Result<(), ExpressionError> {
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

impl<S, O, C> LocalState<S, O, C> {
    /// Places this whole expression on an existing factor, without evaluating it.
    pub fn at<I>(&self, factor: &HilbertSpace<I>) -> Result<State<I, S, O, C>, ExpressionError> {
        State::place(self, factor)
    }
}
impl<O, C> LocalOperator<O, C> {
    /// Places this whole expression on an existing factor, without binding it.
    pub fn at<I>(&self, factor: &HilbertSpace<I>) -> Result<Operator<I, O, C>, ExpressionError> {
        Operator::place(self, factor)
    }
}

impl<I, S, O, C> Drop for StateExpr<I, S, O, C> {
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
impl<I, O, C> Drop for OpExpr<I, O, C> {
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

impl<I, O, C> Identity for Operator<I, O, C> {
    fn identity() -> Self {
        Self::identity()
    }
}
