//! A choice of primitives and coefficient type, independent of space and geometry.

use crate::{ExpressionError, HilbertSpace, LocalOperator, LocalState, Operator, State};

/// Selects primitive types and coefficients for constructing local and composite expressions.
///
/// No evaluation, numeric or cloning capabilities are required. Those bounds
/// belong to the operations on the returned expressions. Systems do not own
/// spaces or indices, and do not become part of expression identity: two systems
/// with the same associated types construct the same expression types.
///
/// ```
/// use qmb_core::{QuantumSystem, LocalStateOf};
/// // A custom algebra can opt into construction before providing evaluation.
/// struct Ket;
/// struct Op;
/// struct Coeff;
/// struct Custom;
/// impl QuantumSystem for Custom {
///     type State = Ket;
///     type Operator = Op;
///     type Coeff = Coeff;
/// }
/// let ket: LocalStateOf<Custom> = Custom::state(Ket);
/// let delayed = ket.act(&Custom::operator(Op));
/// let a = qmb_core::HilbertSpace::factor("a");
/// let b = qmb_core::HilbertSpace::factor("b");
/// let state = Custom::state_tensor_at([(&a, Ket), (&b, Ket)])?;
/// let operator = Custom::try_operator_sum([Custom::operator_product_at([(&a, Op)])])?;
/// let delayed = state.act(&operator)?;
/// # Ok::<(), qmb_core::ExpressionError>(())
/// ```
/// Construction and collection live on this trait; operations on an existing
/// expression (scaling, placement, binding and evaluation) remain methods of the
/// expression. Type aliases such as [`StateOf`] are primarily for annotations.
///
/// All `try_` collections accept `ExpressionError`, stop at the first input
/// error, and validate composition only after all inputs succeed. An input
/// error therefore precedes any space errors in earlier successful operands.
/// Only iterator consumption stops: work already done to build an array or
/// vector is not rolled back.
pub trait QuantumSystem: Sized {
    /// State primitive, rather than a state expression.
    type State;
    /// Operator primitive, rather than an operator expression.
    type Operator;
    /// The whole coefficient type, including any complex components.
    type Coeff;

    /// Holds a state primitive without evaluating it.
    fn state(primitive: Self::State) -> LocalStateOf<Self> {
        LocalState::primitive(primitive)
    }

    /// Holds an operator primitive without evaluating it.
    fn operator(primitive: Self::Operator) -> LocalOperatorOf<Self> {
        LocalOperator::primitive(primitive)
    }

    /// Places primitives on leaves and tensors them in input order.
    ///
    /// Repeated indices are rejected; empty input gives the unit state.
    /// The first placement error takes precedence over duplicate indices.
    /// Each target must be a leaf; tensors are not expanded automatically.
    fn state_tensor_at<'a, I: Eq + 'a>(
        assignments: impl IntoIterator<Item = (&'a HilbertSpace<I>, Self::State)>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        Self::state_tensor_at_with_kind(crate::TensorKind::Ordinary, assignments)
    }

    /// Places primitives and tensors them using an explicit rule.
    fn state_tensor_at_with_kind<'a, I: Eq + 'a>(
        kind: crate::TensorKind,
        assignments: impl IntoIterator<Item = (&'a HilbertSpace<I>, Self::State)>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        Self::try_state_tensor_with_kind(
            kind,
            assignments
                .into_iter()
                .map(|(factor, primitive)| Self::state(primitive).at(factor)),
        )
    }
    /// Places primitives and forms a graded tensor.
    fn state_tensor_g_at<'a, I: Eq + 'a>(
        assignments: impl IntoIterator<Item = (&'a HilbertSpace<I>, Self::State)>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        Self::state_tensor_at_with_kind(crate::TensorKind::Graded, assignments)
    }
    /// Forms a tensor state using the specified rule.
    fn state_tensor_with_kind<I: Eq>(
        kind: crate::TensorKind,
        states: impl IntoIterator<Item = StateOf<Self, I>>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        State::tensor_with_kind(kind, states)
    }
    /// Forms a graded tensor state.
    fn state_tensor_g<I: Eq>(
        states: impl IntoIterator<Item = StateOf<Self, I>>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        Self::state_tensor_with_kind(crate::TensorKind::Graded, states)
    }
    /// Collects fallible tensor operands before validating composition.
    fn try_state_tensor_with_kind<I: Eq>(
        kind: crate::TensorKind,
        states: impl IntoIterator<Item = Result<StateOf<Self, I>, ExpressionError>>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        State::try_tensor_with_kind(kind, states)
    }
    /// Collects fallible graded tensor operands, stopping at the first error.
    fn try_state_tensor_g<I: Eq>(
        states: impl IntoIterator<Item = Result<StateOf<Self, I>, ExpressionError>>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        Self::try_state_tensor_with_kind(crate::TensorKind::Graded, states)
    }

    /// Places primitives on leaves and multiplies them in input order.
    ///
    /// Repeated indices are allowed. The rightmost operator acts first.
    /// The result is unbound, including identity for empty input; use `bind`
    /// to select an ambient space. The first placement error stops iteration.
    fn operator_product_at<'a, I: Eq + 'a>(
        assignments: impl IntoIterator<Item = (&'a HilbertSpace<I>, Self::Operator)>,
    ) -> Result<OperatorOf<Self, I>, ExpressionError> {
        Self::try_operator_product(
            assignments
                .into_iter()
                .map(|(factor, primitive)| Self::operator(primitive).at(factor)),
        )
    }

    /// Sums placed states on `space`, including zero for empty input.
    fn state_sum<I: Eq>(
        space: &HilbertSpace<I>,
        states: impl IntoIterator<Item = StateOf<Self, I>>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        State::sum(space, states)
    }

    /// Sums fallible states; the first input error precedes space validation.
    /// Later inputs are not pulled after an error. Empty input is zero on `space`.
    fn try_state_sum<I: Eq>(
        space: &HilbertSpace<I>,
        states: impl IntoIterator<Item = Result<StateOf<Self, I>, ExpressionError>>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        State::try_sum(space, states)
    }

    /// Tensors placed states in input order, rejecting repeated indices.
    /// Empty input gives the unit state.
    fn state_tensor<I: Eq>(
        states: impl IntoIterator<Item = StateOf<Self, I>>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        State::tensor(states)
    }

    /// Collects fallible states before tensoring, with input errors taking priority.
    /// The first input error stops iteration. Empty input gives the unit state.
    fn try_state_tensor<I: Eq>(
        states: impl IntoIterator<Item = Result<StateOf<Self, I>, ExpressionError>>,
    ) -> Result<StateOf<Self, I>, ExpressionError> {
        State::try_tensor(states)
    }

    /// Sums operators, validating bindings only if every operand is bound.
    /// Empty input gives unbound zero.
    fn operator_sum<I: Eq>(
        operators: impl IntoIterator<Item = OperatorOf<Self, I>>,
    ) -> Result<OperatorOf<Self, I>, ExpressionError> {
        Operator::sum(operators)
    }

    /// Sums fallible operators, returning the first input error before validation.
    /// Later inputs are not pulled after an error. Empty input is unbound zero.
    fn try_operator_sum<I: Eq>(
        operators: impl IntoIterator<Item = Result<OperatorOf<Self, I>, ExpressionError>>,
    ) -> Result<OperatorOf<Self, I>, ExpressionError> {
        Operator::try_sum(operators)
    }

    /// Multiplies placed operators in input order; the rightmost acts first.
    /// Bindings follow the same rules as sums; empty input is unbound identity.
    fn operator_product<I: Eq>(
        operators: impl IntoIterator<Item = OperatorOf<Self, I>>,
    ) -> Result<OperatorOf<Self, I>, ExpressionError> {
        Operator::product_all(operators)
    }

    /// Multiplies fallible operators, with input errors preceding binding errors.
    /// The first input error stops iteration. Empty input is unbound identity.
    fn try_operator_product<I: Eq>(
        operators: impl IntoIterator<Item = Result<OperatorOf<Self, I>, ExpressionError>>,
    ) -> Result<OperatorOf<Self, I>, ExpressionError> {
        Operator::try_product_all(operators)
    }
}

/// A local state type annotation for `Q`; construct it with [`QuantumSystem::state`].
pub type LocalStateOf<Q> = LocalState<
    <Q as QuantumSystem>::State,
    <Q as QuantumSystem>::Operator,
    <Q as QuantumSystem>::Coeff,
>;

/// A local operator type annotation for `Q`; use [`QuantumSystem::operator`] to construct it.
pub type LocalOperatorOf<Q> =
    LocalOperator<<Q as QuantumSystem>::Operator, <Q as QuantumSystem>::Coeff>;

/// A composite state type annotation for `Q` and independent index `I`.
/// Construct it with methods such as [`QuantumSystem::state_tensor_at`].
pub type StateOf<Q, I> = State<
    I,
    <Q as QuantumSystem>::State,
    <Q as QuantumSystem>::Operator,
    <Q as QuantumSystem>::Coeff,
>;

/// A composite operator type annotation for `Q` and independent index `I`.
/// Construct it with methods such as [`QuantumSystem::operator_product_at`].
pub type OperatorOf<Q, I> =
    Operator<I, <Q as QuantumSystem>::Operator, <Q as QuantumSystem>::Coeff>;
