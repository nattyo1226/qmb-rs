//! A choice of primitives and coefficient type, independent of space and geometry.

use crate::{LocalOperator, LocalState, Operator, State};

/// Selects primitive types and coefficients for constructing ordinary expressions.
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
/// ```
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
}

/// A local state with the primitives and coefficients chosen by `Q`.
pub type LocalStateOf<Q> = LocalState<
    <Q as QuantumSystem>::State,
    <Q as QuantumSystem>::Operator,
    <Q as QuantumSystem>::Coeff,
>;

/// A local operator with the primitive and coefficients chosen by `Q`.
pub type LocalOperatorOf<Q> =
    LocalOperator<<Q as QuantumSystem>::Operator, <Q as QuantumSystem>::Coeff>;

/// A composite state with system `Q` and independent index type `I`.
pub type StateOf<Q, I> = State<
    I,
    <Q as QuantumSystem>::State,
    <Q as QuantumSystem>::Operator,
    <Q as QuantumSystem>::Coeff,
>;

/// A composite operator with system `Q` and independent index type `I`.
pub type OperatorOf<Q, I> =
    Operator<I, <Q as QuantumSystem>::Operator, <Q as QuantumSystem>::Coeff>;
