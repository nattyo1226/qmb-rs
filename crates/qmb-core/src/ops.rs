//! Primitive capabilities implemented by concrete local algebras.
//!
//! Rules are infallible for valid primitives. Concrete types must reject invalid
//! values at construction. Numeric bounds and algebraic laws are the concrete
//! implementation's responsibility, rather than blanket bounds on these traits.

use crate::{Complex, LinearCombination};

/// The inner product, conjugate-linear on the left and linear on the right.
pub trait InnerProduct<T> {
    /// Returns the inner product of `self` with `rhs`.
    fn inner_product(&self, rhs: &Self) -> Complex<T>;
}

/// Ordered multiplication of two operator primitives.
pub trait Product<T> {
    /// Element type in the product rule result.
    type Term;
    /// Returns `self * rhs`; `rhs` acts first when applied to a state.
    fn product(&self, rhs: &Self) -> LinearCombination<Self::Term, T>;
}

/// A state primitive that can receive the left action of an operator primitive.
///
/// Implementing one action does not enable action by arbitrary operator types:
///
/// ```compile_fail
/// use qmb_core::{Act, LinearCombination};
/// struct Operator;
/// struct State;
/// struct OtherOperator;
/// impl Act<f64> for State {
///     type Operator = Operator;
///     type Term = Self;
///     fn act(&self, _: &Operator) -> LinearCombination<Self, f64> {
///         LinearCombination::primitive(State)
///     }
/// }
/// let _ = State.act(&OtherOperator);
/// ```
pub trait Act<T> {
    /// Operator type accepted by this action.
    type Operator;
    /// Element type in the action rule result.
    type Term;
    /// Returns `operator * self` as a finite linear combination of states.
    fn act(&self, operator: &Self::Operator) -> LinearCombination<Self::Term, T>;
}

/// The adjoint of an operator primitive.
pub trait Adjoint<T> {
    /// Element type in the adjoint rule result.
    type Term;
    /// Returns the adjoint as a finite linear combination.
    fn adjoint(&self) -> LinearCombination<Self::Term, T>;
}

/// The ordinary, unnormalized trace of an operator primitive.
pub trait Trace<T> {
    /// Returns the trace without dividing by the local dimension.
    fn trace(&self) -> Complex<T>;
}
