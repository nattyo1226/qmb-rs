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
pub trait Product<T>: Sized {
    /// Returns `self * rhs`; `rhs` acts first when applied to a state.
    fn product(&self, rhs: &Self) -> LinearCombination<Self, T>;
}

/// The action of an operator primitive on a compatible state primitive.
///
/// Implementing one action does not enable action on arbitrary state types:
///
/// ```compile_fail
/// use qmb_core::{ActOn, LinearCombination};
/// struct Operator;
/// struct State;
/// struct OtherState;
/// impl ActOn<State, f64> for Operator {
///     fn act_on(&self, _: &State) -> LinearCombination<State, f64> {
///         LinearCombination::primitive(State)
///     }
/// }
/// let _ = Operator.act_on(&OtherState);
/// ```
pub trait ActOn<S, T> {
    /// Returns the result of acting on `state`.
    fn act_on(&self, state: &S) -> LinearCombination<S, T>;
}

/// The adjoint of an operator primitive.
pub trait Adjoint<T>: Sized {
    /// Returns the adjoint as a finite linear combination.
    fn adjoint(&self) -> LinearCombination<Self, T>;
}

/// The ordinary, unnormalized trace of an operator primitive.
pub trait Trace<T> {
    /// Returns the trace without dividing by the local dimension.
    fn trace(&self) -> Complex<T>;
}
