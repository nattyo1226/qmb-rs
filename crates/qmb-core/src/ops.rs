//! Dimension queries and algebraic capabilities of primitives and expressions.
//!
//! Rules are infallible for valid primitives. Concrete types must reject invalid
//! values at construction. Numeric bounds and algebraic laws are the concrete
//! implementation's responsibility, rather than blanket bounds on these traits.

use crate::{DimensionOverflow, LinearCombination};
use core::num::NonZeroUsize;

/// The inner product, conjugate-linear on the left and linear on the right.
pub trait InnerProduct<C> {
    /// Returns the inner product of `self` with `rhs`.
    fn inner_product(&self, rhs: &Self) -> C;
}

/// Ordered multiplication of two operator primitives.
pub trait Product<C> {
    /// Element type in the product rule result.
    type Term;
    /// Returns `self * rhs`; `rhs` acts first when applied to a state.
    fn product(&self, rhs: &Self) -> LinearCombination<Self::Term, C>;
}

/// A state primitive that can receive the left action of an operator primitive.
///
/// Implementing one action does not enable action by arbitrary operator types:
///
/// ```compile_fail
/// use num_complex::Complex;
/// use qmb_core::{Act, LinearCombination};
/// struct Operator;
/// struct State;
/// struct OtherOperator;
/// impl Act<Complex<f64>> for State {
///     type Operator = Operator;
///     type Term = Self;
///     fn act(&self, _: &Operator) -> LinearCombination<Self, Complex<f64>> {
///         LinearCombination::primitive(State)
///     }
/// }
/// let _ = State.act(&OtherOperator);
/// ```
pub trait Act<C> {
    /// Operator type accepted by this action.
    type Operator;
    /// Element type in the action rule result.
    type Term;
    /// Returns `operator * self` as a finite linear combination of states.
    fn act(&self, operator: &Self::Operator) -> LinearCombination<Self::Term, C>;
}

/// The adjoint of an operator primitive.
pub trait Adjoint<C> {
    /// Element type in the adjoint rule result.
    type Term;
    /// Returns the adjoint as a finite linear combination.
    fn adjoint(&self) -> LinearCombination<Self::Term, C>;
}

/// The ordinary, unnormalized trace of an operator primitive.
pub trait Trace<C> {
    /// Returns the trace without dividing by the local dimension.
    fn trace(&self) -> C;
}

/// Construction of a value representing the identity operator.
///
/// The value must act as identity and be a multiplicative unit. It need not be
/// part of a chosen basis; basis expansion may return a sum of basis primitives.
/// Both operator primitives and expression types can implement this contract.
pub trait Identity: Sized {
    /// Returns the identity without constructing numeric coefficients.
    fn identity() -> Self;
}

/// Dimension of the linear space spanned by all values of a primitive type.
///
/// For a complete state space this is `d`; for its full operator space it is
/// `d²`. It is not the number of primitive variants or the span of one value.
/// Dimensions must agree with the algebra, action, and complete basis contracts.
pub trait Dimension {
    /// Returns the positive dimension, or an error if it does not fit in usize.
    fn dim() -> Result<NonZeroUsize, DimensionOverflow>;
}
