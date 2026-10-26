//! Independent capabilities of local primitive types.
//!
//! A primitive's type and capability implementations define its interpretation.
//! Storage and lazy syntax require neither of these capabilities.

use crate::DimensionOverflow;
use core::num::NonZeroUsize;

/// Construction of the identity operator as a primitive value.
///
/// The value must act as identity and be a multiplicative unit. It need not be
/// part of a chosen basis; basis expansion may return a sum of basis primitives.
pub trait Identity: Sized {
    /// Returns the identity without constructing numeric coefficients.
    fn identity() -> Self;
}

/// Dimension of the linear space spanned by all values of a primitive type.
///
/// For a complete state space this is `d`; for its full operator space it is
/// `d²`. It is not the number of primitive variants or the span of one value.
/// Dimensions must agree with the algebra, action, and complete basis contracts.
pub trait SpanDimension {
    /// Returns the positive dimension, or an error if it does not fit in usize.
    fn span_dimension() -> Result<NonZeroUsize, DimensionOverflow>;
}
