//! Local basis expansion, indexing, and lazy enumeration.

use crate::{DimensionOverflow, LinearCombination, SpanDimension};

mod coordinates;
pub use coordinates::{BasisError, Coordinates};

/// Expansion, indexing and enumeration of a complete, ordered basis.
///
/// The basis need not be orthonormal, and target primitives need not be
/// independent. Configurable bases must reject invalid configurations at
/// construction. Implementations must expand every valid target primitive into
/// elements of this basis while preserving its mathematical value.
///
/// Numeric bounds belong to each concrete implementation, rather than this
/// trait. All three methods use that implementation's coefficient precision.
pub trait Basis<P, C> {
    /// Expands a primitive into a combination of elements of this basis.
    fn expand(&self, primitive: &P) -> LinearCombination<P, C>;

    /// Returns the zero-based index, or `None` for a non-basis primitive.
    fn index(&self, primitive: &P) -> Option<usize>;

    /// Yields each basis element exactly once in index order.
    ///
    /// Elements are owned; allocation and an exact-size iterator are not required.
    /// Use a qualified call when `C` cannot be inferred from context.
    fn elements(&self) -> impl Iterator<Item = P> + '_;
}

/// Returns the state basis size from the span of state primitive S.
pub fn state_basis_len<S: SpanDimension>() -> Result<usize, DimensionOverflow> {
    Ok(S::span_dimension()?.get())
}

/// Returns the operator basis size from the span of operator primitive O.
///
/// This is already the operator-space dimension, not a value to square again.
pub fn operator_basis_len<O: SpanDimension>() -> Result<usize, DimensionOverflow> {
    Ok(O::span_dimension()?.get())
}
