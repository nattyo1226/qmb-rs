//! Local basis expansion, indexing, and lazy enumeration.

use crate::{DimensionOverflow, LinearCombination, LocalStructure};

/// Expansion, indexing and enumeration of a complete, ordered basis.
///
/// The basis need not be orthonormal, and target primitives need not be
/// independent. Configurable bases must reject invalid configurations at
/// construction. Implementations must expand every valid target primitive into
/// elements of this basis while preserving its mathematical value.
///
/// Numeric bounds belong to each concrete implementation, rather than this
/// trait. All three methods use that implementation's coefficient precision.
pub trait Basis<P, T> {
    /// Expands a primitive into a combination of elements of this basis.
    fn expand(&self, primitive: &P) -> LinearCombination<P, T>;

    /// Returns the zero-based index, or `None` for a non-basis primitive.
    fn index(&self, primitive: &P) -> Option<usize>;

    /// Yields each basis element exactly once in index order.
    ///
    /// Elements are owned; allocation and an exact-size iterator are not required.
    /// Use a qualified call when `T` cannot be inferred from context.
    fn elements(&self) -> impl Iterator<Item = P> + '_;
}

/// Returns the number of elements in any complete state basis for `L`.
pub fn state_basis_len<L>() -> usize
where
    L: LocalStructure,
{
    L::DIMENSION.get()
}

/// Returns the squared state dimension, the size of a complete operator basis.
///
/// # Errors
///
/// Returns [`DimensionOverflow`] if the square does not fit in `usize`.
pub fn operator_basis_len<L>() -> Result<usize, DimensionOverflow>
where
    L: LocalStructure,
{
    let dimension = L::DIMENSION.get();
    dimension.checked_mul(dimension).ok_or(DimensionOverflow)
}
