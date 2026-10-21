//! Local basis expansion, indexing, and lazy enumeration.

use crate::{DimensionOverflow, LinearCombination, LocalStructure};

/// A complete, ordered basis represented using the target primitive type.
///
/// Implementations must preserve mathematical values on expansion. The basis
/// need not be orthonormal, and target primitives need not be independent.
/// Every valid target primitive must be expandable. Configurable basis types
/// must reject invalid configurations at construction.
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

/// A state basis for the indicated local structure.
///
/// Its number of elements must equal [`state_basis_len`].
pub trait StateBasis<L, T>: Basis<L::State, T>
where
    L: LocalStructure<T>,
{
}

/// A basis of the full operator space for the indicated local structure.
///
/// Its number of elements must equal [`operator_basis_len`].
pub trait OperatorBasis<L, T>: Basis<L::Operator, T>
where
    L: LocalStructure<T>,
{
}

/// Returns the number of elements in any complete state basis for `L`.
pub fn state_basis_len<L, T>() -> usize
where
    L: LocalStructure<T>,
{
    L::DIMENSION.get()
}

/// Returns the squared state dimension, the size of a complete operator basis.
///
/// # Errors
///
/// Returns [`DimensionOverflow`] if the square does not fit in `usize`.
pub fn operator_basis_len<L, T>() -> Result<usize, DimensionOverflow>
where
    L: LocalStructure<T>,
{
    let dimension = L::DIMENSION.get();
    dimension.checked_mul(dimension).ok_or(DimensionOverflow)
}
