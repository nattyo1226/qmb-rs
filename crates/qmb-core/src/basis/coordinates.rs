use alloc::{collections::BTreeMap, vec::Vec};

use crate::{DimensionOverflow, ExpressionError};

/// An invalid basis expansion, coordinate, or selected expression space.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum BasisError {
    /// A primitive returned by a basis rule has no index in that basis.
    #[error("basis expansion returned a non-basis primitive")]
    NonBasisElement,
    /// A basis index does not fit the primitive's span dimension.
    #[error("basis index {index} is outside the local dimension {dimension}")]
    IndexOutOfBounds {
        /// Invalid local basis index.
        index: usize,
        /// Local state-space or operator-space dimension.
        dimension: usize,
    },
    /// Enumeration and the basis index function disagree.
    #[error("basis enumeration disagrees with its indices or dimension")]
    InvalidEnumeration,
    /// A required dimension or coordinate cannot be represented as `usize`.
    #[error(transparent)]
    Dimension(#[from] DimensionOverflow),
    /// The expression has no ambient space or refers to an incompatible space.
    #[error(transparent)]
    Expression(#[from] ExpressionError),
}

/// Aggregated sparse coordinates in a borrowed, uniform local basis.
///
/// Each entry's index has one digit per factor, in tensor order. Digits are
/// local basis indices; lexicographic entry order is mixed-radix basis order.
/// A local expression has one digit; the unit space has the empty index.
/// Exact zero coefficients are omitted, without approximate pruning.
///
/// `H` is `()` for local expressions and [`crate::HilbertSpace<I>`] for placed
/// expressions. The latter shares the selected space without cloning indices.
/// Coordinates own their coefficients and borrow only the basis, not the source
/// expression. No dense vector, matrix, or collection of all basis elements is
/// constructed. The basis must retain the same interpretation while borrowed.
pub struct Coordinates<'b, B, C, H = ()> {
    pub(crate) basis: &'b B,
    pub(crate) space: H,
    pub(crate) factor_count: usize,
    pub(crate) dim_local: Option<usize>,
    pub(crate) entries: BTreeMap<Vec<usize>, C>,
}

impl<'b, B, C, H> Coordinates<'b, B, C, H> {
    /// Borrows the local basis applied to every factor.
    pub fn basis(&self) -> &'b B {
        self.basis
    }

    /// Borrows the selected space, or `()` for a local expression.
    pub fn space(&self) -> &H {
        &self.space
    }

    /// Returns the number of indexed factors (one for a local expression).
    pub fn factor_count(&self) -> usize {
        self.factor_count
    }

    /// Returns the local basis size, or `None` on the zero-factor unit space.
    pub fn dim_local(&self) -> Option<usize> {
        self.dim_local
    }

    /// Computes the full coordinate-space dimension, detecting overflow.
    ///
    /// This can fail while sparse entries and basis expansion remain usable.
    pub fn dim(&self) -> Result<usize, DimensionOverflow> {
        (0..self.factor_count).try_fold(1usize, |n, _| {
            n.checked_mul(self.dim_local.expect("nonempty factor space"))
                .ok_or(DimensionOverflow)
        })
    }

    /// Returns the number of stored nonzero coordinates.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns whether all coordinates are exactly zero.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Borrows unique aggregated entries in basis order, including their digits.
    pub fn entries(&self) -> impl ExactSizeIterator<Item = (&[usize], &C)> + '_ {
        self.entries.iter().map(|(index, c)| (index.as_slice(), c))
    }

    /// Converts each stored index to a flat index, detecting overflow per entry.
    ///
    /// This does not allocate a dense vector or require that the total dimension
    /// fit `usize`. The rightmost digit changes fastest; the empty index is 0.
    pub fn flat_entries(
        &self,
    ) -> impl ExactSizeIterator<Item = Result<(usize, &C), DimensionOverflow>> + '_ {
        self.entries().map(|(digits, c)| {
            let index = digits.iter().try_fold(0usize, |n, digit| {
                n.checked_mul(self.dim_local.expect("nonempty coordinate"))
                    .and_then(|n| n.checked_add(*digit))
                    .ok_or(DimensionOverflow)
            })?;
            Ok((index, c))
        })
    }

    /// Moves out the aggregated entries in basis order.
    pub fn into_entries(self) -> impl ExactSizeIterator<Item = (Vec<usize>, C)> {
        self.entries.into_iter()
    }
}
