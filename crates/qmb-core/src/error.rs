//! Errors common to the abstract contracts.

/// The dimension of the full operator space does not fit in `usize`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("operator-space dimension overflows usize")]
pub struct DimensionOverflow;
