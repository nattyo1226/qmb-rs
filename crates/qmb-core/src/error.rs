//! Errors common to the abstract contracts.

/// A state-space or operator-space dimension does not fit in `usize`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("space dimension overflows usize")]
pub struct DimensionOverflow;
