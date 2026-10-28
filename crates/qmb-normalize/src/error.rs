use num_traits::Zero;
use qmb_core::ExpressionError;

/// Failure to normalize a state or an expectation value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum NormalizationError {
    /// The computed self-inner-product is exactly zero, including cancellation.
    #[error("state has zero norm")]
    ZeroNorm,
    /// Placement, binding or contraction failed.
    #[error(transparent)]
    Expression(#[from] ExpressionError),
}

pub(crate) fn nonzero<C: Zero>(value: &C) -> Result<(), NormalizationError> {
    if value.is_zero() {
        Err(NormalizationError::ZeroNorm)
    } else {
        Ok(())
    }
}
