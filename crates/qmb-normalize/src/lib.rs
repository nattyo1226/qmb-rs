//! Norms, state normalization and normalized expectations.
//!
//! Import the corresponding extension traits to use these operations. Squared
//! norms, inner products, expectations and scalar division belong to `qmb-core`.
//! Only square-root-based operations need [`Float`](num_traits::Float) on the associated real value;
//! normalized expectations require neither square roots nor real-part extraction.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod error;
mod expectation;
mod norm;

pub use error::NormalizationError;
pub use expectation::{NormalizedExpectation, NormalizedExpectationIn};
pub use norm::{Norm, Normalize};
