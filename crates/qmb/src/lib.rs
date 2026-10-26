//! Quantum many-body algebras, spaces and site geometries.
//!
//! The default `normalize` feature enables norms and normalization extensions. Disable default
//! features to keep algebraic operations, including scalar division.
//! The legacy `libm` feature is an alias for `normalize`, which uses std.
//!
//! ```
//! use num_complex::Complex;
//! use qmb::prelude::*;
//! use qmb::{Divide, SpinHalfState, SpinHalfOperator};
//! let state = LocalState::<SpinHalfState, SpinHalfOperator, Complex<f64>>::primitive(SpinHalfState::Up);
//! let identity = LocalOperator::identity();
//! assert_eq!(state.norm_squared(), 1.0);
//! // This state has nonzero norm, so its scalar quotient is well-defined.
//! assert_eq!(state.expectation(&identity).divide(state.inner_product(&state)),
//!     Complex::new(1.0, 0.0));
//! # #[cfg(feature = "normalize")] {
//! assert_eq!(state.norm(), 1.0);
//! assert_eq!(state.normalized()?.norm(), 1.0);
//! assert_eq!(state.normalized_expectation(&identity)?, Complex::new(1.0, 0.0));
//! # }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// UUID-backed index generation.
pub use qmb_basic as basic;
/// Abstract algebra, expressions and Hilbert spaces.
pub use qmb_core as core;
/// Site assignments, chains and geometric relations.
pub use qmb_geometry as geometry;
/// Spin primitives, bases and model utilities.
pub use qmb_spin_half as spin_half;

pub use qmb_basic::Id128;
pub use qmb_core::{
    Basis, BasisError, Conjugate, Coordinates, Divide, ExpressionError, HilbertSpace, Identity,
    ImaginaryUnit, LocalOperator, LocalState, Operator, RealPart, SpanDimension, State,
};
pub use qmb_geometry::{Boundary, Chain, ChainError, NearestNeighborRelation};
pub use qmb_spin_half::{SpinHalfOperator, SpinHalfState};

pub mod space;
pub use space::ToSpace;

/// Norm, normalization and normalized expectation extensions (uses std).
#[cfg(feature = "normalize")]
pub use qmb_normalize as normalize;
#[cfg(feature = "normalize")]
pub use qmb_normalize::{
    Norm, NormalizationError, Normalize, NormalizedExpectation, NormalizedExpectationIn,
};

/// Common expression types and optional normalization extension traits.
pub mod prelude {
    pub use crate::{Basis, HilbertSpace, LocalOperator, LocalState, Operator, State, ToSpace};
    #[cfg(feature = "normalize")]
    pub use qmb_normalize::{Norm, Normalize, NormalizedExpectation, NormalizedExpectationIn};
}
