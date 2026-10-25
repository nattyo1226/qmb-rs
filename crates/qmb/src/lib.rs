//! Quantum many-body algebras, spaces and site geometries.
//!
//! The default `libm` feature enables square-root-based norms. Disable default
//! features to keep only algebraic operations and normalized expectations.

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
    Complex, ExpressionError, HilbertSpace, LocalOperator, LocalState, LocalStructure,
    NormalizationError, Operator, State,
};
pub use qmb_geometry::{Boundary, Chain, Geometry, GeometryError, NearestNeighborRelation};
pub use qmb_spin_half::{SpinHalf, SpinHalfOperator, SpinHalfState};
