use qmb_core::{BasisError, DimensionOverflow, ExpressionError, SpaceError};
use std::collections::TryReserveError;

/// Failure to construct or interpret a numerical export.
#[derive(Debug, thiserror::Error)]
pub enum ExportError {
    /// Invalid basis expansion or enumeration.
    #[error("{0}")]
    Basis(#[from] BasisError),
    /// Incompatible, unbound, or invalid expression space.
    #[error("{0}")]
    Expression(#[from] ExpressionError),
    /// Invalid tensor factor space.
    #[error("{0}")]
    Space(#[from] SpaceError),
    /// An axis dimension or flat coordinate overflows usize.
    #[error("{0}")]
    Dimension(#[from] DimensionOverflow),
    /// The output element count, byte size or compressed pointer length overflows.
    #[error("numerical export storage size overflows")]
    SizeOverflow,
    /// An allocation owned by this export could not be reserved.
    #[error("{0}")]
    Allocation(#[from] TryReserveError),
    /// The container shape or entry positions disagree with their metadata.
    #[error("numerical export shape disagrees with its data")]
    InvalidShape,
    /// An ndarray constructor rejected the shape.
    #[cfg(feature = "ndarray")]
    #[error("{0}")]
    NdarrayShape(#[from] ndarray::ShapeError),
    /// A sprs constructor rejected the compressed structure.
    #[cfg(feature = "sprs")]
    #[error("{0}")]
    SparseStructure(#[from] sprs::errors::StructureError),
}
