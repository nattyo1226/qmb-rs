//! Basis-aware exports into owned vectors, ndarray arrays and sprs containers.
//!
//! Dense data always has a public `Vec` export path. Optional `ndarray` exports
//! reuse that path; optional `sprs` exports use sparse coordinates directly.
//! Results retain the basis and selected factor space needed to interpret them.
//!
//! ```
//! use num_complex::Complex;
//! use qmb_core::LocalState;
//! use qmb_export::{Export, Exported};
//! use qmb_spin_half::{SpinHalfState as S, SpinHalfOperator as O, SzBasis};
//! let psi = LocalState::<S,O,Complex<f64>>::primitive(S::Down);
//! let output: Exported<Vec<Complex<f64>>, _> = (&psi).export(&SzBasis)?;
//! assert_eq!(output.metadata().shape(), [2]);
//! assert_eq!(output.data().as_slice(), &[Complex::new(0.,0.),Complex::new(1.,0.)]);
//! # #[cfg(feature="ndarray")] {
//! let array: Exported<qmb_export::ndarray::Array1<_>, _> = output.export(())?;
//! assert_eq!(array.data().len(),2);
//! # }
//! # Ok::<(), qmb_export::ExportError>(())
//! ```
//!
//! Coefficient precision cannot change implicitly during export:
//!
//! ```compile_fail
//! use num_complex::Complex;
//! use qmb_core::LocalState;
//! use qmb_export::{Export, Exported};
//! use qmb_spin_half::{SpinHalfState as S, SpinHalfOperator as O, SzBasis};
//! let psi = LocalState::<S,O,Complex<f32>>::primitive(S::Up);
//! let output: Exported<Vec<Complex<f64>>, _> = (&psi).export(&SzBasis).unwrap();
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

#[cfg(feature = "ndarray")]
mod backend_ndarray;
#[cfg(feature = "sprs")]
mod backend_sprs;
mod data;
mod error;
mod matrix;
mod vector;

pub use data::{Exported, Layout, MatrixArgs, MatrixMetadata, VectorMetadata};
pub use error::ExportError;

#[cfg(feature = "ndarray")]
pub use ndarray;
#[cfg(feature = "sprs")]
pub use sprs;

/// Ability to export into `T` with the operation-specific arguments `Args`.
///
/// Expressions implement this for shared references. Owned coordinate and Vec
/// results can instead move their data into the target, without coefficient clones.
pub trait Export<T, Args = ()> {
    /// Failure of basis interpretation, space selection, sizing or construction.
    type Error;
    /// Converts this source without implicitly changing coefficient precision.
    fn export(self, args: Args) -> Result<T, Self::Error>;
}
