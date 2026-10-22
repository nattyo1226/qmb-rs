//! Abstract local-algebra contracts and finite linear combinations.
//!
//! Coefficients use [`Complex<T>`](Complex). Numeric bounds belong to individual
//! operations; primitive types need not be tied to a numeric precision.
//! This crate uses `alloc` but does not require `std` or a math-function backend.
//! Local expressions share immutable nodes and evaluate using primitive rules.
//!
//! ```
//! use qmb_core::{Complex, LinearCombination};
//!
//! let sum = LinearCombination::<_, f64>::primitive("p")
//!     + LinearCombination::term(Complex::new(0.0, 0.0), "p");
//! assert_eq!(sum.len(), 2); // Duplicate primitives and zero coefficients remain.
//! let scaled = sum.scaled(Complex::new(0.0, 1.0));
//! assert_eq!(scaled.terms().next().unwrap().coefficient, Complex::new(0.0, 1.0));
//! ```
//!
//! Different coefficient precisions cannot be mixed implicitly:
//!
//! ```compile_fail
//! use qmb_core::LinearCombination;
//! let single = LinearCombination::<u8, f32>::primitive(0);
//! let double = LinearCombination::<u8, f64>::primitive(0);
//! let _ = single + double;
//! ```
//!
//! Scaling also preserves the coefficient precision:
//!
//! ```compile_fail
//! use qmb_core::{Complex, LinearCombination};
//! let single = LinearCombination::<u8, f32>::primitive(0);
//! let _ = single.scaled(Complex::<f64>::new(2.0, 0.0));
//! ```

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

extern crate alloc;

#[cfg(test)]
extern crate std;

pub mod basis;
pub mod error;
pub mod expr;
pub mod linear_combination;
pub mod local_structure;
pub mod ops;

pub use basis::{Basis, operator_basis_len, state_basis_len};
pub use error::DimensionOverflow;
pub use expr::{LocalOperator, LocalState};
pub use linear_combination::{IntoTerms, LinearCombination, Term};
pub use local_structure::{Identity, LocalStructure};
pub use num_complex::Complex;
pub use num_traits;
pub use ops::{Act, Adjoint, InnerProduct, Product, Trace};
