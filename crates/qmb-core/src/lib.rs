//! Abstract algebra contracts, shared expressions and ordinary tensor spaces.
//!
//! Expressions hold the whole coefficient `C`, including [`Complex<T>`](num_complex::Complex).
//! Numeric bounds belong to individual operations; primitive types need not be
//! tied to a numeric precision.
//! This crate uses `alloc` but does not require `std` or a math-function backend.
//! Local and composite expressions share immutable nodes and evaluate using
//! primitive rules. [`HilbertSpace`] holds ordered factors with generic indices;
//! index generation and concrete algebras belong to separate crates.
//!
//! ```
//! use num_complex::Complex;
//! use qmb_core::LinearCombination;
//!
//! let sum = LinearCombination::<_, Complex<f64>>::primitive("p")
//!     + LinearCombination::term(Complex::new(0.0, 0.0), "p");
//! assert_eq!(sum.len(), 2); // Duplicate primitives and zero coefficients remain.
//! let scaled = sum.scaled(Complex::new(0.0, 1.0));
//! assert_eq!(scaled.terms().next().unwrap().coefficient, Complex::new(0.0, 1.0));
//! ```
//!
//! Different coefficient precisions cannot be mixed implicitly:
//!
//! ```compile_fail
//! use num_complex::Complex;
//! use qmb_core::LinearCombination;
//! let single = LinearCombination::<u8, Complex<f32>>::primitive(0);
//! let double = LinearCombination::<u8, Complex<f64>>::primitive(0);
//! let _ = single + double;
//! ```
//!
//! Scaling also preserves the coefficient precision:
//!
//! ```compile_fail
//! use num_complex::Complex;
//! use qmb_core::LinearCombination;
//! let single = LinearCombination::<u8, Complex<f32>>::primitive(0);
//! let _ = single.scaled(Complex::<f64>::new(2.0, 0.0));
//! ```

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

extern crate alloc;

#[cfg(test)]
extern crate self as qmb_core;
#[cfg(test)]
extern crate std;
#[cfg(test)]
#[path = "../tests/support/coefficient.rs"]
mod test_coefficient;

pub mod coefficient;

pub mod basis;
pub mod error;
pub mod expr;
pub mod linear_combination;
pub mod ops;
pub mod primitive;
pub mod space;
pub mod system;

pub use basis::{Basis, BasisError, Coordinates, operator_basis_len, state_basis_len};
pub use coefficient::{Conjugate, Divide, ImaginaryUnit, RealPart};
pub use error::DimensionOverflow;
pub use expr::{ExpressionError, LocalOperator, LocalState, Operator, State};
pub use linear_combination::{IntoTerms, LinearCombination, Term};
pub use ops::{Act, Adjoint, InnerProduct, Product, Trace};
pub use primitive::{Identity, SpanDimension};
pub use space::{HilbertSpace, SpaceError};
pub use system::{LocalOperatorOf, LocalStateOf, OperatorOf, QuantumSystem, StateOf};
