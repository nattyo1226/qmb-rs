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

//! Scaling accepts inputs convertible into the coefficient type selected at
//! construction. Numeric literals can be real while the stored C stays complex:
//!
//! ```
//! use num_complex::Complex;
//! use qmb_core::LocalOperator;
//! let op = LocalOperator::<u8, Complex<f32>>::primitive(0).scaled(0.5);
//! let same = op.scaled(Complex::new(1.0, 0.0));
//! ```
//!
//! Choosing C explicitly (or through [`QuantumSystem`]) replaces inference from
//! the scalar argument. Remove redundant `.into()` calls at the call site and
//! specify inputs such as `Complex::<f32>::default()` when needed.
//!
//! Real inputs also preserve precision; f64 is not narrowed to `Complex<f32>`:
//!
//! ```compile_fail
//! use num_complex::Complex;
//! use qmb_core::LocalState;
//! let state = LocalState::<u8, u8, Complex<f32>>::primitive(0);
//! let _ = state.scaled(0.5_f64);
//! ```
//!
//! ```compile_fail
//! use num_complex::Complex;
//! use qmb_core::{HilbertSpace, LocalOperator};
//! let factor = HilbertSpace::factor(0);
//! let op = LocalOperator::<u8, Complex<f32>>::primitive(0).at(&factor).unwrap();
//! let _ = op.scaled(Complex::<f64>::new(0.5, 0.0));
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
mod evaluate;

pub mod basis;
pub mod error;
pub mod expr;
pub mod ops;
pub mod space;
pub mod system;

pub use basis::{Basis, BasisError, Coordinates, operator_basis_len, state_basis_len};
pub use coefficient::{Conjugate, Divide, ImaginaryUnit, RealPart};
pub use error::DimensionOverflow;
pub use expr::parts::{IntoTerms, LinearCombination, Term};
pub use expr::{ExpressionError, LocalOperator, LocalState, Operator, State};
pub use ops::{Act, Adjoint, Dimension, Identity, InnerProduct, Product, Trace};
pub use space::{HilbertSpace, SpaceError};
pub use system::{LocalOperatorOf, LocalStateOf, OperatorOf, QuantumSystem, StateOf};
