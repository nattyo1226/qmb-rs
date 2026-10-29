//! Single-mode fermions, su(2)-normalized Majorana operators, and Hubbard chains.
//!
//! States are ordered `Empty, Occupied`. With `a = Annihilate`, the Hermite
//! operators are `X = (a + a†)/2`, `Y = (a - a†)/(2i)`, and
//! `Z = -2iXY = I/2 - n`. Thus `[X,Y] = iZ` and `X² = Y² = Z² = I/4`.
//! X and Y are half the usual Majorana operators; Z is an even bilinear,
//! not a third odd Majorana. The parity operator `P = I - 2n = 2Z` is retained
//! independently. Eight primitives span a four-dimensional operator space.
//!
//! Primitive types have no coefficient precision. [`Fermion<C>`] selects it
//! once, while [`qmb_core::ExpressionFactory`] constructs shared core expressions.
//! Numeric bounds belong to individual operations, not expression construction.
//! Complete operator bases are [`LadderBasis`] and [`MajoranaBasis`].
//!
//! ```
//! use num_complex::Complex;
//! use qmb_fermion::{Fermion, FermionOperator as O, FermionState as S};
//! type F = Fermion<Complex<f64>>;
//! let occupied = F::state(S::Occupied);
//! let a = F::operator(O::Annihilate);
//! assert_eq!(F::state(S::Empty).inner_product(&occupied.act(&a)), Complex::new(1., 0.));
//! assert_eq!(a.anticommutator(&a.adjoint()).trace(), Complex::new(2., 0.));
//! ```
//!
//! Multiple modes use **graded** spaces and the core's common evaluator. No
//! mode indices or Jordan–Wigner strings are stored in the local primitives.
//!
//! ```
//! use num_complex::Complex;
//! use qmb_core::HilbertSpace;
//! use qmb_fermion::{Fermion, FermionOperator as O, FermionState as S};
//! type F = Fermion<Complex<f32>>;
//! let up = HilbertSpace::factor("up");
//! let down = HilbertSpace::factor("down");
//! let space = HilbertSpace::tensor_all_g([up.clone(), down.clone()])?;
//! let double = F::product_state(&space, [S::Occupied, S::Occupied])?;
//! let up_only = F::product_state(&space, [S::Occupied, S::Empty])?;
//! let removed = double.act(&F::operator(O::Annihilate).at(&down)?)?;
//! assert_eq!(up_only.inner_product(&removed)?, Complex::new(-1., 0.));
//! # Ok::<(), qmb_core::ExpressionError>(())
//! ```
//!
//! Different precisions cannot mix implicitly:
//!
//! ```compile_fail
//! use num_complex::Complex;
//! use qmb_fermion::{Fermion, FermionOperator as O, FermionState as S};
//! let ket = Fermion::<Complex<f32>>::state(S::Empty);
//! let op = Fermion::<Complex<f64>>::operator(O::Create);
//! let _ = ket.act(&op);
//! ```

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

extern crate alloc;

pub mod basis;
mod coefficient;
pub mod models;
pub mod operator;
pub mod state;

pub use basis::{LadderBasis, MajoranaBasis, OccupationBasis};
pub use models::UtilityError;
pub use operator::FermionOperator;
pub use state::FermionState;

/// Selects fermion primitives and a coefficient type, without storing values.
///
/// `state` and `operator` require no trait import. No default
/// precision, numeric capability, or cloning capability is required here.
pub type Fermion<C> = qmb_core::ExpressionFactory<FermionState, FermionOperator, C>;
