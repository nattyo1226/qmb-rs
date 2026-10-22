//! Spin-1/2 primitives, local rules, and ordered bases.
//!
//! With `hbar = 1`, the Cartesian spins are half the Pauli matrices and
//! `SPlus = Sx + i Sy`, `SMinus = Sx - i Sy`. States are ordered `Up, Down`.
//! Primitive types carry no numeric precision; the coefficient type is inferred
//! from the result or specified through a qualified trait call.
//!
//! ```
//! use qmb_core::{Act, Complex, LinearCombination, Product};
//! use qmb_spin_half::{SpinHalfOperator as Op, SpinHalfState as State};
//!
//! let raised: LinearCombination<State, f64> = State::Down.act(&Op::SPlus);
//! assert_eq!(raised, LinearCombination::primitive(State::Up));
//!
//! let squared = <Op as Product<f32>>::product(&Op::Sx, &Op::Sx);
//! assert_eq!(squared, LinearCombination::term(Complex::new(0.25, 0.0), Op::I));
//! ```
//!
//! Local products retain natural primitives; basis expansion is a separate step:
//!
//! ```
//! use qmb_core::{Basis, Complex, LinearCombination, Product};
//! use qmb_spin_half::{CartesianBasis, SpinHalfOperator as Op};
//!
//! let product: LinearCombination<Op, f64> = Op::I.product(&Op::SPlus);
//! assert_eq!(product, LinearCombination::primitive(Op::SPlus));
//! let expanded = <CartesianBasis as Basis<Op, f64>>::expand(&CartesianBasis, &Op::SPlus);
//! assert_eq!(expanded, LinearCombination::primitive(Op::Sx)
//!     + LinearCombination::term(Complex::new(0.0, 1.0), Op::Sy));
//! let elements: Vec<_> = <CartesianBasis as Basis<Op, f64>>::elements(&CartesianBasis).collect();
//! assert_eq!(elements, [Op::I, Op::Sx, Op::Sy, Op::Sz]);
//! ```
//!
//! Numeric bounds are operation-specific. The numeric type must represent the
//! required real fractions; trait bounds alone cannot guarantee their laws.
//! No `Float`, `Copy`, `Clone`, `Num`, casts, or math backend is required.
//! This crate uses `qmb-core`'s finite linear combinations, and does not implement
//! expression evaluation or operator syntax.

//!
//! Shared local expressions construct actions lazily and evaluate without a basis:
//!
//! ```
//! use qmb_core::{Complex, LocalState, LocalOperator};
//! use qmb_spin_half::{SpinHalf, SpinHalfState as State, SpinHalfOperator as Op};
//! let up = LocalState::<SpinHalf, f64>::primitive(State::Up);
//! let sx = LocalOperator::primitive(Op::Sx);
//! let down = LocalState::primitive(State::Down);
//! let result = up.act(&sx).evaluate();
//! assert_eq!(down.inner_product(&result), Complex::new(0.5, 0.0));
//! assert_eq!(sx.product(&sx).trace(), Complex::new(0.5, 0.0));
//! ```

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod basis;
pub mod local_structure;
pub mod operator;
pub mod state;

mod coefficient;

pub use basis::{CartesianBasis, LadderBasis, SzBasis};
pub use local_structure::SpinHalf;
pub use operator::SpinHalfOperator;
pub use state::SpinHalfState;
