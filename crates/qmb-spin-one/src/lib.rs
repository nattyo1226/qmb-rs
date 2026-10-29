//! Spin-1 states, dipole and quadrupole operators, and ordered local bases.
//!
//! With `hbar = 1`, states are ordered `Plus, Zero, Minus`, with `Sz`
//! eigenvalues `1, 0, -1`. [`SpinOneState::Zero`] is a normalized state,
//! whereas a zero vector is an empty linear combination. The fifteen operator
//! values span a nine-dimensional space; their structural equality does not
//! identify different linear combinations representing the same operator.
//!
//! [`CartesianBasis`] and [`LadderBasis`] are complete, orthogonal operator
//! bases, without Hilbert–Schmidt normalization. Products of elements of the
//! same basis stay in that basis. Products of shared elements `I, Sz, Q0` stay
//! in their common span; mixed Cartesian/ladder products use Cartesian elements.
//! Identity multiplication preserves the other primitive.
//!
//! ```
//! use num_complex::Complex;
//! use qmb_core::{Product, QuantumSystem};
//! use qmb_spin_one::{SpinOne, SpinOneOperator as Op, SpinOneState as Ket};
//! type Spin = SpinOne<Complex<f64>>;
//! let raised = Spin::state(Ket::Zero).act(&Spin::operator(Op::SPlus));
//! let amplitude = Spin::state(Ket::Plus).inner_product(&raised);
//! assert!((amplitude.re - 2.0_f64.sqrt()).abs() < 1e-12);
//! let squared = <Op as Product<Complex<f64>>>::product(&Op::SPlus, &Op::SPlus);
//! assert_eq!(squared, qmb_core::LinearCombination::primitive(Op::Q2Plus));
//! ```
//!
//! State action uses `Complex<R>` with `R: num_traits::Float`. Other local
//! rules require only their individual scalar capabilities: no `Float`,
//! coefficient `Clone`/`Copy`, casts or parsing. Scalars must represent the
//! required fractions and imaginary unit; integer division does not satisfy
//! that mathematical contract. Construction needs no numeric capabilities.
//! This crate permits `std`, forbids unsafe code, and uses core's shared
//! expressions, evaluation, tensor composition and basis coordinates.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod basis;
mod coefficient;
pub mod operator;
pub mod state;

pub use basis::{CartesianBasis, LadderBasis, SzBasis};
pub use operator::SpinOneOperator;
pub use state::SpinOneState;

/// Selects spin-1 primitives and a coefficient type without storing values.
///
/// Import [`qmb_core::QuantumSystem`] for the shared construction methods.
/// No precision is selected implicitly, and no numeric bound applies here.
///
/// ```
/// use num_complex::Complex;
/// use qmb_core::{HilbertSpace, QuantumSystem};
/// use qmb_spin_one::{SpinOne, SpinOneOperator as Op, SpinOneState as Ket};
/// type Spin = SpinOne<Complex<f32>>;
/// let a = HilbertSpace::factor("a");
/// let b = HilbertSpace::factor("b");
/// let psi = Spin::state_tensor_at([(&a, Ket::Plus), (&b, Ket::Minus)])?;
/// let zz = Spin::operator_product_at([(&a, Op::Sz), (&b, Op::Sz)])?;
/// assert_eq!(psi.expectation(&zz)?, Complex::new(-1.0, 0.0));
/// # Ok::<(), qmb_core::ExpressionError>(())
/// ```
///
/// Different precisions cannot be mixed implicitly:
///
/// ```compile_fail
/// use num_complex::Complex;
/// use qmb_core::QuantumSystem;
/// use qmb_spin_one::{SpinOne, SpinOneOperator as Op, SpinOneState as Ket};
/// let psi = SpinOne::<Complex<f32>>::state(Ket::Zero);
/// let op = SpinOne::<Complex<f64>>::operator(Op::Sz);
/// let _ = psi.act(&op);
/// ```
pub struct SpinOne<C>(core::marker::PhantomData<fn() -> C>);

impl<C> qmb_core::QuantumSystem for SpinOne<C> {
    type State = SpinOneState;
    type Operator = SpinOneOperator;
    type Coeff = C;
}
