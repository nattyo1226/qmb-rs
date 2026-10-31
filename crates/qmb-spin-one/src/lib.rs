//! Spin-1 states, dipole and quadrupole operators, and ordered local bases.
//!
//! With `hbar = 1`, physical spin `Jz` has eigenvalues `1, 0, -1` in
//! state order `Plus, Zero, Minus`. Cartesian primitives are `S = sqrt(3/2) J`.
//! [`SpinOneState::Zero`] is a normalized state,
//! whereas a zero vector is an empty linear combination. The fifteen operator
//! values span a nine-dimensional space; their structural equality does not
//! identify different linear combinations representing the same operator.
//!
//! [`CartesianBasis`] and [`LadderBasis`] are complete, orthogonal operator
//! bases. Cartesian elements are orthonormal under `Tr(A†B)/3`; the ladder
//! basis has squared norms 1 for `I, Sz, Q0` and 2 for its six ladder elements.
//! Products of elements of the
//! same basis stay in that basis. Products of shared elements `I, Sz, Q0` stay
//! in their common span; mixed Cartesian/ladder products use Cartesian elements.
//! Identity multiplication preserves the other primitive.
//!
//! ```
//! use num_complex::Complex;
//! use qmb_core::Product;
//! use qmb_spin_one::{SpinOne, SpinOneOperator as Op, SpinOneState as Ket};
//! type Spin = SpinOne<Complex<f64>>;
//! let raised = Spin::state(Ket::Zero).act(&Spin::operator(Op::Sp));
//! let amplitude = Spin::state(Ket::Plus).inner_product(&raised);
//! assert!((amplitude.re - 3.0_f64.sqrt()).abs() < 1e-12);
//! let squared = <Op as Product<Complex<f64>>>::product(&Op::Sp, &Op::Sp);
//! assert_eq!(squared, qmb_core::LinearCombination::term(Complex::new((3.0_f64/2.0).sqrt(), 0.0), Op::Q2p));
//! ```
//!
//! Products and state actions use `Complex<R>` with `R: num_traits::Float`.
//! Basis conversions, adjoints and traces require only their individual capabilities: no `Float`,
//! coefficient `Clone`/`Copy`, casts or parsing. Scalars must represent the
//! required fractions and imaginary unit; integer division does not satisfy
//! that mathematical contract. Construction needs no numeric capabilities.
//! This crate permits `std`, forbids unsafe code, and uses core's shared
//! expressions, evaluation, tensor composition and basis coordinates.
//!
//! [`models::SpinOneModels::xyz`], [`models::SpinOneModels::heisenberg`] and [`models::SpinOneModels::blbq`] construct
//! chain Hamiltonians on caller-owned site spaces. Positive bilinear coupling
//! is ferromagnetic, following the spin-half model convention.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod basis;
mod coefficient;
pub mod models;
pub mod operator;
pub mod state;

pub use basis::{CartesianBasis, LadderBasis, SzBasis};
pub use models::UtilityError;
pub use operator::SpinOneOperator;
pub use state::SpinOneState;

/// Selects spin-1 primitives and a coefficient type without storing values.
///
/// Shared construction methods are inherent and require no trait import.
/// No precision is selected implicitly, and no numeric bound applies here.
///
/// ```
/// use num_complex::Complex;
/// use qmb_core::HilbertSpace;
/// use qmb_spin_one::{SpinOne, SpinOneOperator as Op, SpinOneState as Ket};
/// type Spin = SpinOne<Complex<f32>>;
/// let a = HilbertSpace::factor("a");
/// let b = HilbertSpace::factor("b");
/// let space = a.tensor(&b)?;
/// let psi = Spin::product_state(&space, [Ket::Plus, Ket::Minus])?;
/// let zz = Spin::operator_product_at([(&a, Op::Sz), (&b, Op::Sz)])?;
/// assert!((psi.expectation(&zz)?.re + 1.5).abs() < 1e-5);
/// # Ok::<(), qmb_core::ExpressionError>(())
/// ```
///
/// Different precisions cannot be mixed implicitly:
///
/// ```compile_fail
/// use num_complex::Complex;
/// use qmb_spin_one::{SpinOne, SpinOneOperator as Op, SpinOneState as Ket};
/// let psi = SpinOne::<Complex<f32>>::state(Ket::Zero);
/// let op = SpinOne::<Complex<f64>>::operator(Op::Sz);
/// let _ = psi.act(&op);
/// ```
pub type SpinOne<C> = qmb_core::ExpressionFactory<SpinOneState, SpinOneOperator, C>;
