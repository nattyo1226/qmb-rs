//! Spin-1 states and direct operator action in the orthonormal Sz basis.

use crate::{
    SpinOneOperator,
    coefficient::{spin_scale, sqrt_three, sqrt_two, two},
};
use num_complex::Complex;
use num_traits::{Float, One, Zero};
use qmb_core::{Act, Dimension, DimensionOverflow, InnerProduct, LinearCombination};

/// The Sz eigenstates in the order `Plus, Zero, Minus`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpinOneState {
    /// Magnetic quantum number +1.
    Plus,
    /// Magnetic quantum number 0; a normalized state, not a zero vector.
    Zero,
    /// Magnetic quantum number -1.
    Minus,
}

impl<C: Zero + One> InnerProduct<C> for SpinOneState {
    fn inner_product(&self, rhs: &Self) -> C {
        if self == rhs { C::one() } else { C::zero() }
    }
}

impl<R: Float> Act<Complex<R>> for SpinOneState {
    type Operator = SpinOneOperator;
    type Term = Self;
    fn act(&self, operator: &SpinOneOperator) -> LinearCombination<Self, Complex<R>> {
        use LinearCombination as Lc;
        use SpinOneOperator::*;
        use SpinOneState::*;
        let real = |value| Complex::new(value, R::zero());
        let imaginary = |value| Complex::new(R::zero(), value);
        let half_root = || sqrt_three::<R>() / two::<R>();
        match (*operator, *self) {
            (I, state) => Lc::primitive(state),
            (Q0, state @ (Plus | Minus)) => Lc::term(real(sqrt_two::<R>() / two::<R>()), state),
            (Sz, Plus) | (Dxy, Minus) => Lc::term(real(spin_scale()), Plus),
            (Sz, Minus) => Lc::term(real(-spin_scale::<R>()), Minus),
            (Dxy, Plus) => Lc::term(real(spin_scale()), Minus),
            (Q0, Zero) => Lc::term(real(-sqrt_two::<R>()), Zero),
            (Sp, Zero) | (Q1p, Zero) => Lc::term(real(sqrt_three()), Plus),
            (Sp, Minus) | (Sm, Plus) | (Q1m, Plus) => Lc::term(real(sqrt_three()), Zero),
            (Sm, Zero) => Lc::term(real(sqrt_three()), Minus),
            (Q1p, Minus) => Lc::term(real(-sqrt_three::<R>()), Zero),
            (Q1m, Zero) => Lc::term(real(-sqrt_three::<R>()), Minus),
            (Sx, Plus | Minus) | (Qzx, Plus) => Lc::term(real(half_root()), Zero),
            (Qzx, Minus) => Lc::term(real(-half_root()), Zero),
            (Sx, Zero) => Lc::term(real(half_root()), Plus) + Lc::term(real(half_root()), Minus),
            (Qzx, Zero) => Lc::term(real(half_root()), Plus) + Lc::term(real(-half_root()), Minus),
            (Sy, Plus) | (Qyz, Plus | Minus) => Lc::term(imaginary(half_root()), Zero),
            (Sy, Minus) => Lc::term(imaginary(-half_root()), Zero),
            (Sy, Zero) => {
                Lc::term(imaginary(-half_root()), Plus) + Lc::term(imaginary(half_root()), Minus)
            }
            (Qyz, Zero) => {
                Lc::term(imaginary(-half_root()), Plus) + Lc::term(imaginary(-half_root()), Minus)
            }
            (Qxy, Plus) => Lc::term(imaginary(spin_scale()), Minus),
            (Qxy, Minus) => Lc::term(imaginary(-spin_scale::<R>()), Plus),
            (Q2p, Minus) => Lc::term(real(two::<R>() * spin_scale::<R>()), Plus),
            (Q2m, Plus) => Lc::term(real(two::<R>() * spin_scale::<R>()), Minus),
            (Sz | Dxy | Qxy, Zero)
            | (Sp | Q1p, Plus)
            | (Sm | Q1m, Minus)
            | (Q2p, Plus | Zero)
            | (Q2m, Zero | Minus) => Lc::zero(),
        }
    }
}

impl Dimension for SpinOneState {
    fn dim() -> Result<core::num::NonZeroUsize, DimensionOverflow> {
        Ok(core::num::NonZeroUsize::new(3).unwrap())
    }
}
impl qmb_core::HomogeneousParity for SpinOneState {
    fn parity(&self) -> qmb_core::Parity {
        qmb_core::Parity::Even
    }
}
