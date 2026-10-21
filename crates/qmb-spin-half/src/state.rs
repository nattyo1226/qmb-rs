//! Spin-1/2 states, their orthonormal inner product, and operator actions.

use core::ops::{Div, Neg};

use qmb_core::{
    Act, Complex, InnerProduct, LinearCombination,
    num_traits::{One, Zero},
};

use crate::{
    SpinHalfOperator,
    coefficient::{half, imaginary, real},
};

/// The eigenstates of `Sz`, independent of coefficient precision.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpinHalfState {
    /// The state with `Sz` eigenvalue `+1/2`, first in the state basis.
    Up,
    /// The state with `Sz` eigenvalue `-1/2`, second in the state basis.
    Down,
}

impl<T: Zero + One> InnerProduct<T> for SpinHalfState {
    fn inner_product(&self, rhs: &Self) -> Complex<T> {
        real(if self == rhs { T::one() } else { T::zero() })
    }
}

impl<T> Act<SpinHalfOperator, T> for SpinHalfState
where
    T: Zero + One + Div<Output = T> + Neg<Output = T>,
{
    fn act(&self, operator: &SpinHalfOperator) -> LinearCombination<Self, T> {
        use LinearCombination as Lc;
        use SpinHalfOperator::*;
        use SpinHalfState::*;

        match (*operator, *self) {
            (I, state) => Lc::primitive(state),
            (Sx, Up) => Lc::term(real(half()), Down),
            (Sx, Down) => Lc::term(real(half()), Up),
            (Sy, Up) => Lc::term(imaginary(half()), Down),
            (Sy, Down) => Lc::term(imaginary(-half::<T>()), Up),
            (Sz, Up) => Lc::term(real(half()), Up),
            (Sz, Down) => Lc::term(real(-half::<T>()), Down),
            (SPlus, Up) | (SMinus, Down) => Lc::zero(),
            (SPlus, Down) => Lc::primitive(Up),
            (SMinus, Up) => Lc::primitive(Down),
        }
    }
}
