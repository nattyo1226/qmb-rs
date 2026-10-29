//! Spin-1/2 states, their orthonormal inner product, and operator actions.

use core::ops::{Add, Div, Neg};

use num_traits::{One, Zero};
use qmb_core::{Act, Dimension, DimensionOverflow, ImaginaryUnit, InnerProduct, LinearCombination};

use crate::{
    SpinHalfOperator,
    coefficient::{half, imaginary},
};

/// The eigenstates of `Sz`, independent of coefficient precision.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpinHalfState {
    /// The state with `Sz` eigenvalue `+1/2`, first in the state basis.
    Up,
    /// The state with `Sz` eigenvalue `-1/2`, second in the state basis.
    Down,
}

impl<C: Zero + One> InnerProduct<C> for SpinHalfState {
    fn inner_product(&self, rhs: &Self) -> C {
        if self == rhs { C::one() } else { C::zero() }
    }
}

impl<C> Act<C> for SpinHalfState
where
    C: One + Add<Output = C> + Div<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    type Operator = SpinHalfOperator;
    type Term = Self;

    fn act(&self, operator: &SpinHalfOperator) -> LinearCombination<Self, C> {
        use LinearCombination as Lc;
        use SpinHalfOperator::*;
        use SpinHalfState::*;

        match (*operator, *self) {
            (I, state) => Lc::primitive(state),
            (Sx, Up) => Lc::term(half(), Down),
            (Sx, Down) => Lc::term(half(), Up),
            (Sy, Up) => Lc::term(imaginary(half()), Down),
            (Sy, Down) => Lc::term(imaginary(-half::<C>()), Up),
            (Sz, Up) => Lc::term(half(), Up),
            (Sz, Down) => Lc::term(-half::<C>(), Down),
            (SPlus, Up) | (SMinus, Down) => Lc::zero(),
            (SPlus, Down) => Lc::primitive(Up),
            (SMinus, Up) => Lc::primitive(Down),
        }
    }
}

impl Dimension for SpinHalfState {
    fn dim() -> Result<core::num::NonZeroUsize, DimensionOverflow> {
        Ok(core::num::NonZeroUsize::new(2).unwrap())
    }
}

impl qmb_core::HomogeneousParity for SpinHalfState {
    fn parity(&self) -> qmb_core::Parity {
        qmb_core::Parity::Even
    }
}
