//! Spin primitives and direct product, adjoint, and trace rules.

use core::ops::{Add, Div, Neg};

use num_traits::{One, Zero};
use qmb_core::{
    Adjoint, Dimension, DimensionOverflow, Identity, ImaginaryUnit, LinearCombination, Product,
    Trace,
};

use crate::coefficient::{half, imaginary, quarter};

/// The identity, Cartesian spins, and ladder operators with `hbar = 1`.
///
/// These six primitives span the four-dimensional operator space and are
/// linearly dependent. Structural equality distinguishes their variants.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpinHalfOperator {
    /// The identity operator.
    I,
    /// Half the Pauli X matrix.
    Sx,
    /// Half the Pauli Y matrix.
    Sy,
    /// Half the Pauli Z matrix.
    Sz,
    /// The raising operator `Sx + i Sy`.
    SPlus,
    /// The lowering operator `Sx - i Sy`.
    SMinus,
}

impl<C> Product<C> for SpinHalfOperator
where
    C: One + Add<Output = C> + Div<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    type Term = Self;

    fn product(&self, rhs: &Self) -> LinearCombination<Self, C> {
        use LinearCombination as Lc;
        use SpinHalfOperator::*;

        match (*self, *rhs) {
            (I, op) | (op, I) => Lc::primitive(op),
            (Sx, Sx) | (Sy, Sy) | (Sz, Sz) => Lc::term(quarter(), I),
            (Sx, Sy) => Lc::term(imaginary(half()), Sz),
            (Sy, Sx) => Lc::term(imaginary(-half::<C>()), Sz),
            (Sy, Sz) => Lc::term(imaginary(half()), Sx),
            (Sz, Sy) => Lc::term(imaginary(-half::<C>()), Sx),
            (Sz, Sx) => Lc::term(imaginary(half()), Sy),
            (Sx, Sz) => Lc::term(imaginary(-half::<C>()), Sy),
            (Sx, SPlus) | (SMinus, Sx) => Lc::term(quarter(), I) + Lc::term(-half::<C>(), Sz),
            (Sx, SMinus) | (SPlus, Sx) => Lc::term(quarter(), I) + Lc::term(half(), Sz),
            (Sy, SPlus) => {
                Lc::term(imaginary(quarter()), I) + Lc::term(imaginary(-half::<C>()), Sz)
            }
            (Sy, SMinus) => {
                Lc::term(imaginary(-quarter::<C>()), I) + Lc::term(imaginary(-half::<C>()), Sz)
            }
            (SPlus, Sy) => Lc::term(imaginary(quarter()), I) + Lc::term(imaginary(half()), Sz),
            (SMinus, Sy) => {
                Lc::term(imaginary(-quarter::<C>()), I) + Lc::term(imaginary(half()), Sz)
            }
            (Sz, SPlus) => Lc::term(half(), SPlus),
            (Sz, SMinus) => Lc::term(-half::<C>(), SMinus),
            (SPlus, Sz) => Lc::term(-half::<C>(), SPlus),
            (SMinus, Sz) => Lc::term(half(), SMinus),
            (SPlus, SPlus) | (SMinus, SMinus) => Lc::zero(),
            (SPlus, SMinus) => Lc::term(half(), I) + Lc::primitive(Sz),
            (SMinus, SPlus) => Lc::term(half(), I) + Lc::term(-C::one(), Sz),
        }
    }
}

impl<C: One> Adjoint<C> for SpinHalfOperator {
    type Term = Self;

    fn adjoint(&self) -> LinearCombination<Self, C> {
        LinearCombination::primitive(match self {
            Self::SPlus => Self::SMinus,
            Self::SMinus => Self::SPlus,
            _ => *self,
        })
    }
}

impl<C: Zero + One> Trace<C> for SpinHalfOperator {
    fn trace(&self) -> C {
        if *self == Self::I {
            C::one() + C::one()
        } else {
            C::zero()
        }
    }
}

impl Identity for SpinHalfOperator {
    fn identity() -> Self {
        Self::I
    }
}
impl Dimension for SpinHalfOperator {
    fn dim() -> Result<core::num::NonZeroUsize, DimensionOverflow> {
        Ok(core::num::NonZeroUsize::new(4).unwrap())
    }
}
