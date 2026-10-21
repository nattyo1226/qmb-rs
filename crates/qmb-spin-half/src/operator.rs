//! Spin primitives and direct product, adjoint, and trace rules.

use core::ops::{Div, Neg};

use qmb_core::{
    Adjoint, Complex, LinearCombination, Product, Trace,
    num_traits::{One, Zero},
};

use crate::coefficient::{half, imaginary, quarter, real};

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

impl<T> Product<T> for SpinHalfOperator
where
    T: Zero + One + Div<Output = T> + Neg<Output = T>,
{
    fn product(&self, rhs: &Self) -> LinearCombination<Self, T> {
        use LinearCombination as Lc;
        use SpinHalfOperator::*;

        match (*self, *rhs) {
            (I, op) | (op, I) => Lc::primitive(op),
            (Sx, Sx) | (Sy, Sy) | (Sz, Sz) => Lc::term(real(quarter()), I),
            (Sx, Sy) => Lc::term(imaginary(half()), Sz),
            (Sy, Sx) => Lc::term(imaginary(-half::<T>()), Sz),
            (Sy, Sz) => Lc::term(imaginary(half()), Sx),
            (Sz, Sy) => Lc::term(imaginary(-half::<T>()), Sx),
            (Sz, Sx) => Lc::term(imaginary(half()), Sy),
            (Sx, Sz) => Lc::term(imaginary(-half::<T>()), Sy),
            (Sx, SPlus) | (SMinus, Sx) => {
                Lc::term(real(quarter()), I) + Lc::term(real(-half::<T>()), Sz)
            }
            (Sx, SMinus) | (SPlus, Sx) => Lc::term(real(quarter()), I) + Lc::term(real(half()), Sz),
            (Sy, SPlus) => {
                Lc::term(imaginary(quarter()), I) + Lc::term(imaginary(-half::<T>()), Sz)
            }
            (Sy, SMinus) => {
                Lc::term(imaginary(-quarter::<T>()), I) + Lc::term(imaginary(-half::<T>()), Sz)
            }
            (SPlus, Sy) => Lc::term(imaginary(quarter()), I) + Lc::term(imaginary(half()), Sz),
            (SMinus, Sy) => {
                Lc::term(imaginary(-quarter::<T>()), I) + Lc::term(imaginary(half()), Sz)
            }
            (Sz, SPlus) => Lc::term(real(half()), SPlus),
            (Sz, SMinus) => Lc::term(real(-half::<T>()), SMinus),
            (SPlus, Sz) => Lc::term(real(-half::<T>()), SPlus),
            (SMinus, Sz) => Lc::term(real(half()), SMinus),
            (SPlus, SPlus) | (SMinus, SMinus) => Lc::zero(),
            (SPlus, SMinus) => Lc::term(real(half()), I) + Lc::primitive(Sz),
            (SMinus, SPlus) => Lc::term(real(half()), I) + Lc::term(real(-T::one()), Sz),
        }
    }
}

impl<T: Zero + One> Adjoint<T> for SpinHalfOperator {
    fn adjoint(&self) -> LinearCombination<Self, T> {
        LinearCombination::primitive(match self {
            Self::SPlus => Self::SMinus,
            Self::SMinus => Self::SPlus,
            _ => *self,
        })
    }
}

impl<T: Zero + One> Trace<T> for SpinHalfOperator {
    fn trace(&self) -> Complex<T> {
        real(if *self == Self::I {
            T::one() + T::one()
        } else {
            T::zero()
        })
    }
}
