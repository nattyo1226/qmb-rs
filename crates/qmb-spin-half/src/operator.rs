//! Spin primitives and direct product, adjoint, and trace rules.

use core::ops::{Add, Div, Neg};

use num_traits::{One, Zero};
use qmb_core::{
    Adjoint, Dimension, DimensionOverflow, Identity, ImaginaryUnit, LinearCombination, Product,
    Trace,
};

use crate::coefficient::{half, imaginary};

/// The identity, Pauli matrices, and spin-1/2 raising and lowering operators.
///
/// These six primitives span the four-dimensional operator space and are
/// linearly dependent. Structural equality distinguishes their variants.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpinHalfOperator {
    /// The identity operator.
    I,
    /// The Pauli X matrix.
    Sx,
    /// The Pauli Y matrix.
    Sy,
    /// The Pauli Z matrix.
    Sz,
    /// The raising operator `(Sx + i Sy)/2`.
    Sp,
    /// The lowering operator `(Sx - i Sy)/2`.
    Sm,
}

impl<C> Product<C> for SpinHalfOperator
where
    C: One + Add<Output = C> + Div<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    type Term = Self;
    fn product(&self, rhs: &Self) -> LinearCombination<Self, C> {
        use LinearCombination as Lc;
        use SpinHalfOperator::*;
        // Direct rules: analytically zero terms are omitted; each result keeps basis order.
        match (*self, *rhs) {
            (I, op) | (op, I) => Lc::primitive(op),
            (Sx, Sx) | (Sy, Sy) | (Sz, Sz) => Lc::primitive(I),
            (Sx, Sy) => Lc::term(imaginary(C::one()), Sz),
            (Sy, Sx) => Lc::term(imaginary(-C::one()), Sz),
            (Sy, Sz) => Lc::term(imaginary(C::one()), Sx),
            (Sz, Sy) => Lc::term(imaginary(-C::one()), Sx),
            (Sz, Sx) => Lc::term(imaginary(C::one()), Sy),
            (Sx, Sz) => Lc::term(imaginary(-C::one()), Sy),
            (Sx, Sp) | (Sm, Sx) | (Sm, Sp) => Lc::term(half::<C>(), I) + Lc::term(-half::<C>(), Sz),
            (Sx, Sm) | (Sp, Sx) | (Sp, Sm) => Lc::term(half::<C>(), I) + Lc::term(half::<C>(), Sz),
            (Sy, Sp) => Lc::term(imaginary(half::<C>()), I) + Lc::term(imaginary(-half::<C>()), Sz),
            (Sy, Sm) => {
                Lc::term(imaginary(-half::<C>()), I) + Lc::term(imaginary(-half::<C>()), Sz)
            }
            (Sp, Sy) => Lc::term(imaginary(half::<C>()), I) + Lc::term(imaginary(half::<C>()), Sz),
            (Sm, Sy) => Lc::term(imaginary(-half::<C>()), I) + Lc::term(imaginary(half::<C>()), Sz),
            (Sz, Sp) => Lc::primitive(Sp),
            (Sz, Sm) => Lc::term(-C::one(), Sm),
            (Sp, Sz) => Lc::term(-C::one(), Sp),
            (Sm, Sz) => Lc::primitive(Sm),
            (Sp, Sp) | (Sm, Sm) => Lc::zero(),
        }
    }
}

impl<C: One> Adjoint<C> for SpinHalfOperator {
    type Term = Self;

    fn adjoint(&self) -> LinearCombination<Self, C> {
        LinearCombination::primitive(match self {
            Self::Sp => Self::Sm,
            Self::Sm => Self::Sp,
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

impl qmb_core::HomogeneousParity for SpinHalfOperator {
    fn parity(&self) -> qmb_core::Parity {
        qmb_core::Parity::Even
    }
}
impl qmb_core::ParityOperator for SpinHalfOperator {
    fn parity_operator() -> Self {
        Self::I
    }
}
