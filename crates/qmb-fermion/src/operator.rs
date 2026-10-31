//! Fermion and Majorana primitives and their local multiplication rules.

use crate::coefficient::imaginary;
use core::ops::{Add, Neg};
use num_traits::{One, Zero};
use qmb_core::{
    Adjoint, Dimension, DimensionOverflow, HomogeneousParity, Identity, ImaginaryUnit,
    LinearCombination, Parity, ParityOperator, Product, Trace,
};

/// Seven homogeneous primitives spanning the four-dimensional operator space.
///
/// Mx and My are standard Majorana operators, with `[Mx,My]=2iMz` and
/// `Mx²=My²=Mz²=I`. Mz is the even parity operator `I-2n`.
/// Structural equality distinguishes the linearly dependent primitive values.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FermionOperator {
    /// The identity operator.
    I,
    /// Annihilation `a = |0><1|`, with odd parity.
    Annihilate,
    /// Creation `a† = |1><0|`, with odd parity.
    Create,
    /// Occupation number `n = a†a`.
    Number,
    /// The first Majorana: `Mx = a+a†`, odd and Hermite.
    Mx,
    /// The second Majorana: `My = (a-a†)/i`, odd and Hermite.
    My,
    /// Parity `Mz = -i Mx My = I-2n`, a Hermitian even bilinear.
    Mz,
}

impl<C> Product<C> for FermionOperator
where
    C: One + Add<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    type Term = Self;
    fn product(&self, rhs: &Self) -> LinearCombination<Self, C> {
        use FermionOperator::{
            Annihilate as A, Create as D, I, Mx as X, My as Y, Mz as Z, Number as N,
        };
        use LinearCombination as Lc;
        // Direct rules: analytically zero terms are omitted; each result keeps basis order.
        match (*self, *rhs) {
            (I, op) | (op, I) => Lc::primitive(op),
            (A, A) | (D, D) | (D, N) | (N, A) => Lc::zero(),
            (A, D) | (A, X) | (X, D) => Lc::primitive(I) + Lc::term(-C::one(), N),
            (A, N) | (X, N) | (Z, A) => Lc::primitive(A),
            (D, A) | (N, N) | (D, X) | (X, A) => Lc::primitive(N),
            (N, D) | (D, Z) | (N, X) => Lc::primitive(D),
            (A, Y) => Lc::term(imaginary(C::one()), I) + Lc::term(imaginary(-C::one()), N),
            (A, Z) => Lc::term(-C::one(), A),
            (Z, D) => Lc::term(-C::one(), D),
            (D, Y) => Lc::term(imaginary(-C::one()), N),
            (N, Y) => Lc::term(imaginary(C::one()), D),
            (N, Z) | (Z, N) => Lc::term(-C::one(), N),
            (X, X) | (Y, Y) | (Z, Z) => Lc::primitive(I),
            (X, Y) => Lc::term(imaginary(C::one()), Z),
            (X, Z) => Lc::term(imaginary(-C::one()), Y),
            (Y, A) => Lc::term(imaginary(C::one()), N),
            (Y, D) => Lc::term(imaginary(-C::one()), I) + Lc::term(imaginary(C::one()), N),
            (Y, N) => Lc::term(imaginary(-C::one()), A),
            (Y, X) => Lc::term(imaginary(-C::one()), Z),
            (Y, Z) => Lc::term(imaginary(C::one()), X),
            (Z, X) => Lc::term(imaginary(C::one()), Y),
            (Z, Y) => Lc::term(imaginary(-C::one()), X),
        }
    }
}

impl<C: One> Adjoint<C> for FermionOperator {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, C> {
        LinearCombination::primitive(match self {
            Self::Annihilate => Self::Create,
            Self::Create => Self::Annihilate,
            _ => *self,
        })
    }
}

impl<C: Zero + One> Trace<C> for FermionOperator {
    fn trace(&self) -> C {
        match self {
            Self::I => C::one() + C::one(),
            Self::Number => C::one(),
            _ => C::zero(),
        }
    }
}

impl Identity for FermionOperator {
    fn identity() -> Self {
        Self::I
    }
}
impl Dimension for FermionOperator {
    fn dim() -> Result<core::num::NonZeroUsize, DimensionOverflow> {
        Ok(core::num::NonZeroUsize::new(4).unwrap())
    }
}
impl HomogeneousParity for FermionOperator {
    fn parity(&self) -> Parity {
        match self {
            Self::Annihilate | Self::Create | Self::Mx | Self::My => Parity::Odd,
            _ => Parity::Even,
        }
    }
}
impl ParityOperator for FermionOperator {
    fn parity_operator() -> Self {
        Self::Mz
    }
}
