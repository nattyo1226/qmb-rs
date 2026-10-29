//! Fermion and Majorana primitives and their local multiplication rules.

use crate::coefficient::{half, imaginary, quarter};
use core::ops::{Add, Div, Neg};
use num_traits::{One, Zero};
use qmb_core::{
    Adjoint, Dimension, DimensionOverflow, HomogeneousParity, Identity, ImaginaryUnit,
    LinearCombination, Parity, ParityOperator, Product, Trace,
};

/// Eight homogeneous primitives spanning the four-dimensional operator space.
///
/// X and Y are half the standard Majorana operators, so `[X,Y]=iZ` and
/// `X²=Y²=Z²=I/4`, matching spin-half normalization. Z is an even bilinear.
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
    /// Half the first Majorana: `X = (a+a†)/2`, odd and Hermite.
    MajoranaX,
    /// Half the second Majorana: `Y = (a-a†)/(2i)`, odd and Hermite.
    MajoranaY,
    /// Hermite even bilinear `Z = -2iXY = I/2-n`, not an odd Majorana.
    MajoranaZ,
    /// The parity operator `P = I-2n = 2Z`.
    Parity,
}

impl<C> Product<C> for FermionOperator
where
    C: One + Add<Output = C> + Div<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    type Term = Self;
    fn product(&self, rhs: &Self) -> LinearCombination<Self, C> {
        use FermionOperator::{
            Annihilate as A, Create as D, I, MajoranaX as X, MajoranaY as Y, MajoranaZ as Z,
            Number as N, Parity as P,
        };
        use LinearCombination as Lc;
        match (*self, *rhs) {
            (I, op) | (op, I) => Lc::primitive(op),
            (A, A) | (D, D) | (D, N) | (N, A) => Lc::zero(),
            (A, D) => Lc::primitive(I) + Lc::term(-C::one(), N),
            (A, N) | (P, A) => Lc::primitive(A),
            (A, P) => Lc::term(-C::one(), A),
            (D, A) | (N, N) => Lc::primitive(N),
            (D, P) | (N, D) => Lc::primitive(D),
            (N, P) | (P, N) => Lc::term(-C::one(), N),
            (P, D) => Lc::term(-C::one(), D),
            (P, P) => Lc::primitive(I),
            (A, X) | (X, D) => Lc::term(half(), I) + Lc::term(-half::<C>(), N),
            (A, Y) => Lc::term(imaginary(half()), I) + Lc::term(imaginary(-half::<C>()), N),
            (A, Z) => Lc::term(-half::<C>(), A),
            (Z, D) => Lc::term(-half::<C>(), D),
            (D, X) | (X, A) => Lc::term(half(), N),
            (D, Y) => Lc::term(imaginary(-half::<C>()), N),
            (D, Z) | (N, X) => Lc::term(half(), D),
            (N, Y) => Lc::term(imaginary(half()), D),
            (N, Z) | (Z, N) => Lc::term(-half::<C>(), N),
            (P, X) => Lc::term(imaginary(C::one()), Y),
            (Y, P) => Lc::term(imaginary(C::one()), X),
            (P, Y) => Lc::term(imaginary(-C::one()), X),
            (X, P) => Lc::term(imaginary(-C::one()), Y),
            (P, Z) | (Z, P) => Lc::term(half(), I),
            (X, N) | (Z, A) => Lc::term(half(), A),
            (X, X) | (Y, Y) | (Z, Z) => Lc::term(quarter(), I),
            (X, Y) => Lc::term(imaginary(half()), Z),
            (X, Z) => Lc::term(imaginary(-half::<C>()), Y),
            (Y, A) => Lc::term(imaginary(half()), N),
            (Y, D) => Lc::term(imaginary(-half::<C>()), I) + Lc::term(imaginary(half()), N),
            (Y, N) => Lc::term(imaginary(-half::<C>()), A),
            (Y, X) => Lc::term(imaginary(-half::<C>()), Z),
            (Y, Z) => Lc::term(imaginary(half()), X),
            (Z, X) => Lc::term(imaginary(half()), Y),
            (Z, Y) => Lc::term(imaginary(-half::<C>()), X),
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
            Self::Annihilate | Self::Create | Self::MajoranaX | Self::MajoranaY => Parity::Odd,
            _ => Parity::Even,
        }
    }
}
impl ParityOperator for FermionOperator {
    fn parity_operator() -> Self {
        Self::Parity
    }
}
