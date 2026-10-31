//! Ordered spin-1 state and dipole/quadrupole operator bases.

use crate::{
    SpinOneOperator as O, SpinOneState as S,
    coefficient::{half, imaginary},
};
use core::ops::{Add, Div, Neg};
use num_traits::One;
use qmb_core::{Basis, ImaginaryUnit, LinearCombination};

/// The orthonormal state basis in the order `Plus, Zero, Minus`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct SzBasis;

/// The orthogonal Hermitian basis `I, Sx, Sy, Sz, Dxy, Q0, Qxy, Qyz, Qzx`.
///
/// Orthonormal under the dimension-normalized inner product `Tr(A† B)/3`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct CartesianBasis;

/// The orthogonal basis `I, Sp, Sm, Sz, Q0, Q1p, Q1m, Q2p, Q2m`.
///
/// Dimension-normalized squared norms are 1 for I, Sz and Q0, and 2
/// for the six ladder elements. No spherical-tensor phase factors are applied.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct LadderBasis;

impl<C: One> Basis<S, C> for SzBasis {
    fn expand(&self, primitive: &S) -> LinearCombination<S, C> {
        LinearCombination::primitive(*primitive)
    }
    fn index(&self, primitive: &S) -> Option<usize> {
        Some(match primitive {
            S::Plus => 0,
            S::Zero => 1,
            S::Minus => 2,
        })
    }
    fn elements(&self) -> impl Iterator<Item = S> + '_ {
        [S::Plus, S::Zero, S::Minus].into_iter()
    }
}

impl<C: One + Neg<Output = C> + ImaginaryUnit> Basis<O, C> for CartesianBasis {
    fn expand(&self, primitive: &O) -> LinearCombination<O, C> {
        use LinearCombination as Lc;
        use O::*;
        match primitive {
            Sp => Lc::primitive(Sx) + Lc::term(imaginary(C::one()), Sy),
            Sm => Lc::primitive(Sx) + Lc::term(imaginary(-C::one()), Sy),
            Q1p => Lc::primitive(Qzx) + Lc::term(imaginary(C::one()), Qyz),
            Q1m => Lc::primitive(Qzx) + Lc::term(imaginary(-C::one()), Qyz),
            Q2p => Lc::primitive(Dxy) + Lc::term(imaginary(C::one()), Qxy),
            Q2m => Lc::primitive(Dxy) + Lc::term(imaginary(-C::one()), Qxy),
            _ => Lc::primitive(*primitive),
        }
    }
    fn index(&self, primitive: &O) -> Option<usize> {
        use O::*;
        match primitive {
            I => Some(0),
            Sx => Some(1),
            Sy => Some(2),
            Sz => Some(3),
            Dxy => Some(4),
            Q0 => Some(5),
            Qxy => Some(6),
            Qyz => Some(7),
            Qzx => Some(8),
            Sp | Sm | Q1p | Q1m | Q2p | Q2m => None,
        }
    }
    fn elements(&self) -> impl Iterator<Item = O> + '_ {
        use O::*;
        [I, Sx, Sy, Sz, Dxy, Q0, Qxy, Qyz, Qzx].into_iter()
    }
}

impl<C> Basis<O, C> for LadderBasis
where
    C: One + Add<Output = C> + Div<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    fn expand(&self, primitive: &O) -> LinearCombination<O, C> {
        use LinearCombination as Lc;
        use O::*;
        match primitive {
            Sx => Lc::term(half(), Sp) + Lc::term(half(), Sm),
            Sy => Lc::term(imaginary(-half::<C>()), Sp) + Lc::term(imaginary(half()), Sm),
            Qzx => Lc::term(half(), Q1p) + Lc::term(half(), Q1m),
            Qyz => Lc::term(imaginary(-half::<C>()), Q1p) + Lc::term(imaginary(half()), Q1m),
            Dxy => Lc::term(half(), Q2p) + Lc::term(half(), Q2m),
            Qxy => Lc::term(imaginary(-half::<C>()), Q2p) + Lc::term(imaginary(half()), Q2m),
            _ => Lc::primitive(*primitive),
        }
    }
    fn index(&self, primitive: &O) -> Option<usize> {
        use O::*;
        match primitive {
            I => Some(0),
            Sp => Some(1),
            Sm => Some(2),
            Sz => Some(3),
            Q0 => Some(4),
            Q1p => Some(5),
            Q1m => Some(6),
            Q2p => Some(7),
            Q2m => Some(8),
            Sx | Sy | Dxy | Qxy | Qyz | Qzx => None,
        }
    }
    fn elements(&self) -> impl Iterator<Item = O> + '_ {
        use O::*;
        [I, Sp, Sm, Sz, Q0, Q1p, Q1m, Q2p, Q2m].into_iter()
    }
}
