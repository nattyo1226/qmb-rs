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
/// Squared Hilbert–Schmidt norms are 3 for I, 6 for Q0, and 2 otherwise.
/// These elements are not rescaled to unit norm.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct CartesianBasis;

/// The orthogonal basis `I, SPlus, SMinus, Sz, Q0, Q1Plus, Q1Minus, Q2Plus, Q2Minus`.
///
/// Squared Hilbert–Schmidt norms are 3 for I, 2 for Sz, 6 for Q0, and 4
/// otherwise. No spherical-tensor normalization or phase factors are applied.
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
            SPlus => Lc::primitive(Sx) + Lc::term(imaginary(C::one()), Sy),
            SMinus => Lc::primitive(Sx) + Lc::term(imaginary(-C::one()), Sy),
            Q1Plus => Lc::primitive(Qzx) + Lc::term(imaginary(C::one()), Qyz),
            Q1Minus => Lc::primitive(Qzx) + Lc::term(imaginary(-C::one()), Qyz),
            Q2Plus => Lc::primitive(Dxy) + Lc::term(imaginary(C::one()), Qxy),
            Q2Minus => Lc::primitive(Dxy) + Lc::term(imaginary(-C::one()), Qxy),
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
            SPlus | SMinus | Q1Plus | Q1Minus | Q2Plus | Q2Minus => None,
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
            Sx => Lc::term(half(), SPlus) + Lc::term(half(), SMinus),
            Sy => Lc::term(imaginary(-half::<C>()), SPlus) + Lc::term(imaginary(half()), SMinus),
            Qzx => Lc::term(half(), Q1Plus) + Lc::term(half(), Q1Minus),
            Qyz => Lc::term(imaginary(-half::<C>()), Q1Plus) + Lc::term(imaginary(half()), Q1Minus),
            Dxy => Lc::term(half(), Q2Plus) + Lc::term(half(), Q2Minus),
            Qxy => Lc::term(imaginary(-half::<C>()), Q2Plus) + Lc::term(imaginary(half()), Q2Minus),
            _ => Lc::primitive(*primitive),
        }
    }
    fn index(&self, primitive: &O) -> Option<usize> {
        use O::*;
        match primitive {
            I => Some(0),
            SPlus => Some(1),
            SMinus => Some(2),
            Sz => Some(3),
            Q0 => Some(4),
            Q1Plus => Some(5),
            Q1Minus => Some(6),
            Q2Plus => Some(7),
            Q2Minus => Some(8),
            Sx | Sy | Dxy | Qxy | Qyz | Qzx => None,
        }
    }
    fn elements(&self) -> impl Iterator<Item = O> + '_ {
        use O::*;
        [I, SPlus, SMinus, Sz, Q0, Q1Plus, Q1Minus, Q2Plus, Q2Minus].into_iter()
    }
}
