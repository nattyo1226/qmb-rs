//! Ordered state, Cartesian-spin, and ladder-operator bases.

use core::ops::{Add, Div, Neg};

use num_traits::One;
use qmb_core::{Basis, ImaginaryUnit, LinearCombination};

use crate::{
    SpinHalfOperator, SpinHalfState,
    coefficient::{half, imaginary},
};

/// The orthonormal `Sz` eigenbasis in the order `Up, Down`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct SzBasis;

/// The operator basis in the order `I, Sx, Sy, Sz`.
///
/// Its spin elements are Pauli matrices, orthonormal under `Tr(A†B)/2`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct CartesianBasis;

/// The operator basis in the order `I, Sp, Sm, Sz`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct LadderBasis;

impl<C: One> Basis<SpinHalfState, C> for SzBasis {
    fn expand(&self, primitive: &SpinHalfState) -> LinearCombination<SpinHalfState, C> {
        LinearCombination::primitive(*primitive)
    }

    fn index(&self, primitive: &SpinHalfState) -> Option<usize> {
        Some(match primitive {
            SpinHalfState::Up => 0,
            SpinHalfState::Down => 1,
        })
    }

    fn elements(&self) -> impl Iterator<Item = SpinHalfState> + '_ {
        [SpinHalfState::Up, SpinHalfState::Down].into_iter()
    }
}

impl<C: One + Add<Output = C> + Div<Output = C> + Neg<Output = C> + ImaginaryUnit>
    Basis<SpinHalfOperator, C> for CartesianBasis
{
    fn expand(&self, primitive: &SpinHalfOperator) -> LinearCombination<SpinHalfOperator, C> {
        use SpinHalfOperator::*;
        match primitive {
            Sp => {
                LinearCombination::term(half(), Sx) + LinearCombination::term(imaginary(half()), Sy)
            }
            Sm => {
                LinearCombination::term(half(), Sx)
                    + LinearCombination::term(imaginary(-half::<C>()), Sy)
            }
            _ => LinearCombination::primitive(*primitive),
        }
    }

    fn index(&self, primitive: &SpinHalfOperator) -> Option<usize> {
        use SpinHalfOperator::*;
        match primitive {
            I => Some(0),
            Sx => Some(1),
            Sy => Some(2),
            Sz => Some(3),
            Sp | Sm => None,
        }
    }

    fn elements(&self) -> impl Iterator<Item = SpinHalfOperator> + '_ {
        use SpinHalfOperator::*;
        [I, Sx, Sy, Sz].into_iter()
    }
}

impl<C> Basis<SpinHalfOperator, C> for LadderBasis
where
    C: One + Neg<Output = C> + ImaginaryUnit,
{
    fn expand(&self, primitive: &SpinHalfOperator) -> LinearCombination<SpinHalfOperator, C> {
        use SpinHalfOperator::*;
        match primitive {
            Sx => LinearCombination::primitive(Sp) + LinearCombination::primitive(Sm),
            Sy => {
                LinearCombination::term(imaginary(-C::one()), Sp)
                    + LinearCombination::term(imaginary(C::one()), Sm)
            }
            _ => LinearCombination::primitive(*primitive),
        }
    }

    fn index(&self, primitive: &SpinHalfOperator) -> Option<usize> {
        use SpinHalfOperator::*;
        match primitive {
            I => Some(0),
            Sp => Some(1),
            Sm => Some(2),
            Sz => Some(3),
            Sx | Sy => None,
        }
    }

    fn elements(&self) -> impl Iterator<Item = SpinHalfOperator> + '_ {
        use SpinHalfOperator::*;
        [I, Sp, Sm, Sz].into_iter()
    }
}
