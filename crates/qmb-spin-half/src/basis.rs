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
/// Its spin elements are half the Pauli matrices, not normalized Pauli operators.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct CartesianBasis;

/// The operator basis in the order `I, SPlus, SMinus, Sz`.
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

impl<C: One + Neg<Output = C> + ImaginaryUnit> Basis<SpinHalfOperator, C> for CartesianBasis {
    fn expand(&self, primitive: &SpinHalfOperator) -> LinearCombination<SpinHalfOperator, C> {
        use SpinHalfOperator::*;
        match primitive {
            SPlus => {
                LinearCombination::primitive(Sx) + LinearCombination::term(imaginary(C::one()), Sy)
            }
            SMinus => {
                LinearCombination::primitive(Sx) + LinearCombination::term(imaginary(-C::one()), Sy)
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
            SPlus | SMinus => None,
        }
    }

    fn elements(&self) -> impl Iterator<Item = SpinHalfOperator> + '_ {
        use SpinHalfOperator::*;
        [I, Sx, Sy, Sz].into_iter()
    }
}

impl<C> Basis<SpinHalfOperator, C> for LadderBasis
where
    C: One + Add<Output = C> + Div<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    fn expand(&self, primitive: &SpinHalfOperator) -> LinearCombination<SpinHalfOperator, C> {
        use SpinHalfOperator::*;
        match primitive {
            Sx => LinearCombination::term(half(), SPlus) + LinearCombination::term(half(), SMinus),
            Sy => {
                LinearCombination::term(imaginary(-half::<C>()), SPlus)
                    + LinearCombination::term(imaginary(half()), SMinus)
            }
            _ => LinearCombination::primitive(*primitive),
        }
    }

    fn index(&self, primitive: &SpinHalfOperator) -> Option<usize> {
        use SpinHalfOperator::*;
        match primitive {
            I => Some(0),
            SPlus => Some(1),
            SMinus => Some(2),
            Sz => Some(3),
            Sx | Sy => None,
        }
    }

    fn elements(&self) -> impl Iterator<Item = SpinHalfOperator> + '_ {
        use SpinHalfOperator::*;
        [I, SPlus, SMinus, Sz].into_iter()
    }
}
