//! Ordered state, Cartesian-spin, and ladder-operator bases.

use core::ops::{Div, Neg};

use qmb_core::{
    Basis, LinearCombination,
    num_traits::{One, Zero},
};

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

impl<T: Zero + One> Basis<SpinHalfState, T> for SzBasis {
    fn expand(&self, primitive: &SpinHalfState) -> LinearCombination<SpinHalfState, T> {
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

impl<T: Zero + One + Neg<Output = T>> Basis<SpinHalfOperator, T> for CartesianBasis {
    fn expand(&self, primitive: &SpinHalfOperator) -> LinearCombination<SpinHalfOperator, T> {
        use SpinHalfOperator::*;
        match primitive {
            SPlus => {
                LinearCombination::primitive(Sx) + LinearCombination::term(imaginary(T::one()), Sy)
            }
            SMinus => {
                LinearCombination::primitive(Sx) + LinearCombination::term(imaginary(-T::one()), Sy)
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

impl<T> Basis<SpinHalfOperator, T> for LadderBasis
where
    T: Zero + One + Div<Output = T> + Neg<Output = T>,
{
    fn expand(&self, primitive: &SpinHalfOperator) -> LinearCombination<SpinHalfOperator, T> {
        use crate::coefficient::real;
        use SpinHalfOperator::*;
        match primitive {
            Sx => {
                LinearCombination::term(real(half()), SPlus)
                    + LinearCombination::term(real(half()), SMinus)
            }
            Sy => {
                LinearCombination::term(imaginary(-half::<T>()), SPlus)
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
