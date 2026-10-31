//! Complete ordered occupation, ladder, and Majorana bases.

use crate::{
    FermionOperator, FermionState,
    coefficient::{half, imaginary},
};
use core::ops::{Add, Div, Neg};
use num_traits::One;
use qmb_core::{Basis, ImaginaryUnit, LinearCombination};

/// The orthonormal state basis `[Empty, Occupied]`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct OccupationBasis;

/// The operator basis `[I, Annihilate, Create, Number]`.
///
/// Complete, but not orthogonal: `Tr(I† Number) = 1`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct LadderBasis;

/// The Hermite operator basis `[I, Mx, My, Mz]`.
///
/// Orthonormal under the dimension-normalized inner product `Tr(A† B)/2`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct MajoranaBasis;

impl<C: One> Basis<FermionState, C> for OccupationBasis {
    fn expand(&self, primitive: &FermionState) -> LinearCombination<FermionState, C> {
        LinearCombination::primitive(*primitive)
    }
    fn index(&self, primitive: &FermionState) -> Option<usize> {
        Some(match primitive {
            FermionState::Empty => 0,
            FermionState::Occupied => 1,
        })
    }
    fn elements(&self) -> impl Iterator<Item = FermionState> + '_ {
        [FermionState::Empty, FermionState::Occupied].into_iter()
    }
}

impl<C> Basis<FermionOperator, C> for LadderBasis
where
    C: One + Add<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    fn expand(&self, primitive: &FermionOperator) -> LinearCombination<FermionOperator, C> {
        use FermionOperator::*;
        use LinearCombination as Lc;
        match primitive {
            Mx => Lc::primitive(Annihilate) + Lc::primitive(Create),
            My => {
                Lc::term(imaginary(-C::one()), Annihilate) + Lc::term(imaginary(C::one()), Create)
            }
            Mz => Lc::primitive(I) + Lc::term(-(C::one() + C::one()), Number),
            _ => Lc::primitive(*primitive),
        }
    }
    fn index(&self, primitive: &FermionOperator) -> Option<usize> {
        use FermionOperator::*;
        match primitive {
            I => Some(0),
            Annihilate => Some(1),
            Create => Some(2),
            Number => Some(3),
            _ => None,
        }
    }
    fn elements(&self) -> impl Iterator<Item = FermionOperator> + '_ {
        use FermionOperator::*;
        [I, Annihilate, Create, Number].into_iter()
    }
}

impl<C> Basis<FermionOperator, C> for MajoranaBasis
where
    C: One + Add<Output = C> + Div<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    fn expand(&self, primitive: &FermionOperator) -> LinearCombination<FermionOperator, C> {
        use FermionOperator::*;
        use LinearCombination as Lc;
        match primitive {
            Annihilate => Lc::term(half(), Mx) + Lc::term(imaginary(half()), My),
            Create => Lc::term(half(), Mx) + Lc::term(imaginary(-half::<C>()), My),
            Number => Lc::term(half(), I) + Lc::term(-half::<C>(), Mz),
            _ => Lc::primitive(*primitive),
        }
    }
    fn index(&self, primitive: &FermionOperator) -> Option<usize> {
        use FermionOperator::*;
        match primitive {
            I => Some(0),
            Mx => Some(1),
            My => Some(2),
            Mz => Some(3),
            _ => None,
        }
    }
    fn elements(&self) -> impl Iterator<Item = FermionOperator> + '_ {
        use FermionOperator::*;
        [I, Mx, My, Mz].into_iter()
    }
}
