//! The fixed two-dimensional spin-1/2 local structure.

use core::num::NonZeroUsize;

use qmb_core::{
    Identity, LinearCombination, LocalStructure,
    num_traits::{One, Zero},
};

use crate::{SpinHalfOperator, SpinHalfState};

/// Spin-1/2 with `hbar = 1`, with no runtime definition parameters.
///
/// All supported precisions use the same primitives and dimension. The numeric
/// type must represent the real constants `1/2` and `1/4` appropriately; integer
/// division does not satisfy this contract. Trait bounds specify operations,
/// not a proof of the required numeric laws.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct SpinHalf;

impl LocalStructure for SpinHalf {
    type State = SpinHalfState;
    type Operator = SpinHalfOperator;
    const DIMENSION: NonZeroUsize = NonZeroUsize::new(2).unwrap();
}

impl<T: Zero + One> Identity<T> for SpinHalf {
    fn identity() -> LinearCombination<Self::Operator, T> {
        LinearCombination::primitive(SpinHalfOperator::I)
    }
}
