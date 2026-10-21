//! Spin-1/2 state primitives and their orthonormal inner product.

use qmb_core::{
    Complex, InnerProduct,
    num_traits::{One, Zero},
};

use crate::coefficient::real;

/// The eigenstates of `Sz`, independent of coefficient precision.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpinHalfState {
    /// The state with `Sz` eigenvalue `+1/2`, first in the state basis.
    Up,
    /// The state with `Sz` eigenvalue `-1/2`, second in the state basis.
    Down,
}

impl<T: Zero + One> InnerProduct<T> for SpinHalfState {
    fn inner_product(&self, rhs: &Self) -> Complex<T> {
        real(if self == rhs { T::one() } else { T::zero() })
    }
}
