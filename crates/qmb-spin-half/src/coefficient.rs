//! Constant coefficient construction without cloning, casts, or parsing.

use core::ops::{Add, Div};

use qmb_core::{
    Complex,
    num_traits::{One, Zero},
};

pub(crate) fn real<T: Zero>(value: T) -> Complex<T> {
    Complex::new(value, T::zero())
}

pub(crate) fn imaginary<T: Zero>(value: T) -> Complex<T> {
    Complex::new(T::zero(), value)
}

pub(crate) fn half<T: One + Add<Output = T> + Div<Output = T>>() -> T {
    T::one() / (T::one() + T::one())
}

pub(crate) fn quarter<T: One + Add<Output = T> + Div<Output = T>>() -> T {
    T::one() / (T::one() + T::one() + T::one() + T::one())
}
