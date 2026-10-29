//! Constant coefficients constructed without casts, parsing, or cloning.

use core::ops::{Add, Div, Mul};
use num_traits::One;
use qmb_core::ImaginaryUnit;

pub(crate) fn half<C: One + Add<Output = C> + Div<Output = C>>() -> C {
    C::one() / (C::one() + C::one())
}

pub(crate) fn quarter<C: One + Add<Output = C> + Div<Output = C>>() -> C {
    C::one() / (C::one() + C::one() + C::one() + C::one())
}

pub(crate) fn imaginary<C: ImaginaryUnit + Mul<Output = C>>(value: C) -> C {
    C::imaginary_unit() * value
}
