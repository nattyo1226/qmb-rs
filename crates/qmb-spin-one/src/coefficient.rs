//! Small constant coefficients, constructed in the selected scalar type.

use core::ops::{Add, Div, Mul};
use num_traits::{Float, One};
use qmb_core::ImaginaryUnit;

pub(crate) fn imaginary<C: ImaginaryUnit + Mul<Output = C>>(value: C) -> C {
    C::imaginary_unit() * value
}
pub(crate) fn two<C: One + Add<Output = C>>() -> C {
    C::one() + C::one()
}
pub(crate) fn three<C: One + Add<Output = C>>() -> C {
    C::one() + C::one() + C::one()
}
pub(crate) fn four<C: One + Add<Output = C>>() -> C {
    two::<C>() + two::<C>()
}
pub(crate) fn half<C: One + Add<Output = C> + Div<Output = C>>() -> C {
    C::one() / two()
}
pub(crate) fn two_thirds<C: One + Add<Output = C> + Div<Output = C>>() -> C {
    two::<C>() / three::<C>()
}
pub(crate) fn four_ninths<C: One + Add<Output = C> + Div<Output = C>>() -> C {
    four::<C>() / (three::<C>() + three::<C>() + three::<C>())
}
pub(crate) fn sqrt_two<R: Float>() -> R {
    two::<R>().sqrt()
}
pub(crate) fn sqrt_three<R: Float>() -> R {
    three::<R>().sqrt()
}
pub(crate) fn spin_scale<R: Float>() -> R {
    (three::<R>() / two::<R>()).sqrt()
}
