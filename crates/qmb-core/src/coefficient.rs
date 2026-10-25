//! Internal complex arithmetic with operation-specific scalar bounds.

use core::ops::{Add, Mul, Neg, Sub};
use num_traits::Zero;

use crate::Complex;

pub(crate) fn conjugate<T: Neg<Output = T>>(value: Complex<T>) -> Complex<T> {
    Complex::new(value.re, -value.im)
}

// Complex<T>'s multiplication implementation also requires Num, including
// division and remainder. These operations only need a ring's arithmetic.
pub(crate) fn multiply<T>(a: Complex<T>, b: Complex<T>) -> Complex<T>
where
    T: Clone + Add<Output = T> + Sub<Output = T> + Mul<Output = T>,
{
    Complex::new(
        a.re.clone() * b.re.clone() - a.im.clone() * b.im.clone(),
        a.re * b.im + a.im * b.re,
    )
}

pub(crate) fn add<T: Add<Output = T>>(a: Complex<T>, b: Complex<T>) -> Complex<T> {
    Complex::new(a.re + b.re, a.im + b.im)
}

pub(crate) fn zero<T: Zero>() -> Complex<T> {
    Complex::new(T::zero(), T::zero())
}
