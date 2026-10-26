//! Independent coefficient capabilities, requested only by operations using them.
//!
//! Coefficients are commutative scalars. Implementations must preserve the unit,
//! distributive, conjugation and real-embedding laws required by their algebra,
//! up to the precision of the numeric representation. Storage has no such bound.

use core::ops::{Add, Div, Mul, Neg, Sub};
use num_complex::Complex;
use num_traits::{One, Zero};

/// Conjugation from a shared reference, independent of other numeric abilities.
///
/// Conjugation must be involutive and preserve addition and multiplication.
pub trait Conjugate {
    /// Result of conjugation; evaluators require this to be the coefficient type.
    type Output;
    /// Returns the complex conjugate without modifying the input.
    fn conjugate(&self) -> Self::Output;
}

/// Extraction of the real component, used by squared norms.
pub trait RealPart {
    /// Real value associated with this coefficient representation.
    type Real;
    /// Moves out the real component, without clamping or approximation.
    fn real_part(self) -> Self::Real;
}

/// Construction of the imaginary unit, required by particular local algebras.
///
/// The returned value must square to minus one. It is not required for storing
/// expressions or evaluating an algebra whose rules use only real coefficients.
pub trait ImaginaryUnit {
    /// Returns the scalar `i`.
    fn imaginary_unit() -> Self;
}

/// Division of scalar coefficients.
///
/// This separate capability lets a representation preserve its numerical division
/// algorithm. In particular, complex division by a real self-inner-product must
/// not unnecessarily square that denominator. Implementations for ordinary real
/// scalars can forward to [`Div`]. The caller checks for an exactly zero divisor.
pub trait Divide<Rhs = Self> {
    /// Result of division.
    type Output;
    /// Divides by a nonzero coefficient.
    fn divide(self, rhs: Rhs) -> Self::Output;
}

impl<T: Clone + Neg<Output = T>> Conjugate for Complex<T> {
    type Output = Self;
    fn conjugate(&self) -> Self {
        Self::new(self.re.clone(), -self.im.clone())
    }
}
impl<T> RealPart for Complex<T> {
    type Real = T;
    fn real_part(self) -> T {
        self.re
    }
}
impl<T: Zero + One> ImaginaryUnit for Complex<T> {
    fn imaginary_unit() -> Self {
        Self::new(T::zero(), T::one())
    }
}
impl<T> Divide for Complex<T>
where
    T: Clone + Zero + Sub<Output = T> + Mul<Output = T> + Div<Output = T>,
{
    type Output = Self;
    fn divide(self, rhs: Self) -> Self {
        if rhs.im.is_zero() {
            return Self::new(self.re / rhs.re.clone(), self.im / rhs.re);
        }
        let denominator = rhs.re.clone() * rhs.re.clone() + rhs.im.clone() * rhs.im.clone();
        Self::new(
            (self.re.clone() * rhs.re.clone() + self.im.clone() * rhs.im.clone())
                / denominator.clone(),
            (self.im * rhs.re - self.re * rhs.im) / denominator,
        )
    }
}

macro_rules! real {
    ($($real:ty),+) => {$(
        impl Conjugate for $real {
            type Output = Self;
            fn conjugate(&self) -> Self { *self }
        }
        impl RealPart for $real {
            type Real = Self;
            fn real_part(self) -> Self { self }
        }
        impl Divide for $real {
            type Output = Self;
            fn divide(self, rhs: Self) -> Self { self / rhs }
        }
    )+};
}
real!(f32, f64);

pub(crate) fn conjugate<C: Conjugate<Output = C>>(value: C) -> C {
    value.conjugate()
}
pub(crate) fn multiply<C: Mul<Output = C>>(a: C, b: C) -> C {
    a * b
}
pub(crate) fn add<C: Add<Output = C>>(a: C, b: C) -> C {
    a + b
}
pub(crate) fn zero<C: Zero>() -> C {
    C::zero()
}
