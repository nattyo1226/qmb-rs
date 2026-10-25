//! Derived state operations using contraction, without a basis or matrix.

use super::ExpressionError;
use crate::{
    Act, Adjoint, Complex, HilbertSpace, InnerProduct, LocalOperator, LocalState, LocalStructure,
    Operator, State,
};
use core::ops::{Div, Mul, Neg, Sub};
use num_traits::{One, Zero};

/// Failure to normalize a state or an expectation value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum NormalizationError {
    /// The computed self-inner-product is exactly zero, including cancellation.
    #[error("state has zero norm")]
    ZeroNorm,
    /// Placement, binding or contraction failed.
    #[error(transparent)]
    Expression(#[from] ExpressionError),
}

fn nonzero<T: Zero>(value: &Complex<T>) -> Result<(), NormalizationError> {
    if value.re.is_zero() && value.im.is_zero() {
        Err(NormalizationError::ZeroNorm)
    } else {
        Ok(())
    }
}

// Avoid imposing Num (and remainder) on square-root-free derived operations.
fn quotient<T>(a: Complex<T>, b: Complex<T>) -> Complex<T>
where
    T: Clone + Zero + Sub<Output = T> + Mul<Output = T> + Div<Output = T>,
{
    // Valid self-inner-products are real. Keeping this path also avoids squaring
    // a very large/small real denominator. The other path retains complex division.
    if b.im.is_zero() {
        return Complex::new(a.re / b.re.clone(), a.im / b.re);
    }
    let denominator = b.re.clone() * b.re.clone() + b.im.clone() * b.im.clone();
    Complex::new(
        (a.re.clone() * b.re.clone() + a.im.clone() * b.im.clone()) / denominator.clone(),
        (a.im * b.re - a.re * b.im) / denominator,
    )
}

impl<L: LocalStructure, T> LocalState<L, T>
where
    T: Clone + Zero + Sub<Output = T> + Mul<Output = T> + Neg<Output = T>,
    L::State: Act<T, Operator = L::Operator, Term = L::State> + InnerProduct<T>,
    L::Operator: Adjoint<T, Term = L::Operator>,
{
    /// Returns the real part of the self-inner-product, without clamping.
    pub fn norm_squared(&self) -> T {
        self.inner_product(self).re
    }

    /// Divides the expectation by the full complex self-inner-product.
    pub fn normalized_expectation(
        &self,
        operator: &LocalOperator<L, T>,
    ) -> Result<Complex<T>, NormalizationError>
    where
        T: Div<Output = T>,
    {
        let denominator = self.inner_product(self);
        nonzero(&denominator)?;
        Ok(quotient(self.expectation(operator), denominator))
    }

    /// Returns the square root of the norm squared. Zero has norm zero.
    #[cfg(feature = "libm")]
    pub fn norm(&self) -> T
    where
        T: num_traits::Float,
    {
        self.norm_squared().sqrt()
    }

    /// Returns a shared expression scaled by the inverse norm.
    #[cfg(feature = "libm")]
    pub fn normalized(&self) -> Result<Self, NormalizationError>
    where
        T: num_traits::Float,
    {
        let denominator = self.inner_product(self);
        nonzero(&denominator)?;
        Ok(self.scaled(Complex::new(T::one() / denominator.re.sqrt(), T::zero())))
    }
}

impl<I: Eq, L: LocalStructure, T> State<I, L, T>
where
    T: Clone + Zero + One + Sub<Output = T> + Neg<Output = T>,
    L::State: Act<T, Operator = L::Operator, Term = L::State> + InnerProduct<T>,
    L::Operator: Adjoint<T, Term = L::Operator>,
{
    /// Returns the real part of the self-inner-product, without clamping.
    pub fn norm_squared(&self) -> Result<T, ExpressionError> {
        Ok(self.inner_product(self)?.re)
    }

    /// Divides the expectation by the full complex self-inner-product.
    pub fn normalized_expectation(
        &self,
        operator: &Operator<I, L, T>,
    ) -> Result<Complex<T>, NormalizationError>
    where
        T: Div<Output = T>,
    {
        let numerator = self.expectation(operator)?;
        let denominator = self.inner_product(self)?;
        nonzero(&denominator)?;
        Ok(quotient(numerator, denominator))
    }

    /// Computes a normalized expectation with explicit-space precedence.
    pub fn normalized_expectation_in(
        &self,
        operator: &Operator<I, L, T>,
        ambient: &HilbertSpace<I, L>,
    ) -> Result<Complex<T>, NormalizationError>
    where
        T: Div<Output = T>,
    {
        let numerator = self.expectation_in(operator, ambient)?;
        let denominator = self.inner_product(self)?;
        nonzero(&denominator)?;
        Ok(quotient(numerator, denominator))
    }

    /// Returns the square root of the norm squared. Zero has norm zero.
    #[cfg(feature = "libm")]
    pub fn norm(&self) -> Result<T, ExpressionError>
    where
        T: num_traits::Float,
    {
        Ok(self.norm_squared()?.sqrt())
    }

    /// Returns a shared expression scaled by the inverse norm.
    #[cfg(feature = "libm")]
    pub fn normalized(&self) -> Result<Self, NormalizationError>
    where
        T: num_traits::Float,
    {
        let denominator = self.inner_product(self)?;
        nonzero(&denominator)?;
        Ok(self.scaled(Complex::new(T::one() / denominator.re.sqrt(), T::zero())))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{LinearCombination, LocalStructure};
    use core::{num::NonZeroUsize, ops::Add};

    // No Num, Float, Copy, PartialEq, or remainder operation. Derived algebraic
    // operations must not acquire a square-root backend or broad numeric bound.
    #[derive(Clone)]
    struct Scalar(f64);
    macro_rules! binary {
        ($trait:ident, $method:ident, $op:tt) => {
            impl $trait for Scalar {
                type Output=Self;
                fn $method(self,rhs:Self)->Self { Self(self.0 $op rhs.0) }
            }
        }
    }
    binary!(Add,add,+);
    binary!(Sub,sub,-);
    binary!(Mul,mul,*);
    binary!(Div,div,/);
    impl Neg for Scalar {
        type Output = Self;
        fn neg(self) -> Self {
            Self(-self.0)
        }
    }
    impl Zero for Scalar {
        fn zero() -> Self {
            Self(0.0)
        }
        fn is_zero(&self) -> bool {
            self.0 == 0.0
        }
    }
    impl One for Scalar {
        fn one() -> Self {
            Self(1.0)
        }
    }
    struct Primitive;
    struct Structure;
    impl LocalStructure for Structure {
        type State = Primitive;
        type Operator = Primitive;
        const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
    }
    fn term() -> LinearCombination<Primitive, Scalar> {
        LinearCombination::term(Complex::new(Scalar(1.0), Scalar(0.0)), Primitive)
    }
    impl Act<Scalar> for Primitive {
        type Operator = Self;
        type Term = Self;
        fn act(&self, _: &Self) -> LinearCombination<Self, Scalar> {
            term()
        }
    }
    impl Adjoint<Scalar> for Primitive {
        type Term = Self;
        fn adjoint(&self) -> LinearCombination<Self, Scalar> {
            term()
        }
    }
    impl InnerProduct<Scalar> for Primitive {
        fn inner_product(&self, _: &Self) -> Complex<Scalar> {
            Complex::new(Scalar(1.0), Scalar(0.0))
        }
    }
    #[test]
    fn algebraic_norms_and_expectations_require_no_float_or_num() {
        let state = LocalState::<Structure, Scalar>::primitive(Primitive)
            .scaled(Complex::new(Scalar(2.0), Scalar(1.0)));
        let op = LocalOperator::primitive(Primitive);
        assert_eq!(state.norm_squared().0, 5.0);
        let value = state.normalized_expectation(&op).unwrap();
        assert_eq!((value.re.0, value.im.0), (1.0, 0.0));
        let factor = HilbertSpace::factor(7);
        let state = state.at(&factor).unwrap();
        let op = op.at(&factor).unwrap();
        assert_eq!(state.norm_squared().unwrap().0, 5.0);
        let value = state.normalized_expectation(&op).unwrap();
        assert_eq!((value.re.0, value.im.0), (1.0, 0.0));
    }
    #[test]
    fn quotient_retains_the_imaginary_part_of_the_denominator() {
        assert_eq!(
            quotient(Complex::new(3.0, 4.0), Complex::new(1.0, 2.0)),
            Complex::new(2.2, -0.4)
        );
    }
}
