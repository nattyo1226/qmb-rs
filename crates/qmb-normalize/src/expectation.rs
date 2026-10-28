//! Expectations normalized by the full self-inner-product.

use crate::{NormalizationError, error::nonzero};
use core::ops::Mul;
use num_traits::{One, Zero};
use qmb_core::{
    Act, Adjoint, Conjugate, Divide, HilbertSpace, InnerProduct, LocalOperator, LocalState,
    Operator, State,
};

/// An expectation divided by the full self-inner-product.
///
/// This operation rejects an exactly zero self-inner-product. It does not need
/// a square root, a real projection, or real-to-coefficient conversion.
pub trait NormalizedExpectation<O> {
    /// Coefficient returned by the contraction.
    type Output;
    /// Computes `<self | operator | self> / <self | self>`.
    fn normalized_expectation(&self, operator: &O) -> Result<Self::Output, NormalizationError>;
}

/// A normalized expectation with explicit ambient-space precedence.
pub trait NormalizedExpectationIn<O, H> {
    /// Coefficient returned by the contraction.
    type Output;
    /// Uses `ambient` in preference to the operator's stored binding.
    fn normalized_expectation_in(
        &self,
        operator: &O,
        ambient: &H,
    ) -> Result<Self::Output, NormalizationError>;
}

impl<S, O, C> NormalizedExpectation<LocalOperator<O, C>> for LocalState<S, O, C>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C> + Divide<Output = C>,
    S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
    O: Adjoint<C, Term = O>,
{
    type Output = C;
    fn normalized_expectation(
        &self,
        operator: &LocalOperator<O, C>,
    ) -> Result<C, NormalizationError> {
        let denominator = self.inner_product(self);
        nonzero(&denominator)?;
        Ok(self.expectation(operator).divide(denominator))
    }
}

impl<I: Eq, S, O, C> NormalizedExpectation<Operator<I, O, C>> for State<I, S, O, C>
where
    C: Clone + Zero + One + Conjugate<Output = C> + Divide<Output = C>,
    S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
    O: Adjoint<C, Term = O>,
{
    type Output = C;
    fn normalized_expectation(
        &self,
        operator: &Operator<I, O, C>,
    ) -> Result<C, NormalizationError> {
        let numerator = self.expectation(operator)?;
        let denominator = self.inner_product(self)?;
        nonzero(&denominator)?;
        Ok(numerator.divide(denominator))
    }
}

impl<I: Eq, S, O, C> NormalizedExpectationIn<Operator<I, O, C>, HilbertSpace<I>>
    for State<I, S, O, C>
where
    C: Clone + Zero + One + Conjugate<Output = C> + Divide<Output = C>,
    S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
    O: Adjoint<C, Term = O>,
{
    type Output = C;
    fn normalized_expectation_in(
        &self,
        operator: &Operator<I, O, C>,
        ambient: &HilbertSpace<I>,
    ) -> Result<C, NormalizationError> {
        let numerator = self.expectation_in(operator, ambient)?;
        let denominator = self.inner_product(self)?;
        nonzero(&denominator)?;
        Ok(numerator.divide(denominator))
    }
}
