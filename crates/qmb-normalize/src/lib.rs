//! Norms, state normalization and normalized expectations.
//!
//! Import the corresponding extension traits to use these operations. Squared
//! norms, inner products, expectations and scalar division belong to `qmb-core`.
//! Only square-root-based operations need [`Float`] on the associated real value;
//! normalized expectations require neither square roots nor real-part extraction.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::ops::Mul;
use num_traits::{Float, One, Zero};
use qmb_core::{
    Act, Adjoint, Conjugate, Divide, ExpressionError, HilbertSpace, InnerProduct, LocalOperator,
    LocalState, LocalStructure, Operator, RealPart, State,
};

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

fn nonzero<C: Zero>(value: &C) -> Result<(), NormalizationError> {
    if value.is_zero() {
        Err(NormalizationError::ZeroNorm)
    } else {
        Ok(())
    }
}

/// A norm derived from the self-inner-product, with no clamping or zero rejection.
pub trait Norm {
    /// A real norm for local states, or a fallible real norm for composite states.
    type Output;
    /// Computes `sqrt(Re(<self | self>))`; zero has norm zero.
    fn norm(&self) -> Self::Output;
}

/// Normalization by a positive real scale, preserving the state's phase.
pub trait Normalize: Sized {
    /// Returns a shared expression scaled by its inverse norm.
    ///
    /// Rejects an exactly zero self-inner-product, including cancellation. No
    /// approximate zero test, negative-value correction or special-value policy
    /// is added.
    fn normalized(&self) -> Result<Self, NormalizationError>;
}

impl<L: LocalStructure, C> Norm for LocalState<L, C>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C> + RealPart,
    C::Real: Float,
    L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
    L::Operator: Adjoint<C, Term = L::Operator>,
{
    type Output = C::Real;
    fn norm(&self) -> Self::Output {
        self.norm_squared().sqrt()
    }
}

impl<L: LocalStructure, C> Normalize for LocalState<L, C>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C> + RealPart + From<C::Real>,
    C::Real: Float,
    L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
    L::Operator: Adjoint<C, Term = L::Operator>,
{
    fn normalized(&self) -> Result<Self, NormalizationError> {
        let denominator = self.inner_product(self);
        nonzero(&denominator)?;
        Ok(self.scaled(C::from(C::Real::one() / denominator.real_part().sqrt())))
    }
}

impl<I: Eq, L: LocalStructure, C> Norm for State<I, L, C>
where
    C: Clone + Zero + One + Conjugate<Output = C> + RealPart,
    C::Real: Float,
    L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
    L::Operator: Adjoint<C, Term = L::Operator>,
{
    type Output = Result<C::Real, ExpressionError>;
    fn norm(&self) -> Self::Output {
        Ok(self.norm_squared()?.sqrt())
    }
}

impl<I: Eq, L: LocalStructure, C> Normalize for State<I, L, C>
where
    C: Clone + Zero + One + Conjugate<Output = C> + RealPart + From<C::Real>,
    C::Real: Float,
    L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
    L::Operator: Adjoint<C, Term = L::Operator>,
{
    fn normalized(&self) -> Result<Self, NormalizationError> {
        let denominator = self.inner_product(self)?;
        nonzero(&denominator)?;
        Ok(self.scaled(C::from(C::Real::one() / denominator.real_part().sqrt())))
    }
}

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

impl<L: LocalStructure, C> NormalizedExpectation<LocalOperator<L, C>> for LocalState<L, C>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C> + Divide<Output = C>,
    L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
    L::Operator: Adjoint<C, Term = L::Operator>,
{
    type Output = C;
    fn normalized_expectation(
        &self,
        operator: &LocalOperator<L, C>,
    ) -> Result<C, NormalizationError> {
        let denominator = self.inner_product(self);
        nonzero(&denominator)?;
        Ok(self.expectation(operator).divide(denominator))
    }
}

impl<I: Eq, L: LocalStructure, C> NormalizedExpectation<Operator<I, L, C>> for State<I, L, C>
where
    C: Clone + Zero + One + Conjugate<Output = C> + Divide<Output = C>,
    L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
    L::Operator: Adjoint<C, Term = L::Operator>,
{
    type Output = C;
    fn normalized_expectation(
        &self,
        operator: &Operator<I, L, C>,
    ) -> Result<C, NormalizationError> {
        let numerator = self.expectation(operator)?;
        let denominator = self.inner_product(self)?;
        nonzero(&denominator)?;
        Ok(numerator.divide(denominator))
    }
}

impl<I: Eq, L: LocalStructure, C> NormalizedExpectationIn<Operator<I, L, C>, HilbertSpace<I, L>>
    for State<I, L, C>
where
    C: Clone + Zero + One + Conjugate<Output = C> + Divide<Output = C>,
    L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
    L::Operator: Adjoint<C, Term = L::Operator>,
{
    type Output = C;
    fn normalized_expectation_in(
        &self,
        operator: &Operator<I, L, C>,
        ambient: &HilbertSpace<I, L>,
    ) -> Result<C, NormalizationError> {
        let numerator = self.expectation_in(operator, ambient)?;
        let denominator = self.inner_product(self)?;
        nonzero(&denominator)?;
        Ok(numerator.divide(denominator))
    }
}
