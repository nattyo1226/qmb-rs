//! Norms and square-root-based state normalization.

use crate::{NormalizationError, error::nonzero};
use core::ops::Mul;
use num_traits::{Float, One, Zero};
use qmb_core::{
    Act, Adjoint, Conjugate, ExpressionError, InnerProduct, LocalState, RealPart, State,
};

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

impl<S, O, C> Norm for LocalState<S, O, C>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C> + RealPart,
    C::Real: Float,
    S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
    O: Adjoint<C, Term = O>,
{
    type Output = C::Real;
    fn norm(&self) -> Self::Output {
        self.norm_squared().sqrt()
    }
}

impl<S, O, C> Normalize for LocalState<S, O, C>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C> + RealPart + From<C::Real>,
    C::Real: Float,
    S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
    O: Adjoint<C, Term = O>,
{
    fn normalized(&self) -> Result<Self, NormalizationError> {
        let denominator = self.inner_product(self);
        nonzero(&denominator)?;
        Ok(self.scaled(C::from(C::Real::one() / denominator.real_part().sqrt())))
    }
}

impl<I: Eq, S, O, C> Norm for State<I, S, O, C>
where
    C: Clone + Zero + One + Conjugate<Output = C> + RealPart,
    C::Real: Float,
    S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
    O: Adjoint<C, Term = O> + qmb_core::HomogeneousParity + qmb_core::ParityOperator,
{
    type Output = Result<C::Real, ExpressionError>;
    fn norm(&self) -> Self::Output {
        Ok(self.norm_squared()?.sqrt())
    }
}

impl<I: Eq, S, O, C> Normalize for State<I, S, O, C>
where
    C: Clone + Zero + One + Conjugate<Output = C> + RealPart + From<C::Real>,
    C::Real: Float,
    S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
    O: Adjoint<C, Term = O> + qmb_core::HomogeneousParity + qmb_core::ParityOperator,
{
    fn normalized(&self) -> Result<Self, NormalizationError> {
        let denominator = self.inner_product(self)?;
        nonzero(&denominator)?;
        Ok(self.scaled(C::from(C::Real::one() / denominator.real_part().sqrt())))
    }
}
