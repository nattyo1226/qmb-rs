//! Ordinary traces and Hilbert–Schmidt contractions.

use crate::{
    Adjoint, Conjugate, Divide, ExpressionError, HilbertSpace, Identity, LocalOperator, Operator,
    Product, Trace,
    coefficient::{add, multiply, zero},
    expr::traversal::{LinearTerms, OrdinaryOpTerms, operator_view},
};
use core::ops::Mul;
use num_traits::{One, Zero};

impl<O, C> LocalOperator<O, C> {
    /// Computes the ordinary, unnormalized local trace without a chosen basis.
    pub fn trace(&self) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        O: Identity,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        let evaluated = self.evaluate();
        let mut result = zero();
        for (coefficient, primitive) in LinearTerms::new(&evaluated, operator_view) {
            let mut value = primitive.trace();
            if let Some(c) = coefficient {
                value = multiply(c, value);
            }
            result = add(result, value);
        }
        result
    }

    /// Computes the normalized Hilbert–Schmidt inner product `Tr(self† rhs)/d`.
    ///
    /// The local dimension is obtained as `Tr(I)` in the coefficient type. It
    /// must be nonzero and division by it must represent mathematical division.
    pub fn hs_inner_product(&self, rhs: &Self) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C> + Divide<Output = C>,
        O: Identity,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        (self.adjoint() * rhs).trace().divide(O::identity().trace())
    }
}

fn one<C: One>() -> C {
    C::one()
}

impl<I: Eq, O, C> Operator<I, O, C> {
    /// Computes the ordinary full trace using the existing binding.
    pub fn trace(&self) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Identity + crate::HomogeneousParity + crate::ParityOperator,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        let ambient = self.select_space(None, None)?;
        self.trace_on(ambient)
    }

    /// Computes full trace on an explicit space, overriding an old binding.
    pub fn trace_in(&self, ambient: &HilbertSpace<I>) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Identity + crate::HomogeneousParity + crate::ParityOperator,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        let ambient = self.select_space(Some(ambient), None)?;
        self.trace_on(ambient)
    }

    fn trace_on(&self, ambient: &HilbertSpace<I>) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Identity + crate::HomogeneousParity + crate::ParityOperator,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        self.trace_on_with(ambient, |value| value)
    }

    // Both contractions use ordinary local traces after grading has been
    // resolved. Only the HS caller supplies division; ordinary trace keeps its
    // original coefficient bounds.
    fn trace_on_with(
        &self,
        ambient: &HilbertSpace<I>,
        local_trace: impl Fn(C) -> C,
    ) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Identity + crate::HomogeneousParity + crate::ParityOperator,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        let mut result = zero();
        for (coefficient, factors) in OrdinaryOpTerms::new(&self.expression, ambient) {
            let mut value = one();
            for index in ambient.indices() {
                let local_trace_value = match factors
                    .iter()
                    .find(|(factor, _)| factor.factor_index() == Some(index))
                {
                    Some((_, local)) => local.trace(),
                    None => LocalOperator::<O, C>::identity().trace(),
                };
                value = multiply(value, local_trace(local_trace_value));
            }
            if let Some(c) = coefficient {
                value = multiply(c, value);
            }
            result = add(result, value);
        }
        Ok(result)
    }

    fn normalized_trace_on(&self, ambient: &HilbertSpace<I>) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C> + Divide<Output = C>,
        O: Identity + crate::HomogeneousParity + crate::ParityOperator,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        let dimension = O::identity().trace();
        self.trace_on_with(ambient, |value| value.divide(dimension.clone()))
    }

    /// Computes the normalized Hilbert–Schmidt inner product `Tr(self† rhs)/D`.
    ///
    /// Each ordinary local trace is divided by `Tr(I)` before contraction, so
    /// no integer total dimension or large unnormalized trace is constructed.
    /// Local dimensions must be nonzero and support mathematical division.
    pub fn hs_inner_product(&self, rhs: &Self) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C> + Divide<Output = C>,
        O: Identity + crate::HomogeneousParity + crate::ParityOperator,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        let product = (self.adjoint() * rhs)?;
        let ambient = product.select_space(None, None)?;
        product.normalized_trace_on(ambient)
    }

    /// Computes `Tr(self† rhs)/D` on an explicit ambient space.
    ///
    /// Uses local trace normalization, with the same coefficient requirements
    /// as [`Self::hs_inner_product`]. The explicit space overrides old bindings.
    pub fn hs_inner_product_in(
        &self,
        rhs: &Self,
        ambient: &HilbertSpace<I>,
    ) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C> + Divide<Output = C>,
        O: Identity + crate::HomogeneousParity + crate::ParityOperator,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        (self.bind(ambient)?.adjoint() * &rhs.bind(ambient)?)?.normalized_trace_on(ambient)
    }
}
