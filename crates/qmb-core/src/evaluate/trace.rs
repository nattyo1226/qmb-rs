//! Ordinary traces and Hilbert–Schmidt contractions.

use crate::{
    Adjoint, Conjugate, ExpressionError, HilbertSpace, Identity, LocalOperator, Operator, Product,
    Trace,
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
        let evaluated = self.evaluated();
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

    /// Computes `Tr(self† rhs)`, without division by the local dimension.
    pub fn hilbert_schmidt_inner_product(&self, rhs: &Self) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        O: Identity,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        self.adjoint().product(rhs).trace()
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
        let mut result = zero();
        for (coefficient, factors) in OrdinaryOpTerms::new(&self.expression, ambient) {
            let mut value = one();
            for index in ambient.indices() {
                let local_trace = match factors
                    .iter()
                    .find(|(factor, _)| factor.factor_index() == Some(index))
                {
                    Some((_, local)) => local.trace(),
                    None => LocalOperator::<O, C>::identity().trace(),
                };
                value = multiply(value, local_trace);
            }
            if let Some(c) = coefficient {
                value = multiply(c, value);
            }
            result = add(result, value);
        }
        Ok(result)
    }

    /// Computes the unnormalized Hilbert-Schmidt inner product `Tr(self† rhs)`.
    pub fn hilbert_schmidt_inner_product(&self, rhs: &Self) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Identity + crate::HomogeneousParity + crate::ParityOperator,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        self.adjoint().product(rhs)?.trace()
    }

    /// Computes Hilbert-Schmidt inner product on an explicit ambient space.
    pub fn hilbert_schmidt_inner_product_in(
        &self,
        rhs: &Self,
        ambient: &HilbertSpace<I>,
    ) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Identity + crate::HomogeneousParity + crate::ParityOperator,
        O: Product<C, Term = O> + Adjoint<C, Term = O> + Trace<C>,
    {
        self.bind(ambient)?
            .adjoint()
            .product(&rhs.bind(ambient)?)?
            .trace()
    }
}
