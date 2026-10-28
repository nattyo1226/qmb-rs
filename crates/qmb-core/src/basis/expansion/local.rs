//! Local basis expansion and sparse coordinate extraction.

use super::accumulate;
use crate::{
    Act, Adjoint, Basis, BasisError, Conjugate, Coordinates, Dimension, Identity, LocalOperator,
    LocalState, Product,
    expr::traversal::{LinearTerms, operator_view, state_view},
};
use alloc::{collections::BTreeMap, vec};
use core::ops::Mul;
use num_traits::Zero;

type LocalExpansion<V, C> = BTreeMap<usize, (C, V)>;
fn local_expansion<'a, P: 'a, C, B, V>(
    terms: impl Iterator<Item = (Option<C>, &'a P)>,
    basis: &B,
    dimension: usize,
    make: impl Fn(P) -> V,
) -> Result<LocalExpansion<V, C>, BasisError>
where
    B: Basis<P, C>,
    C: Clone + Zero + Mul<Output = C>,
{
    let mut result = BTreeMap::new();
    for (coefficient, primitive) in terms {
        for term in basis.expand(primitive) {
            let index = basis
                .index(&term.operand)
                .ok_or(BasisError::NonBasisElement)?;
            if index >= dimension {
                return Err(BasisError::IndexOutOfBounds { index, dimension });
            }
            let c = match &coefficient {
                None => term.coefficient,
                Some(c) => c.clone() * term.coefficient,
            };
            accumulate(&mut result, index, c, make(term.operand));
        }
    }
    Ok(result)
}

pub(super) fn local_state<S, O, C, B, V>(
    state: &LocalState<S, O, C>,
    basis: &B,
    dimension: usize,
    make: impl Fn(S) -> V,
) -> Result<LocalExpansion<V, C>, BasisError>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
    S: Act<C, Operator = O, Term = S>,
    O: Adjoint<C, Term = O>,
    B: Basis<S, C>,
{
    let state = state.evaluate();
    local_expansion(LinearTerms::new(&state, state_view), basis, dimension, make)
}

pub(super) fn local_operator<O, C, B, V>(
    op: &LocalOperator<O, C>,
    basis: &B,
    dimension: usize,
    make: impl Fn(O) -> V,
) -> Result<LocalExpansion<V, C>, BasisError>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
    O: Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
    B: Basis<O, C>,
{
    let op = op.evaluated();
    local_expansion(LinearTerms::new(&op, operator_view), basis, dimension, make)
}

fn local_coordinates<B, C>(
    basis: &B,
    dimension: usize,
    terms: LocalExpansion<(), C>,
) -> Coordinates<'_, B, C> {
    Coordinates {
        basis,
        space: (),
        factor_count: 1,
        local_dimension: Some(dimension),
        entries: terms.into_iter().map(|(i, (c, ()))| (vec![i], c)).collect(),
    }
}

impl<S, O, C> LocalState<S, O, C> {
    /// Expands and aggregates into ordinary local state syntax in `basis`.
    ///
    /// The input is unchanged. No basis metadata, dense coordinates, operator
    /// multiplication capability, or primitive cloning is required.
    pub fn expand_in<B>(&self, basis: &B) -> Result<Self, BasisError>
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        S: Dimension + Act<C, Operator = O, Term = S>,
        O: Adjoint<C, Term = O>,
        B: Basis<S, C>,
    {
        let terms = local_state(self, basis, S::dim()?.get(), Self::primitive)?;
        Ok(Self::sum(
            terms.into_values().map(|(c, state)| state.into_scaled(c)),
        ))
    }

    /// Obtains aggregated sparse basis coordinates, without building an expanded tree.
    pub fn coordinates_in<'b, B>(&self, basis: &'b B) -> Result<Coordinates<'b, B, C>, BasisError>
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        S: Dimension + Act<C, Operator = O, Term = S>,
        O: Adjoint<C, Term = O>,
        B: Basis<S, C>,
    {
        let dimension = S::dim()?.get();
        Ok(local_coordinates(
            basis,
            dimension,
            local_state(self, basis, dimension, |_| ())?,
        ))
    }
}

impl<O, C> LocalOperator<O, C> {
    /// Evaluates local products and adjoints, then expands into ordinary syntax.
    ///
    /// Equal basis indices are aggregated without comparing primitives. Exact
    /// zeros are removed; no tolerance-based pruning is performed.
    pub fn expand_in<B>(&self, basis: &B) -> Result<Self, BasisError>
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        O: Dimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
        B: Basis<O, C>,
    {
        let terms = local_operator(self, basis, O::dim()?.get(), Self::primitive)?;
        Ok(Self::sum(
            terms.into_values().map(|(c, op)| op.into_scaled(c)),
        ))
    }

    /// Obtains coordinates in an operator basis, not a matrix acting on states.
    pub fn coordinates_in<'b, B>(&self, basis: &'b B) -> Result<Coordinates<'b, B, C>, BasisError>
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        O: Dimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
        B: Basis<O, C>,
    {
        let dimension = O::dim()?.get();
        Ok(local_coordinates(
            basis,
            dimension,
            local_operator(self, basis, dimension, |_| ())?,
        ))
    }
}
