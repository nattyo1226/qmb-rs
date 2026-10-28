//! Tensor-product basis expansion and sparse coordinates.

use super::{
    accumulate,
    local::{local_operator, local_state},
};
use crate::basis::enumeration::{leaves, local_size};
use crate::{
    Act, Adjoint, Basis, BasisError, Conjugate, Coordinates, Dimension, HilbertSpace, Identity,
    LocalOperator, LocalState, Operator, Product, State,
    expr::{
        composite::{OpExpr, OpNode, StateExpr, StateNode},
        parts::{ProductExpr, ScaleExpr, SumExpr},
        traversal::OpTerms,
    },
};
use alloc::{collections::BTreeMap, vec, vec::Vec};
use core::ops::Mul;
use num_traits::{One, Zero};

type Expansion<V, C> = BTreeMap<Vec<usize>, (C, V)>;

// Coordinates carry no primitive payload. Expansion carries only cheap local
// expression handles, so neither path needs Clone/Eq/Ord on primitives.
trait TensorValue: Clone {
    fn unit() -> Self;
    fn tensor(&self, rhs: &Self) -> Self;
}
impl TensorValue for () {
    fn unit() -> Self {}
    fn tensor(&self, _: &Self) -> Self {}
}
impl<E: Clone> TensorValue for Vec<E> {
    fn unit() -> Self {
        Vec::new()
    }
    fn tensor(&self, rhs: &Self) -> Self {
        self.iter().chain(rhs).cloned().collect()
    }
}

fn tensor<V: TensorValue, C: Clone + Zero + Mul<Output = C>>(
    lhs: Expansion<V, C>,
    rhs: Expansion<V, C>,
) -> Expansion<V, C> {
    let mut result = BTreeMap::new();
    for (a, (c, x)) in lhs {
        for (b, (d, y)) in &rhs {
            let mut digits = a.clone();
            digits.extend(b);
            accumulate(&mut result, digits, c.clone() * d.clone(), x.tensor(y));
        }
    }
    result
}

enum StateTask<'a, I, S, O, C> {
    Visit(&'a StateExpr<I, S, O, C>),
    Scale(&'a C),
    Sum(usize),
    Tensor(usize),
}

fn state_expansion<I: Eq, S, O, C, B, V: TensorValue>(
    state: &State<I, S, O, C>,
    basis: &B,
    dimension: Option<usize>,
    make: impl Fn(S) -> V,
) -> Result<Expansion<V, C>, BasisError>
where
    C: Clone + Zero + One + Conjugate<Output = C>,
    S: Act<C, Operator = O, Term = S>,
    O: Adjoint<C, Term = O>,
    B: Basis<S, C>,
{
    let state = state.evaluate()?;
    let mut tasks = vec![StateTask::Visit(&state.expression)];
    let mut values: Vec<Expansion<V, C>> = Vec::new();
    while let Some(task) = tasks.pop() {
        match task {
            StateTask::Visit(node) => match node.node() {
                StateNode::Local(_, local) => {
                    let terms = local_state(local, basis, dimension.expect("local factor"), &make)?;
                    values.push(terms.into_iter().map(|(i, v)| (vec![i], v)).collect());
                }
                StateNode::Unit => {
                    values.push(BTreeMap::from([(Vec::new(), (C::one(), V::unit()))]))
                }
                StateNode::Scale(s) => {
                    tasks.push(StateTask::Scale(&s.coefficient));
                    tasks.push(StateTask::Visit(&s.operand));
                }
                StateNode::Sum(s) | StateNode::Tensor(s) => {
                    tasks.push(if matches!(node.node(), StateNode::Sum(_)) {
                        StateTask::Sum(s.len())
                    } else {
                        StateTask::Tensor(s.len())
                    });
                    tasks.extend(s.iter().rev().map(StateTask::Visit));
                }
                StateNode::Act(_) => unreachable!("evaluated state"),
            },
            StateTask::Scale(c) => {
                let value = values.last_mut().expect("scale operand");
                for (d, _) in value.values_mut() {
                    *d = c.clone() * core::mem::replace(d, C::zero());
                }
                value.retain(|_, (d, _)| !d.is_zero());
            }
            StateTask::Sum(count) => {
                let mut sum = BTreeMap::new();
                for term in values.drain(values.len() - count..) {
                    for (digits, (c, v)) in term {
                        accumulate(&mut sum, digits, c, v);
                    }
                }
                values.push(sum);
            }
            StateTask::Tensor(count) => {
                let factors = values.split_off(values.len() - count);
                let mut product = BTreeMap::from([(Vec::new(), (C::one(), V::unit()))]);
                for factor in factors {
                    product = tensor(product, factor);
                }
                values.push(product);
            }
        }
    }
    let mut result = values.pop().expect("one expansion");
    result.retain(|_, (c, _)| !c.is_zero());
    Ok(result)
}

fn operator_expansion<I: Eq, O, C, B, V: TensorValue>(
    op: &Operator<I, O, C>,
    basis: &B,
    factors: &[&HilbertSpace<I>],
    dimension: Option<usize>,
    make: impl Fn(O) -> V,
) -> Result<Expansion<V, C>, BasisError>
where
    C: Clone + Zero + One + Conjugate<Output = C>,
    O: Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
    B: Basis<O, C>,
{
    let mut result = BTreeMap::new();
    for (coefficient, locals) in OpTerms::new(&op.expression) {
        let mut product =
            BTreeMap::from([(Vec::new(), (coefficient.unwrap_or_else(C::one), V::unit()))]);
        for factor in factors {
            let identity;
            let local = match locals
                .iter()
                .find(|(f, _)| f.factor_index() == factor.factor_index())
            {
                Some((_, local)) => local,
                None => {
                    identity = LocalOperator::identity();
                    &identity
                }
            };
            let terms = local_operator(local, basis, dimension.expect("operator factor"), &make)?;
            product = tensor(
                product,
                terms.into_iter().map(|(i, v)| (vec![i], v)).collect(),
            );
        }
        for (digits, (c, v)) in product {
            accumulate(&mut result, digits, c, v);
        }
    }
    Ok(result)
}

impl<I: Eq, S, O, C> State<I, S, O, C> {
    /// Expands in one local basis applied uniformly in the space's factor order.
    ///
    /// Returns ordinary state syntax on the original shared space. No global
    /// dimension, primitive clone, or primitive operator product is required.
    pub fn expand_in<B>(&self, basis: &B) -> Result<Self, BasisError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        S: Dimension + Act<C, Operator = O, Term = S>,
        O: Adjoint<C, Term = O>,
        B: Basis<S, C>,
    {
        let factors = leaves(&self.space);
        let dimension = local_size::<S>(factors.len())?;
        let terms = state_expansion(self, basis, dimension, |p| vec![LocalState::primitive(p)])?;
        let expressions = terms
            .into_values()
            .map(|(c, locals)| {
                let nodes: SumExpr<_> = locals
                    .into_iter()
                    .zip(&factors)
                    .map(|(local, factor)| {
                        StateExpr::new(StateNode::Local((*factor).clone(), local))
                    })
                    .collect();
                let term = match nodes {
                    SumExpr::Empty => StateExpr::new(StateNode::Unit),
                    SumExpr::One(node) => node,
                    nodes => StateExpr::new(StateNode::Tensor(nodes)),
                };
                term.scaled(c)
            })
            .collect();
        Ok(Self {
            space: self.space.clone(),
            expression: StateExpr::sum(expressions),
        })
    }

    /// Obtains sparse coordinates with the original space and borrowed local basis.
    pub fn coordinates_in<'b, B>(
        &self,
        basis: &'b B,
    ) -> Result<Coordinates<'b, B, C, HilbertSpace<I>>, BasisError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        S: Dimension + Act<C, Operator = O, Term = S>,
        O: Adjoint<C, Term = O>,
        B: Basis<S, C>,
    {
        let count = self.space.indices().count();
        let dimension = local_size::<S>(count)?;
        let terms = state_expansion(self, basis, dimension, |_| ())?;
        Ok(Coordinates {
            basis,
            space: self.space.clone(),
            factor_count: count,
            local_dimension: dimension,
            entries: terms.into_iter().map(|(k, (c, ()))| (k, c)).collect(),
        })
    }
}

impl<I: Eq, O, C> Operator<I, O, C> {
    /// Expands on the existing binding; an unbound operator returns an error.
    pub fn expand_in<B>(&self, basis: &B) -> Result<Self, BasisError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Dimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
        B: Basis<O, C>,
    {
        self.expand_in_space(basis, self.select_space(None, None)?)
    }

    /// Expands on an explicit space, overriding any previous binding.
    ///
    /// Spectators are expanded as local identities in the same basis, even when
    /// identity itself is not a basis primitive. Zero retains the selected space.
    pub fn expand_in_space<B>(
        &self,
        basis: &B,
        ambient: &HilbertSpace<I>,
    ) -> Result<Self, BasisError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Dimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
        B: Basis<O, C>,
    {
        let ambient = self.select_space(Some(ambient), None)?;
        let factors = leaves(ambient);
        let dimension = local_size::<O>(factors.len())?;
        let terms = operator_expansion(self, basis, &factors, dimension, |p| {
            vec![LocalOperator::primitive(p)]
        })?;
        let expressions = terms
            .into_values()
            .map(|(c, locals)| {
                let mut iter = locals
                    .into_iter()
                    .zip(&factors)
                    .map(|(local, factor)| OpExpr::new(OpNode::Local((*factor).clone(), local)));
                let product = match iter.next() {
                    None => OpExpr::new(OpNode::Identity),
                    Some(first) => iter.fold(first, |lhs, rhs| {
                        OpExpr::new(OpNode::Product(ProductExpr { lhs, rhs }))
                    }),
                };
                OpExpr::new(OpNode::Scale(ScaleExpr {
                    operand: product,
                    coefficient: c,
                }))
            })
            .collect();
        Ok(Self {
            ambient: Some(ambient.clone()),
            expression: OpExpr::new(OpNode::Sum(expressions)),
        })
    }

    /// Obtains operator-basis coordinates using the existing binding.
    pub fn coordinates_in<'b, B>(
        &self,
        basis: &'b B,
    ) -> Result<Coordinates<'b, B, C, HilbertSpace<I>>, BasisError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Dimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
        B: Basis<O, C>,
    {
        self.coordinates_in_space(basis, self.select_space(None, None)?)
    }

    /// Obtains operator-basis coordinates on an explicit space, overriding a binding.
    pub fn coordinates_in_space<'b, B>(
        &self,
        basis: &'b B,
        ambient: &HilbertSpace<I>,
    ) -> Result<Coordinates<'b, B, C, HilbertSpace<I>>, BasisError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        O: Dimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
        B: Basis<O, C>,
    {
        let ambient = self.select_space(Some(ambient), None)?;
        let factors = leaves(ambient);
        let dimension = local_size::<O>(factors.len())?;
        let terms = operator_expansion(self, basis, &factors, dimension, |_| ())?;
        Ok(Coordinates {
            basis,
            space: ambient.clone(),
            factor_count: factors.len(),
            local_dimension: dimension,
            entries: terms.into_iter().map(|(k, (c, ()))| (k, c)).collect(),
        })
    }
}
