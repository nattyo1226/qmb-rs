//! Basis expansion and sparse coordinates share aggregation, not result trees.

use alloc::{collections::BTreeMap, vec, vec::Vec};
use core::ops::Mul;
use num_traits::{One, Zero};

use super::{
    composite::{OpExpr, OpNode, Operator, State, StateExpr, StateNode},
    composite_eval::OpTerms,
    evaluate::{LinearTerms, operator_view, state_view},
    local::{LocalOperator, LocalState},
    parts::{ScaleExpr, SumExpr},
};
use crate::{
    Act, Adjoint, Basis, BasisError, Conjugate, Coordinates, HilbertSpace, Identity, Product,
    SpanDimension,
};

type LocalExpansion<V, C> = BTreeMap<usize, (C, V)>;
type Expansion<V, C> = BTreeMap<Vec<usize>, (C, V)>;
type BasisElement<T> = Result<T, BasisError>;

// Remove exact zeros before intermediate maps feed later tensor products.
fn accumulate<K: Ord, V, C: Zero>(map: &mut BTreeMap<K, (C, V)>, key: K, c: C, value: V) {
    use alloc::collections::btree_map::Entry;
    match map.entry(key) {
        Entry::Vacant(entry) => {
            if !c.is_zero() {
                entry.insert((c, value));
            }
        }
        Entry::Occupied(mut entry) => {
            let existing = &mut entry.get_mut().0;
            *existing = core::mem::replace(existing, C::zero()) + c;
            if existing.is_zero() {
                entry.remove();
            }
        }
    }
}

fn local_expansion<'a, P: 'a, C, B, V>(
    terms: impl Iterator<Item = (Option<C>, Option<&'a P>)>,
    basis: &B,
    dimension: usize,
    identity: impl Fn() -> P,
    make: impl Fn(P) -> V,
) -> Result<LocalExpansion<V, C>, BasisError>
where
    B: Basis<P, C>,
    C: Clone + Zero + Mul<Output = C>,
{
    let mut result = BTreeMap::new();
    for (coefficient, primitive) in terms {
        let generated;
        let primitive = match primitive {
            Some(p) => p,
            None => {
                generated = identity();
                &generated
            }
        };
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

fn local_state<S, O, C, B, V>(
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
    local_expansion(
        LinearTerms::new(&state, state_view),
        basis,
        dimension,
        || unreachable!("state primitive"),
        make,
    )
}

fn local_operator<O, C, B, V>(
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
    local_expansion(
        LinearTerms::new(&op, operator_view),
        basis,
        dimension,
        O::identity,
        make,
    )
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
        S: SpanDimension + Act<C, Operator = O, Term = S>,
        O: Adjoint<C, Term = O>,
        B: Basis<S, C>,
    {
        let terms = local_state(self, basis, S::span_dimension()?.get(), Self::primitive)?;
        Ok(Self::sum(
            terms.into_values().map(|(c, state)| state.into_scaled(c)),
        ))
    }

    /// Obtains aggregated sparse basis coordinates, without building an expanded tree.
    pub fn coordinates_in<'b, B>(&self, basis: &'b B) -> Result<Coordinates<'b, B, C>, BasisError>
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        S: SpanDimension + Act<C, Operator = O, Term = S>,
        O: Adjoint<C, Term = O>,
        B: Basis<S, C>,
    {
        let dimension = S::span_dimension()?.get();
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
        O: SpanDimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
        B: Basis<O, C>,
    {
        let terms = local_operator(self, basis, O::span_dimension()?.get(), Self::primitive)?;
        Ok(Self::sum(
            terms.into_values().map(|(c, op)| op.into_scaled(c)),
        ))
    }

    /// Obtains coordinates in an operator basis, not a matrix acting on states.
    pub fn coordinates_in<'b, B>(&self, basis: &'b B) -> Result<Coordinates<'b, B, C>, BasisError>
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        O: SpanDimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
        B: Basis<O, C>,
    {
        let dimension = O::span_dimension()?.get();
        Ok(local_coordinates(
            basis,
            dimension,
            local_operator(self, basis, dimension, |_| ())?,
        ))
    }
}

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

fn leaves<I>(space: &HilbertSpace<I>) -> Vec<&HilbertSpace<I>> {
    let mut pending = vec![space];
    let mut result = Vec::new();
    while let Some(space) = pending.pop() {
        if space.factor_index().is_some() {
            result.push(space);
        } else {
            pending.extend(space.children().rev());
        }
    }
    result
}

fn local_size<P: SpanDimension>(count: usize) -> Result<Option<usize>, BasisError> {
    if count == 0 {
        Ok(None)
    } else {
        Ok(Some(P::span_dimension()?.get()))
    }
}

// An odometer, independent of the representability of the global dimension.
fn basis_indices(count: usize, dimension: Option<usize>) -> impl Iterator<Item = Vec<usize>> {
    let mut digits = vec![0; count];
    let mut done = false;
    core::iter::from_fn(move || {
        if done {
            return None;
        }
        let current = digits.clone();
        done = true;
        for digit in digits.iter_mut().rev() {
            *digit += 1;
            if *digit < dimension.expect("nonempty basis index") {
                done = false;
                break;
            }
            *digit = 0;
        }
        Some(current)
    })
}

fn basis_element<P, C, B: Basis<P, C>>(
    basis: &B,
    index: usize,
    dimension: usize,
) -> Result<P, BasisError> {
    let primitive = basis
        .elements()
        .nth(index)
        .ok_or(BasisError::InvalidEnumeration)?;
    let actual = basis.index(&primitive).ok_or(BasisError::NonBasisElement)?;
    if actual >= dimension {
        return Err(BasisError::IndexOutOfBounds {
            index: actual,
            dimension,
        });
    }
    if actual != index || (index == dimension - 1 && basis.elements().nth(dimension).is_some()) {
        return Err(BasisError::InvalidEnumeration);
    }
    Ok(primitive)
}

impl<I> HilbertSpace<I> {
    /// Lazily enumerates tensor-product basis states in mixed-radix order.
    ///
    /// The same local basis is used on every factor. Enumeration needs neither
    /// index comparison, algebraic evaluation, nor cloning primitives/indices,
    /// and does not require a representable global dimension. The unit space
    /// yields one unit state.
    /// Invalid enumeration is reported when the corresponding element is read.
    pub fn basis_states<S, O, C, B>(
        &self,
        basis: &B,
    ) -> Result<impl Iterator<Item = BasisElement<State<I, S, O, C>>>, BasisError>
    where
        S: SpanDimension,
        B: Basis<S, C>,
    {
        let factors = leaves(self);
        let dimension = local_size::<S>(factors.len())?;
        Ok(basis_indices(factors.len(), dimension).map(move |digits| {
            let nodes = digits
                .into_iter()
                .zip(&factors)
                .map(|(index, factor)| {
                    let primitive = basis_element(basis, index, dimension.expect("basis factor"))?;
                    Ok(StateExpr::new(StateNode::Local(
                        (*factor).clone(),
                        LocalState::primitive(primitive),
                    )))
                })
                .collect::<Result<SumExpr<_>, BasisError>>()?;
            let expression = match nodes {
                SumExpr::Empty => StateExpr::new(StateNode::Unit),
                SumExpr::One(node) => node,
                nodes => StateExpr::new(StateNode::Tensor(nodes)),
            };
            Ok(State {
                space: self.clone(),
                expression,
            })
        }))
    }

    /// Lazily enumerates operator-basis products, bound to this space.
    ///
    /// These are basis operators, not columns of an action matrix. No index
    /// comparison, product, adjoint, identity, coefficient arithmetic, or
    /// primitive clone is needed.
    pub fn basis_operators<O, C, B>(
        &self,
        basis: &B,
    ) -> Result<impl Iterator<Item = Result<Operator<I, O, C>, BasisError>>, BasisError>
    where
        O: SpanDimension,
        B: Basis<O, C>,
    {
        let factors = leaves(self);
        let dimension = local_size::<O>(factors.len())?;
        Ok(basis_indices(factors.len(), dimension).map(move |digits| {
            let mut nodes = Vec::new();
            for (index, factor) in digits.into_iter().zip(&factors) {
                let primitive = basis_element(basis, index, dimension.expect("basis factor"))?;
                nodes.push(OpExpr::new(OpNode::Local(
                    (*factor).clone(),
                    LocalOperator::primitive(primitive),
                )));
            }
            let mut nodes = nodes.into_iter();
            let expression = match nodes.next() {
                None => OpExpr::new(OpNode::Identity),
                Some(first) => nodes.fold(first, |lhs, rhs| {
                    OpExpr::new(OpNode::Product(super::parts::ProductExpr { lhs, rhs }))
                }),
            };
            Ok(Operator {
                ambient: Some(self.clone()),
                expression,
            })
        }))
    }
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
        S: SpanDimension + Act<C, Operator = O, Term = S>,
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
        S: SpanDimension + Act<C, Operator = O, Term = S>,
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
        O: SpanDimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
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
        O: SpanDimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
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
                        OpExpr::new(OpNode::Product(super::parts::ProductExpr { lhs, rhs }))
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
        O: SpanDimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
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
        O: SpanDimension + Identity + Product<C, Term = O> + Adjoint<C, Term = O>,
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
