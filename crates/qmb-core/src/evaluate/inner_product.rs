//! Contraction of evaluated states without distributing independent factors.

use alloc::{collections::BTreeMap, vec, vec::Vec};
use core::{ops::Mul, ptr};
use num_traits::{One, Zero};

use crate::expr::{
    composite::{StateExpr, StateNode},
    traversal::{LinearTerms, state_view},
};
use crate::{
    Act, Adjoint, Conjugate, ExpressionError, HilbertSpace, InnerProduct, LocalOperator,
    LocalState, Operator, RealPart, State,
    coefficient::{add, conjugate, multiply, zero},
};

type Factors<'a, I, S, O, C> = Vec<&'a StateExpr<I, S, O, C>>;
type Pair<'a, I, S, O, C> = (Factors<'a, I, S, O, C>, Factors<'a, I, S, O, C>);

// Per-call structural information, not a persistent support or numeric cache.
// None denotes an empty sum, including tensors annihilated by an empty sum.
// Some(0) denotes a scalar state on the unit space, which is different from zero.
struct FactorCounts<I, S, O, C>(BTreeMap<*const StateNode<I, S, O, C>, Option<usize>>);

impl<I, S, O, C> FactorCounts<I, S, O, C> {
    fn new(lhs: &StateExpr<I, S, O, C>, rhs: &StateExpr<I, S, O, C>) -> Self {
        let mut counts = Self(BTreeMap::new());
        let mut pending = vec![(lhs, false), (rhs, false)];
        while let Some((state, visited)) = pending.pop() {
            let key = ptr::from_ref(state.node());
            if counts.0.contains_key(&key) {
                continue;
            }
            if !visited {
                pending.push((state, true));
                match state.node() {
                    StateNode::Scale(s) => pending.push((&s.operand, false)),
                    StateNode::Sum(s) | StateNode::Tensor(s) => {
                        pending.extend(s.iter().map(|child| (child, false)));
                    }
                    StateNode::Local(_, _) | StateNode::Unit => {}
                    StateNode::Act(_) => unreachable!("evaluated state"),
                }
                continue;
            }
            let count = match state.node() {
                StateNode::Local(_, _) => Some(1),
                StateNode::Unit => Some(0),
                StateNode::Scale(s) => counts.get(&s.operand),
                // All nonempty alternatives have the same ordered factor space.
                StateNode::Sum(s) => s.iter().find_map(|child| counts.get(child)),
                StateNode::Tensor(s) => s
                    .iter()
                    .try_fold(0, |n, child| Some(n + counts.get(child)?)),
                StateNode::Act(_) => unreachable!("evaluated state"),
            };
            counts.0.insert(key, count);
        }
        counts
    }

    fn get(&self, state: &StateExpr<I, S, O, C>) -> Option<usize> {
        self.0[&ptr::from_ref(state.node())]
    }
}

enum Task<'a, I, S, O, C> {
    Contract(Pair<'a, I, S, O, C>),
    Scale(&'a C, bool),
    Sum(usize),
    Product(usize),
}

// Flatten tensor grouping only. Sums stay intact until a contraction needs them.
fn tensor_factors<'a, I, S, O, C>(
    mut pending: Factors<'a, I, S, O, C>,
    tasks: &mut Vec<Task<'a, I, S, O, C>>,
    conjugated: bool,
) -> Factors<'a, I, S, O, C> {
    pending.reverse();
    let mut factors = Vec::new();
    while let Some(state) = pending.pop() {
        match state.node() {
            StateNode::Tensor(s) => pending.extend(s.iter().rev()),
            StateNode::Unit => {}
            StateNode::Scale(s) => {
                tasks.push(Task::Scale(&s.coefficient, conjugated));
                pending.push(&s.operand);
            }
            StateNode::Local(_, _) | StateNode::Sum(_) => factors.push(state),
            StateNode::Act(_) => unreachable!("evaluated state"),
        }
    }
    factors
}

// Root-space compatibility already verifies every Index and its order. Within
// that order, equal prefix lengths identify the same contiguous subsystem.
// Find all common boundaries, even when the two tensor groupings differ.
fn matching_groups<'a, I, S, O, C>(
    mut lhs: Factors<'a, I, S, O, C>,
    mut rhs: Factors<'a, I, S, O, C>,
    counts: &FactorCounts<I, S, O, C>,
) -> Vec<Pair<'a, I, S, O, C>> {
    let mut groups = Vec::new();
    // Scalar sums do not cover any indexed factor. Separate them wherever
    // they occur, so they cannot obscure an otherwise shared tensor boundary.
    lhs.retain(|s| {
        if counts.get(s) == Some(0) {
            groups.push((vec![*s], Vec::new()));
            false
        } else {
            true
        }
    });
    rhs.retain(|s| {
        if counts.get(s) == Some(0) {
            groups.push((Vec::new(), vec![*s]));
            false
        } else {
            true
        }
    });
    let mut lhs = lhs.into_iter().peekable();
    let mut rhs = rhs.into_iter().peekable();
    while lhs.peek().is_some() || rhs.peek().is_some() {
        let left = lhs.next().expect("compatible factor spaces");
        let right = rhs.next().expect("compatible factor spaces");
        let mut a = vec![left];
        let mut b = vec![right];
        let mut na = counts.get(left).expect("nonzero factor");
        let mut nb = counts.get(right).expect("nonzero factor");
        while na != nb {
            if na < nb {
                let next = lhs.next().expect("compatible factor spaces");
                na += counts.get(next).expect("nonzero factor");
                a.push(next);
            } else {
                let next = rhs.next().expect("compatible factor spaces");
                nb += counts.get(next).expect("nonzero factor");
                b.push(next);
            }
        }
        groups.push((a, b));
    }
    groups
}

pub(crate) fn contract<I, S, O, C>(lhs: &StateExpr<I, S, O, C>, rhs: &StateExpr<I, S, O, C>) -> C
where
    C: Clone + Zero + One + Conjugate<Output = C>,
    S: InnerProduct<C>,
{
    let counts = FactorCounts::new(lhs, rhs);
    let mut tasks = vec![Task::Contract((vec![lhs], vec![rhs]))];
    let mut values = Vec::new();
    while let Some(task) = tasks.pop() {
        match task {
            Task::Contract((lhs, rhs)) => {
                if lhs.iter().chain(&rhs).any(|s| counts.get(s).is_none()) {
                    values.push(zero());
                    continue;
                }
                let lhs = tensor_factors(lhs, &mut tasks, true);
                let rhs = tensor_factors(rhs, &mut tasks, false);
                let mut groups = matching_groups(lhs, rhs, &counts);
                match groups.len() {
                    0 => values.push(C::one()),
                    1 => {
                        let (lhs, rhs) = groups.pop().expect("one subsystem");
                        if let ([a], [b]) = (lhs.as_slice(), rhs.as_slice()) {
                            if let (StateNode::Local(_, a), StateNode::Local(_, b)) =
                                (a.node(), b.node())
                            {
                                values.push(a.evaluated_inner_product(b));
                                continue;
                            }
                        }
                        // No further tensor boundary is shared. Distribute one
                        // sum here, leaving any separable surroundings outside.
                        let left_sum = lhs
                            .iter()
                            .position(|s| matches!(s.node(), StateNode::Sum(_)));
                        let (side, other, conjugated, position) = match left_sum {
                            Some(position) => (&lhs, &rhs, true, position),
                            None => {
                                let position = rhs
                                    .iter()
                                    .position(|s| matches!(s.node(), StateNode::Sum(_)))
                                    .expect("inseparable subsystem contains a sum");
                                (&rhs, &lhs, false, position)
                            }
                        };
                        let StateNode::Sum(children) = side[position].node() else {
                            unreachable!()
                        };
                        tasks.push(Task::Sum(children.len()));
                        for child in children.iter().rev() {
                            let mut branch = side.clone();
                            branch[position] = child;
                            let pair = if conjugated {
                                (branch, other.clone())
                            } else {
                                (other.clone(), branch)
                            };
                            tasks.push(Task::Contract(pair));
                        }
                    }
                    count => {
                        tasks.push(Task::Product(count));
                        tasks.extend(groups.into_iter().rev().map(Task::Contract));
                    }
                }
            }
            Task::Scale(coefficient, conjugated) => {
                let value = values.pop().expect("scaled contraction");
                let c = coefficient.clone();
                values.push(multiply(if conjugated { conjugate(c) } else { c }, value));
            }
            Task::Sum(count) | Task::Product(count) => {
                let parts = values.split_off(values.len() - count);
                let value = if matches!(task, Task::Sum(_)) {
                    parts.into_iter().fold(zero(), add)
                } else {
                    parts.into_iter().fold(C::one(), multiply)
                };
                values.push(value);
            }
        }
    }
    values.pop().expect("one contraction result")
}

impl<S, O, C> LocalState<S, O, C> {
    /// Computes the inner product, conjugate-linear on the left.
    ///
    /// Actions are evaluated internally; no basis or explicit prior evaluation
    /// is needed. Coefficients require arithmetic but no division or remainder.
    pub fn inner_product(&self, rhs: &Self) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
        O: Adjoint<C, Term = O>,
    {
        let lhs = self.evaluate();
        let rhs = rhs.evaluate();
        lhs.evaluated_inner_product(&rhs)
    }

    // Both operands must be action-free. Separating this reduction avoids
    // requiring action/adjoint capabilities again in composite contraction.
    pub(crate) fn evaluated_inner_product(&self, rhs: &Self) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        S: InnerProduct<C>,
    {
        let mut result = zero();
        for (a, left) in LinearTerms::new(self, state_view) {
            for (b, right) in LinearTerms::new(rhs, state_view) {
                let mut value = left.inner_product(right);
                if let Some(b) = b {
                    value = multiply(value, b);
                }
                if let Some(a) = a.clone() {
                    value = multiply(conjugate(a), value);
                }
                result = add(result, value);
            }
        }
        result
    }

    /// Returns the real part of the self-inner-product, without clamping.
    pub fn norm_squared(&self) -> C::Real
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C> + RealPart,
        S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
        O: Adjoint<C, Term = O>,
    {
        self.inner_product(self).real_part()
    }

    /// Computes the unnormalized expectation `self† operator self`.
    pub fn expectation(&self, operator: &LocalOperator<O, C>) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
        O: Adjoint<C, Term = O>,
    {
        self.inner_product(&self.act(operator))
    }
}

impl<I: Eq, S, O, C> State<I, S, O, C> {
    /// Computes the conjugate-linear/linear inner product on compatible spaces.
    ///
    /// Contracts matching tensor factors before distributing sums. Tensor
    /// grouping need not agree, and no global product-of-sums expansion is built.
    pub fn inner_product(&self, rhs: &Self) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
        O: Adjoint<C, Term = O>,
    {
        if !self.space.is_compatible_with(&rhs.space) {
            return Err(ExpressionError::IncompatibleSpaces);
        }
        let lhs = self.evaluate()?;
        let rhs = rhs.evaluate()?;
        Ok(super::inner_product::contract(
            &lhs.expression,
            &rhs.expression,
        ))
    }

    /// Returns the real part of the self-inner-product, without clamping.
    pub fn norm_squared(&self) -> Result<C::Real, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C> + RealPart,
        S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
        O: Adjoint<C, Term = O>,
    {
        Ok(self.inner_product(self)?.real_part())
    }

    /// Computes the unnormalized expectation without requiring prior evaluation.
    pub fn expectation(&self, operator: &Operator<I, O, C>) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
        O: Adjoint<C, Term = O>,
    {
        self.inner_product(&self.act(operator)?)
    }

    /// Computes an expectation using an explicit ambient space over any binding.
    pub fn expectation_in(
        &self,
        operator: &Operator<I, O, C>,
        ambient: &HilbertSpace<I>,
    ) -> Result<C, ExpressionError>
    where
        C: Clone + Zero + One + Conjugate<Output = C>,
        S: Act<C, Operator = O, Term = S> + InnerProduct<C>,
        O: Adjoint<C, Term = O>,
    {
        self.inner_product(&self.act_in(operator, ambient)?)
    }
}
