use alloc::{vec, vec::Vec};
use core::ops::{Add, Mul, Neg, Sub};
use num_traits::{One, Zero};

use super::composite::{ExpressionError, OpExpr, OpNode, Operator, State, StateExpr, StateNode};
use crate::{
    Act, Adjoint, Complex, HilbertSpace, Identity, InnerProduct, LocalOperator, LocalStructure,
    Product, Trace,
    coefficient::{add, conjugate, multiply, zero},
};

enum Task<I, L: LocalStructure, T> {
    Evaluate(StateExpr<I, L, T>),
    Apply(StateExpr<I, L, T>, OpExpr<I, L, T>, bool),
    ThenApply(OpExpr<I, L, T>, bool),
    Local(StateExpr<I, L, T>, HilbertSpace<I, L>, LocalOperator<L, T>),
    Scale(Complex<T>),
    Sum(usize),
    Rebuild(StateExpr<I, L, T>, usize),
}

impl<I: Eq, L: LocalStructure, T> State<I, L, T> {
    /// Evaluates delayed actions using local rules, without a basis or matrix.
    ///
    /// An operator product applies its right operand first; primitive `Product`
    /// is not required. Ready subtrees, spectator states and spaces are shared.
    pub fn evaluate(&self) -> Result<Self, ExpressionError>
    where
        T: Clone + Neg<Output = T>,
        L::State: Act<T, Operator = L::Operator, Term = L::State>,
        L::Operator: Adjoint<T, Term = L::Operator>,
    {
        let mut tasks = vec![Task::Evaluate(self.expression.clone())];
        let mut values: Vec<StateExpr<I, L, T>> = Vec::new();
        while let Some(task) = tasks.pop() {
            match task {
                Task::Evaluate(state) => {
                    if !state.pending_action {
                        values.push(state);
                        continue;
                    }
                    match state.node() {
                        StateNode::Local(factor, local) => values.push(StateExpr::new(
                            StateNode::Local(factor.clone(), local.evaluate()),
                        )),
                        StateNode::Unit => unreachable!("unit is ready"),
                        StateNode::Scale(s) => {
                            tasks.push(Task::Rebuild(state.clone(), 1));
                            tasks.push(Task::Evaluate(s.operand.clone()));
                        }
                        StateNode::Sum(s) | StateNode::Tensor(s) => {
                            tasks.push(Task::Rebuild(state.clone(), s.len()));
                            for child in s.iter().rev() {
                                tasks.push(Task::Evaluate(child.clone()));
                            }
                        }
                        StateNode::Act(a) => {
                            tasks.push(Task::ThenApply(a.operator.clone(), false));
                            tasks.push(Task::Evaluate(a.state.clone()));
                        }
                    }
                }
                Task::ThenApply(op, adjoint) => {
                    let state = values.pop().expect("evaluated action input");
                    tasks.push(Task::Apply(state, op, adjoint));
                }
                Task::Apply(state, op, adjoint) => match op.node() {
                    OpNode::Identity => values.push(state),
                    OpNode::Local(factor, local) => tasks.push(Task::Local(
                        state,
                        factor.clone(),
                        if adjoint {
                            local.adjoint()
                        } else {
                            local.clone()
                        },
                    )),
                    OpNode::Scale(s) => {
                        tasks.push(Task::Scale(if adjoint {
                            conjugate(s.coefficient.clone())
                        } else {
                            s.coefficient.clone()
                        }));
                        tasks.push(Task::Apply(state, s.operand.clone(), adjoint));
                    }
                    OpNode::Sum(s) => {
                        tasks.push(Task::Sum(s.len()));
                        for child in s.iter().rev() {
                            tasks.push(Task::Apply(state.clone(), child.clone(), adjoint));
                        }
                    }
                    OpNode::Product(p) => {
                        let (first, second) = if adjoint {
                            (&p.lhs, &p.rhs)
                        } else {
                            (&p.rhs, &p.lhs)
                        };
                        tasks.push(Task::ThenApply(second.clone(), adjoint));
                        tasks.push(Task::Apply(state, first.clone(), adjoint));
                    }
                    OpNode::Adjoint(a) => {
                        tasks.push(Task::Apply(state, a.operand.clone(), !adjoint))
                    }
                },
                Task::Local(state, factor, operator) => match state.node() {
                    StateNode::Local(target, local) => {
                        if target.factor_index() == factor.factor_index() {
                            values.push(StateExpr::new(StateNode::Local(
                                target.clone(),
                                local.act(&operator).evaluate(),
                            )));
                        } else {
                            values.push(state);
                        }
                    }
                    StateNode::Unit => values.push(state),
                    StateNode::Scale(s) => {
                        tasks.push(Task::Rebuild(state.clone(), 1));
                        tasks.push(Task::Local(s.operand.clone(), factor, operator));
                    }
                    StateNode::Sum(s) | StateNode::Tensor(s) => {
                        tasks.push(Task::Rebuild(state.clone(), s.len()));
                        for child in s.iter().rev() {
                            tasks.push(Task::Local(
                                child.clone(),
                                factor.clone(),
                                operator.clone(),
                            ));
                        }
                    }
                    StateNode::Act(_) => unreachable!("local action input is ready"),
                },
                Task::Scale(coefficient) => {
                    let state = values.pop().expect("scale operand");
                    values.push(state.scaled(coefficient));
                }
                Task::Sum(count) => {
                    let parts = values.split_off(values.len() - count);
                    values.push(StateExpr::sum(parts.into_iter().collect()));
                }
                Task::Rebuild(original, count) => {
                    let parts = values.split_off(values.len() - count);
                    match original.node() {
                        StateNode::Scale(s) => {
                            let operand = parts.into_iter().next().expect("scale operand");
                            if operand.same_node(&s.operand) {
                                values.push(original);
                            } else {
                                values.push(operand.scaled(s.coefficient.clone()));
                            }
                        }
                        StateNode::Sum(s) | StateNode::Tensor(s) => {
                            if s.iter().zip(&parts).all(|(a, b)| a.same_node(b)) {
                                values.push(original);
                            } else {
                                let parts = parts.into_iter().collect();
                                values.push(StateExpr::new(
                                    if matches!(original.node(), StateNode::Tensor(_)) {
                                        StateNode::Tensor(parts)
                                    } else {
                                        StateNode::Sum(parts)
                                    },
                                ));
                            }
                        }
                        _ => unreachable!("linear state container"),
                    }
                }
            }
        }
        Ok(Self {
            space: self.space.clone(),
            expression: values.pop().expect("one state result"),
        })
    }

    /// Computes the conjugate-linear/linear inner product on compatible spaces.
    ///
    /// Contracts matching tensor factors before distributing sums. Tensor
    /// grouping need not agree, and no global product-of-sums expansion is built.
    pub fn inner_product(&self, rhs: &Self) -> Result<Complex<T>, ExpressionError>
    where
        T: Clone + Zero + One + Sub<Output = T> + Neg<Output = T>,
        L::State: Act<T, Operator = L::Operator, Term = L::State> + InnerProduct<T>,
        L::Operator: Adjoint<T, Term = L::Operator>,
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

    /// Computes the unnormalized expectation without requiring prior evaluation.
    pub fn expectation(&self, operator: &Operator<I, L, T>) -> Result<Complex<T>, ExpressionError>
    where
        T: Clone + Zero + One + Sub<Output = T> + Neg<Output = T>,
        L::State: Act<T, Operator = L::Operator, Term = L::State> + InnerProduct<T>,
        L::Operator: Adjoint<T, Term = L::Operator>,
    {
        self.inner_product(&self.act(operator)?)
    }

    /// Computes an expectation using an explicit ambient space over any binding.
    pub fn expectation_in(
        &self,
        operator: &Operator<I, L, T>,
        ambient: &HilbertSpace<I, L>,
    ) -> Result<Complex<T>, ExpressionError>
    where
        T: Clone + Zero + One + Sub<Output = T> + Neg<Output = T>,
        L::State: Act<T, Operator = L::Operator, Term = L::State> + InnerProduct<T>,
        L::Operator: Adjoint<T, Term = L::Operator>,
    {
        self.inner_product(&self.act_in(operator, ambient)?)
    }
}

fn one<T: Zero + One>() -> Complex<T> {
    Complex::new(T::one(), T::zero())
}

type OperatorFactors<'a, I, L, T> = Vec<(&'a HilbertSpace<I, L>, LocalOperator<L, T>)>;
struct OpBranch<'a, I, L: LocalStructure, T> {
    pending: Vec<(&'a OpExpr<I, L, T>, bool)>,
    coefficient: Option<Complex<T>>,
    factors: OperatorFactors<'a, I, L, T>,
}
impl<I, L: LocalStructure, T: Clone> Clone for OpBranch<'_, I, L, T> {
    fn clone(&self) -> Self {
        Self {
            pending: self.pending.clone(),
            coefficient: self.coefficient.clone(),
            factors: self.factors.clone(),
        }
    }
}
struct OpTerms<'a, I, L: LocalStructure, T> {
    branches: Vec<OpBranch<'a, I, L, T>>,
}
impl<'a, I, L: LocalStructure, T> OpTerms<'a, I, L, T> {
    fn new(root: &'a OpExpr<I, L, T>) -> Self {
        Self {
            branches: vec![OpBranch {
                pending: vec![(root, false)],
                coefficient: None,
                factors: Vec::new(),
            }],
        }
    }
}
impl<'a, I: Eq, L: LocalStructure, T> Iterator for OpTerms<'a, I, L, T>
where
    T: Clone + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Neg<Output = T>,
{
    type Item = (Option<Complex<T>>, OperatorFactors<'a, I, L, T>);
    fn next(&mut self) -> Option<Self::Item> {
        'branch: while let Some(mut branch) = self.branches.pop() {
            while let Some((op, adjoint)) = branch.pending.pop() {
                match op.node() {
                    OpNode::Identity => {}
                    OpNode::Local(factor, local) => {
                        let local = if adjoint {
                            local.adjoint()
                        } else {
                            local.clone()
                        };
                        if let Some((_, existing)) = branch
                            .factors
                            .iter_mut()
                            .find(|(other, _)| other.factor_index() == factor.factor_index())
                        {
                            *existing = existing.product(&local);
                        } else {
                            branch.factors.push((factor, local));
                        }
                    }
                    OpNode::Scale(s) => {
                        let c = if adjoint {
                            conjugate(s.coefficient.clone())
                        } else {
                            s.coefficient.clone()
                        };
                        branch.coefficient = Some(match branch.coefficient {
                            None => c,
                            Some(a) => multiply(a, c),
                        });
                        branch.pending.push((&s.operand, adjoint));
                    }
                    OpNode::Sum(s) => {
                        for child in s.iter().rev() {
                            let mut next = branch.clone();
                            next.pending.push((child, adjoint));
                            self.branches.push(next);
                        }
                        continue 'branch;
                    }
                    OpNode::Product(p) => {
                        let (lhs, rhs) = if adjoint {
                            (&p.rhs, &p.lhs)
                        } else {
                            (&p.lhs, &p.rhs)
                        };
                        branch.pending.push((rhs, adjoint));
                        branch.pending.push((lhs, adjoint));
                    }
                    OpNode::Adjoint(a) => branch.pending.push((&a.operand, !adjoint)),
                }
            }
            return Some((branch.coefficient, branch.factors));
        }
        None
    }
}

impl<I: Eq, L: LocalStructure, T> Operator<I, L, T> {
    /// Computes the ordinary full trace using the existing binding.
    pub fn trace(&self) -> Result<Complex<T>, ExpressionError>
    where
        T: Clone + Zero + One + Sub<Output = T> + Neg<Output = T>,
        L: Identity<T>,
        L::Operator: Product<T, Term = L::Operator> + Adjoint<T, Term = L::Operator> + Trace<T>,
    {
        let ambient = self.select_space(None, None)?;
        self.trace_on(ambient)
    }

    /// Computes full trace on an explicit space, overriding an old binding.
    pub fn trace_in(&self, ambient: &HilbertSpace<I, L>) -> Result<Complex<T>, ExpressionError>
    where
        T: Clone + Zero + One + Sub<Output = T> + Neg<Output = T>,
        L: Identity<T>,
        L::Operator: Product<T, Term = L::Operator> + Adjoint<T, Term = L::Operator> + Trace<T>,
    {
        let ambient = self.select_space(Some(ambient), None)?;
        self.trace_on(ambient)
    }

    fn trace_on(&self, ambient: &HilbertSpace<I, L>) -> Result<Complex<T>, ExpressionError>
    where
        T: Clone + Zero + One + Sub<Output = T> + Neg<Output = T>,
        L: Identity<T>,
        L::Operator: Product<T, Term = L::Operator> + Adjoint<T, Term = L::Operator> + Trace<T>,
    {
        let mut result = zero();
        for (coefficient, factors) in OpTerms::new(&self.expression) {
            let mut value = one();
            for index in ambient.indices() {
                let local_trace = match factors
                    .iter()
                    .find(|(factor, _)| factor.factor_index() == Some(index))
                {
                    Some((_, local)) => local.trace(),
                    None => LocalOperator::<L, T>::identity().trace(),
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
    pub fn hilbert_schmidt_inner_product(&self, rhs: &Self) -> Result<Complex<T>, ExpressionError>
    where
        T: Clone + Zero + One + Sub<Output = T> + Neg<Output = T>,
        L: Identity<T>,
        L::Operator: Product<T, Term = L::Operator> + Adjoint<T, Term = L::Operator> + Trace<T>,
    {
        self.adjoint().product(rhs)?.trace()
    }

    /// Computes Hilbert-Schmidt inner product on an explicit ambient space.
    pub fn hilbert_schmidt_inner_product_in(
        &self,
        rhs: &Self,
        ambient: &HilbertSpace<I, L>,
    ) -> Result<Complex<T>, ExpressionError>
    where
        T: Clone + Zero + One + Sub<Output = T> + Neg<Output = T>,
        L: Identity<T>,
        L::Operator: Product<T, Term = L::Operator> + Adjoint<T, Term = L::Operator> + Trace<T>,
    {
        self.bind(ambient)?
            .adjoint()
            .product(&rhs.bind(ambient)?)?
            .trace()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::num::NonZeroUsize;

    #[derive(Clone)]
    struct Ring(i64);
    macro_rules! arithmetic {
        ($trait:ident, $method:ident, $op:tt) => {
            impl core::ops::$trait for Ring {
                type Output = Self;
                fn $method(self, rhs: Self) -> Self { Self(self.0 $op rhs.0) }
            }
        };
    }
    arithmetic!(Add, add, +);
    arithmetic!(Sub, sub, -);
    arithmetic!(Mul, mul, *);
    impl Neg for Ring {
        type Output = Self;
        fn neg(self) -> Self {
            Self(-self.0)
        }
    }
    struct Metadata;
    impl LocalStructure for Metadata {
        type State = ();
        type Operator = ();
        const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
    }

    #[test]
    fn operator_term_traversal_needs_no_zero_or_one() {
        // Term traversal only combines stored coefficients and lazy local
        // operators. Generating numeric identities belongs to trace reduction.
        let lhs = Operator::<u8, Metadata, Ring>::scalar(Complex::new(Ring(3), Ring(4)));
        let rhs = Operator::scalar(Complex::new(Ring(2), Ring(-1)));
        let operator = lhs.product(&rhs).unwrap().adjoint();
        let mut terms = OpTerms::new(&operator.expression);
        let (coefficient, factors) = terms.next().unwrap();
        let coefficient = coefficient.unwrap();
        assert_eq!((coefficient.re.0, coefficient.im.0), (10, -5));
        assert!(factors.is_empty());
        assert!(terms.next().is_none());
    }
}
