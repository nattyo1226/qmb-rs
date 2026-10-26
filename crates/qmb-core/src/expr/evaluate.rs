use alloc::{vec, vec::Vec};
use core::ops::Mul;
use num_traits::Zero;

use super::{
    local::{LocalOperator, LocalState, OperatorNode, StateNode},
    parts::SumExpr,
};
use crate::{
    Act, Adjoint, Conjugate, Identity, InnerProduct, LocalStructure, Product, RealPart, Trace,
    coefficient::{add, conjugate, multiply, zero},
};

enum StateTask<L: LocalStructure, C> {
    Evaluate(LocalState<L, C>),
    Apply(LocalState<L, C>, LocalOperator<L, C>, bool),
    ThenApply(LocalOperator<L, C>, bool),
    Scale(C),
    Sum(usize),
}

impl<L: LocalStructure, C> LocalState<L, C> {
    /// Evaluates actions using local rules and linearity, without choosing a basis.
    ///
    /// No operator `Product` capability is required: `(AB) psi` applies `B`,
    /// then `A`. Ready subtrees are shared; duplicate and zero terms remain.
    pub fn evaluate(&self) -> Self
    where
        C: Clone + Conjugate<Output = C>,
        L::State: Act<C, Operator = L::Operator, Term = L::State>,
        L::Operator: Adjoint<C, Term = L::Operator>,
    {
        let mut tasks = vec![StateTask::Evaluate(self.clone())];
        let mut values: Vec<Self> = Vec::new();
        while let Some(task) = tasks.pop() {
            match task {
                StateTask::Evaluate(state) => {
                    if !state.has_action {
                        values.push(state);
                        continue;
                    }
                    match state.node() {
                        StateNode::Primitive(_) => unreachable!("primitive is ready"),
                        StateNode::Scale(s) => {
                            tasks.push(StateTask::Scale(s.coefficient.clone()));
                            tasks.push(StateTask::Evaluate(s.operand.clone()));
                        }
                        StateNode::Sum(s) => {
                            tasks.push(StateTask::Sum(s.len()));
                            for child in s.iter().rev() {
                                tasks.push(StateTask::Evaluate(child.clone()));
                            }
                        }
                        StateNode::Act(a) => {
                            tasks.push(StateTask::ThenApply(a.operator.clone(), false));
                            tasks.push(StateTask::Evaluate(a.state.clone()));
                        }
                    }
                }
                StateTask::ThenApply(operator, adjoint) => {
                    let state = values.pop().expect("evaluated action operand");
                    tasks.push(StateTask::Apply(state, operator, adjoint));
                }
                StateTask::Apply(state, operator, adjoint) => match operator.node() {
                    OperatorNode::Identity => values.push(state),
                    OperatorNode::Scale(s) => {
                        let c = s.coefficient.clone();
                        tasks.push(StateTask::Scale(if adjoint { conjugate(c) } else { c }));
                        tasks.push(StateTask::Apply(state, s.operand.clone(), adjoint));
                    }
                    OperatorNode::Sum(s) => {
                        tasks.push(StateTask::Sum(s.len()));
                        for op in s.iter().rev() {
                            tasks.push(StateTask::Apply(state.clone(), op.clone(), adjoint));
                        }
                    }
                    OperatorNode::Product(p) => {
                        let (first, second) = if adjoint {
                            (&p.lhs, &p.rhs)
                        } else {
                            (&p.rhs, &p.lhs)
                        };
                        tasks.push(StateTask::ThenApply(second.clone(), adjoint));
                        tasks.push(StateTask::Apply(state, first.clone(), adjoint));
                    }
                    OperatorNode::Adjoint(a) => {
                        tasks.push(StateTask::Apply(state, a.operand.clone(), !adjoint))
                    }
                    OperatorNode::Primitive(op) if adjoint => {
                        let op = LocalOperator::from_linear_combination(
                            op.adjoint().map(LocalOperator::primitive),
                        );
                        tasks.push(StateTask::Apply(state, op, false));
                    }
                    OperatorNode::Primitive(op) => match state.node() {
                        StateNode::Primitive(primitive) => values.push(
                            Self::from_linear_combination(primitive.act(op).map(Self::primitive)),
                        ),
                        StateNode::Scale(s) => {
                            tasks.push(StateTask::Scale(s.coefficient.clone()));
                            tasks.push(StateTask::Apply(
                                s.operand.clone(),
                                operator.clone(),
                                false,
                            ));
                        }
                        StateNode::Sum(s) => {
                            tasks.push(StateTask::Sum(s.len()));
                            for child in s.iter().rev() {
                                tasks.push(StateTask::Apply(
                                    child.clone(),
                                    operator.clone(),
                                    false,
                                ));
                            }
                        }
                        StateNode::Act(_) => unreachable!("action input has been evaluated"),
                    },
                },
                StateTask::Scale(coefficient) => {
                    let value = values.pop().expect("scale operand");
                    values.push(value.into_scaled(coefficient));
                }
                StateTask::Sum(count) => {
                    let children = values.split_off(values.len() - count);
                    values.push(Self::sum(children));
                }
            }
        }
        values.pop().expect("one state result")
    }

    /// Computes the inner product, conjugate-linear on the left.
    ///
    /// Actions are evaluated internally; no basis or explicit prior evaluation
    /// is needed. Coefficients require arithmetic but no division or remainder.
    pub fn inner_product(&self, rhs: &Self) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
        L::Operator: Adjoint<C, Term = L::Operator>,
    {
        let lhs = self.evaluate();
        let rhs = rhs.evaluate();
        lhs.evaluated_inner_product(&rhs)
    }

    // Both operands must be action-free. Separating this reduction avoids
    // requiring action/adjoint capabilities again in composite contraction.
    pub(super) fn evaluated_inner_product(&self, rhs: &Self) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        L::State: InnerProduct<C>,
    {
        let mut result = zero();
        for (a, left) in LinearTerms::new(self, state_view) {
            for (b, right) in LinearTerms::new(rhs, state_view) {
                let mut value = left
                    .expect("state primitive")
                    .inner_product(right.expect("state primitive"));
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
        L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
        L::Operator: Adjoint<C, Term = L::Operator>,
    {
        self.inner_product(self).real_part()
    }

    /// Computes the unnormalized expectation `self† operator self`.
    pub fn expectation(&self, operator: &LocalOperator<L, C>) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        L::State: Act<C, Operator = L::Operator, Term = L::State> + InnerProduct<C>,
        L::Operator: Adjoint<C, Term = L::Operator>,
    {
        self.inner_product(&self.act(operator))
    }
}

// Both state and operator scalar reductions use this same linear walker.
// None as coefficient means an implicit one, so generating it needs no One.
enum LinearView<'a, E, P, C> {
    Leaf(Option<&'a P>),
    Scale(&'a C, &'a E),
    Sum(&'a SumExpr<E>),
}
struct LinearTerms<'a, E, P, C> {
    pending: Vec<(&'a E, Option<C>)>,
    view: fn(&'a E) -> LinearView<'a, E, P, C>,
}
impl<'a, E, P, C> LinearTerms<'a, E, P, C> {
    fn new(root: &'a E, view: fn(&'a E) -> LinearView<'a, E, P, C>) -> Self {
        Self {
            pending: vec![(root, None)],
            view,
        }
    }
}
impl<'a, E, P: 'a, C> Iterator for LinearTerms<'a, E, P, C>
where
    C: Clone + Mul<Output = C>,
{
    type Item = (Option<C>, Option<&'a P>);
    fn next(&mut self) -> Option<Self::Item> {
        while let Some((expression, coefficient)) = self.pending.pop() {
            match (self.view)(expression) {
                LinearView::Leaf(leaf) => return Some((coefficient, leaf)),
                LinearView::Scale(c, operand) => {
                    let coefficient = Some(match coefficient {
                        None => c.clone(),
                        Some(a) => multiply(a, c.clone()),
                    });
                    self.pending.push((operand, coefficient));
                }
                LinearView::Sum(s) => {
                    for child in s.iter().rev() {
                        self.pending.push((child, coefficient.clone()));
                    }
                }
            }
        }
        None
    }
}
fn state_view<L: LocalStructure, C>(
    state: &LocalState<L, C>,
) -> LinearView<'_, LocalState<L, C>, L::State, C> {
    match state.node() {
        StateNode::Primitive(p) => LinearView::Leaf(Some(p)),
        StateNode::Scale(s) => LinearView::Scale(&s.coefficient, &s.operand),
        StateNode::Sum(s) => LinearView::Sum(s),
        StateNode::Act(_) => unreachable!("evaluated state"),
    }
}
fn operator_view<L: LocalStructure, C>(
    op: &LocalOperator<L, C>,
) -> LinearView<'_, LocalOperator<L, C>, L::Operator, C> {
    match op.node() {
        OperatorNode::Primitive(p) => LinearView::Leaf(Some(p)),
        OperatorNode::Identity => LinearView::Leaf(None),
        OperatorNode::Scale(s) => LinearView::Scale(&s.coefficient, &s.operand),
        OperatorNode::Sum(s) => LinearView::Sum(s),
        _ => unreachable!("evaluated operator"),
    }
}

enum OperatorTask<L: LocalStructure, C> {
    Visit(LocalOperator<L, C>, bool),
    Multiply(LocalOperator<L, C>, LocalOperator<L, C>),
    Product,
    Scale(C),
    Sum(usize),
}

impl<L: LocalStructure, C> LocalOperator<L, C> {
    fn evaluated(&self) -> Self
    where
        C: Clone + Mul<Output = C> + Conjugate<Output = C>,
        L::Operator: Product<C, Term = L::Operator> + Adjoint<C, Term = L::Operator>,
    {
        let mut tasks = vec![OperatorTask::Visit(self.clone(), false)];
        let mut values: Vec<Self> = Vec::new();
        while let Some(task) = tasks.pop() {
            match task {
                OperatorTask::Visit(op, adjoint) => match op.node() {
                    OperatorNode::Primitive(p) => {
                        if adjoint {
                            values.push(Self::from_linear_combination(
                                p.adjoint().map(Self::primitive),
                            ));
                        } else {
                            values.push(op);
                        }
                    }
                    OperatorNode::Identity => values.push(op),
                    OperatorNode::Scale(s) => {
                        let c = s.coefficient.clone();
                        tasks.push(OperatorTask::Scale(if adjoint { conjugate(c) } else { c }));
                        tasks.push(OperatorTask::Visit(s.operand.clone(), adjoint));
                    }
                    OperatorNode::Sum(s) => {
                        tasks.push(OperatorTask::Sum(s.len()));
                        for child in s.iter().rev() {
                            tasks.push(OperatorTask::Visit(child.clone(), adjoint));
                        }
                    }
                    OperatorNode::Adjoint(a) => {
                        tasks.push(OperatorTask::Visit(a.operand.clone(), !adjoint))
                    }
                    OperatorNode::Product(p) => {
                        let (lhs, rhs) = if adjoint {
                            (&p.rhs, &p.lhs)
                        } else {
                            (&p.lhs, &p.rhs)
                        };
                        tasks.push(OperatorTask::Product);
                        tasks.push(OperatorTask::Visit(rhs.clone(), adjoint));
                        tasks.push(OperatorTask::Visit(lhs.clone(), adjoint));
                    }
                },
                OperatorTask::Product => {
                    let rhs = values.pop().expect("right product operand");
                    let lhs = values.pop().expect("left product operand");
                    tasks.push(OperatorTask::Multiply(lhs, rhs));
                }
                OperatorTask::Multiply(lhs, rhs) => match (lhs.node(), rhs.node()) {
                    (OperatorNode::Identity, _) => values.push(rhs),
                    (_, OperatorNode::Identity) => values.push(lhs),
                    (OperatorNode::Scale(s), _) => {
                        tasks.push(OperatorTask::Scale(s.coefficient.clone()));
                        tasks.push(OperatorTask::Multiply(s.operand.clone(), rhs));
                    }
                    (OperatorNode::Sum(s), _) => {
                        tasks.push(OperatorTask::Sum(s.len()));
                        for child in s.iter().rev() {
                            tasks.push(OperatorTask::Multiply(child.clone(), rhs.clone()));
                        }
                    }
                    (_, OperatorNode::Scale(s)) => {
                        tasks.push(OperatorTask::Scale(s.coefficient.clone()));
                        tasks.push(OperatorTask::Multiply(lhs, s.operand.clone()));
                    }
                    (_, OperatorNode::Sum(s)) => {
                        tasks.push(OperatorTask::Sum(s.len()));
                        for child in s.iter().rev() {
                            tasks.push(OperatorTask::Multiply(lhs.clone(), child.clone()));
                        }
                    }
                    (OperatorNode::Primitive(a), OperatorNode::Primitive(b)) => {
                        values.push(Self::from_linear_combination(
                            a.product(b).map(Self::primitive),
                        ));
                    }
                    _ => unreachable!("linear operator operands"),
                },
                OperatorTask::Scale(coefficient) => {
                    let operand = values.pop().expect("scale operand");
                    // Combining adjacent evaluated scales avoids quadratic
                    // rewalking of coefficients in long single-term products.
                    let value = if let OperatorNode::Scale(s) = operand.node() {
                        s.operand
                            .scaled(multiply(coefficient, s.coefficient.clone()))
                    } else {
                        operand.into_scaled(coefficient)
                    };
                    values.push(value);
                }
                OperatorTask::Sum(count) => {
                    let children = values.split_off(values.len() - count);
                    values.push(Self::sum(children));
                }
            }
        }
        values.pop().expect("one operator result")
    }

    /// Computes the ordinary, unnormalized local trace without a chosen basis.
    pub fn trace(&self) -> C
    where
        C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
        L: Identity<C>,
        L::Operator: Product<C, Term = L::Operator> + Adjoint<C, Term = L::Operator> + Trace<C>,
    {
        let evaluated = self.evaluated();
        let mut result = zero();
        for (coefficient, primitive) in LinearTerms::new(&evaluated, operator_view) {
            let mut value = if let Some(p) = primitive {
                p.trace()
            } else {
                L::identity().into_iter().fold(zero(), |sum, term| {
                    add(sum, multiply(term.coefficient, term.operand.trace()))
                })
            };
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
        L: Identity<C>,
        L::Operator: Product<C, Term = L::Operator> + Adjoint<C, Term = L::Operator> + Trace<C>,
    {
        self.adjoint().product(rhs).trace()
    }
}
