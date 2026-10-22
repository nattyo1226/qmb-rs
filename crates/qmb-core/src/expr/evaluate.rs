use alloc::{vec, vec::Vec};
use core::ops::{Add, Mul, Neg, Sub};
use num_traits::Zero;

use super::{
    local::{LocalOperator, LocalState, OperatorNode, StateNode},
    parts::SumExpr,
};
use crate::{Act, Adjoint, Complex, Identity, InnerProduct, LocalStructure, Product, Trace};

fn conjugate<T: Neg<Output = T>>(value: Complex<T>) -> Complex<T> {
    Complex::new(value.re, -value.im)
}
// Keep the coefficient bounds at the arithmetic operation. Complex<T>'s Num
// implementation would additionally require division/remainder not used here.
fn multiply<T>(a: Complex<T>, b: Complex<T>) -> Complex<T>
where
    T: Clone + Add<Output = T> + Sub<Output = T> + Mul<Output = T>,
{
    Complex::new(
        a.re.clone() * b.re.clone() - a.im.clone() * b.im.clone(),
        a.re * b.im + a.im * b.re,
    )
}
fn add<T: Add<Output = T>>(a: Complex<T>, b: Complex<T>) -> Complex<T> {
    Complex::new(a.re + b.re, a.im + b.im)
}
fn zero<T: Zero>() -> Complex<T> {
    Complex::new(T::zero(), T::zero())
}

enum StateTask<L: LocalStructure, T> {
    Evaluate(LocalState<L, T>),
    Apply(LocalState<L, T>, LocalOperator<L, T>, bool),
    ThenApply(LocalOperator<L, T>, bool),
    Scale(Complex<T>),
    Sum(usize),
}

impl<L: LocalStructure, T> LocalState<L, T> {
    /// Evaluates actions using local rules and linearity, without choosing a basis.
    ///
    /// No operator `Product` capability is required: `(AB) psi` applies `B`,
    /// then `A`. Ready subtrees are shared; duplicate and zero terms remain.
    pub fn evaluate(&self) -> Self
    where
        T: Clone + Neg<Output = T>,
        L::State: Act<T, Operator = L::Operator, Term = L::State>,
        L::Operator: Adjoint<T, Term = L::Operator>,
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
    pub fn inner_product(&self, rhs: &Self) -> Complex<T>
    where
        T: Clone + Zero + Sub<Output = T> + Mul<Output = T> + Neg<Output = T>,
        L::State: Act<T, Operator = L::Operator, Term = L::State> + InnerProduct<T>,
        L::Operator: Adjoint<T, Term = L::Operator>,
    {
        let lhs = self.evaluate();
        let rhs = rhs.evaluate();
        let mut result = zero();
        for (a, left) in LinearTerms::new(&lhs, state_view) {
            for (b, right) in LinearTerms::new(&rhs, state_view) {
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

    /// Computes the unnormalized expectation `self† operator self`.
    pub fn expectation(&self, operator: &LocalOperator<L, T>) -> Complex<T>
    where
        T: Clone + Zero + Sub<Output = T> + Mul<Output = T> + Neg<Output = T>,
        L::State: Act<T, Operator = L::Operator, Term = L::State> + InnerProduct<T>,
        L::Operator: Adjoint<T, Term = L::Operator>,
    {
        self.inner_product(&self.act(operator))
    }
}

// Both state and operator scalar reductions use this same linear walker.
// None as coefficient means an implicit one, so generating it needs no One.
enum LinearView<'a, E, P, T> {
    Leaf(Option<&'a P>),
    Scale(&'a Complex<T>, &'a E),
    Sum(&'a SumExpr<E>),
}
struct LinearTerms<'a, E, P, T> {
    pending: Vec<(&'a E, Option<Complex<T>>)>,
    view: fn(&'a E) -> LinearView<'a, E, P, T>,
}
impl<'a, E, P, T> LinearTerms<'a, E, P, T> {
    fn new(root: &'a E, view: fn(&'a E) -> LinearView<'a, E, P, T>) -> Self {
        Self {
            pending: vec![(root, None)],
            view,
        }
    }
}
impl<'a, E, P: 'a, T> Iterator for LinearTerms<'a, E, P, T>
where
    T: Clone + Add<Output = T> + Sub<Output = T> + Mul<Output = T>,
{
    type Item = (Option<Complex<T>>, Option<&'a P>);
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
fn state_view<L: LocalStructure, T>(
    state: &LocalState<L, T>,
) -> LinearView<'_, LocalState<L, T>, L::State, T> {
    match state.node() {
        StateNode::Primitive(p) => LinearView::Leaf(Some(p)),
        StateNode::Scale(s) => LinearView::Scale(&s.coefficient, &s.operand),
        StateNode::Sum(s) => LinearView::Sum(s),
        StateNode::Act(_) => unreachable!("evaluated state"),
    }
}
fn operator_view<L: LocalStructure, T>(
    op: &LocalOperator<L, T>,
) -> LinearView<'_, LocalOperator<L, T>, L::Operator, T> {
    match op.node() {
        OperatorNode::Primitive(p) => LinearView::Leaf(Some(p)),
        OperatorNode::Identity => LinearView::Leaf(None),
        OperatorNode::Scale(s) => LinearView::Scale(&s.coefficient, &s.operand),
        OperatorNode::Sum(s) => LinearView::Sum(s),
        _ => unreachable!("evaluated operator"),
    }
}

enum OperatorTask<L: LocalStructure, T> {
    Visit(LocalOperator<L, T>, bool),
    Multiply(LocalOperator<L, T>, LocalOperator<L, T>),
    Product,
    Scale(Complex<T>),
    Sum(usize),
}

impl<L: LocalStructure, T> LocalOperator<L, T> {
    fn evaluated(&self) -> Self
    where
        T: Clone + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Neg<Output = T>,
        L::Operator: Product<T, Term = L::Operator> + Adjoint<T, Term = L::Operator>,
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
    pub fn trace(&self) -> Complex<T>
    where
        T: Clone + Zero + Sub<Output = T> + Mul<Output = T> + Neg<Output = T>,
        L: Identity<T>,
        L::Operator: Product<T, Term = L::Operator> + Adjoint<T, Term = L::Operator> + Trace<T>,
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
    pub fn hilbert_schmidt_inner_product(&self, rhs: &Self) -> Complex<T>
    where
        T: Clone + Zero + Sub<Output = T> + Mul<Output = T> + Neg<Output = T>,
        L: Identity<T>,
        L::Operator: Product<T, Term = L::Operator> + Adjoint<T, Term = L::Operator> + Trace<T>,
    {
        self.adjoint().product(rhs).trace()
    }
}
