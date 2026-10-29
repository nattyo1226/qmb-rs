use alloc::{vec, vec::Vec};
use core::ops::Mul;

use crate::expr::local::{LocalOperator, LocalState, OperatorNode, StateNode};
use crate::{
    Act, Adjoint, Conjugate, Product,
    coefficient::{conjugate, multiply},
};

enum StateTask<S, O, C> {
    Evaluate(LocalState<S, O, C>),
    Apply(LocalState<S, O, C>, LocalOperator<O, C>, bool),
    ThenApply(LocalOperator<O, C>, bool),
    Scale(C),
    Sum(usize),
}

impl<S, O, C> LocalState<S, O, C> {
    /// Evaluates actions using local rules and linearity, without choosing a basis.
    ///
    /// No operator `Product` capability is required: `(AB) psi` applies `B`,
    /// then `A`. Ready subtrees are shared; duplicate and zero terms remain.
    pub fn evaluate(&self) -> Self
    where
        C: Clone + Conjugate<Output = C>,
        S: Act<C, Operator = O, Term = S>,
        O: Adjoint<C, Term = O>,
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
                        let op =
                            LocalOperator::from_terms(op.adjoint().map(LocalOperator::primitive));
                        tasks.push(StateTask::Apply(state, op, false));
                    }
                    OperatorNode::Primitive(op) => match state.node() {
                        StateNode::Primitive(primitive) => {
                            values.push(Self::from_terms(primitive.act(op).map(Self::primitive)))
                        }
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
}

enum OperatorTask<O, C> {
    Visit(LocalOperator<O, C>, bool),
    Multiply(LocalOperator<O, C>, LocalOperator<O, C>),
    Product,
    Scale(C),
    Sum(usize),
}

impl<O, C> LocalOperator<O, C> {
    pub(crate) fn evaluate(&self) -> Self
    where
        C: Clone + Mul<Output = C> + Conjugate<Output = C>,
        O: Product<C, Term = O> + Adjoint<C, Term = O>,
    {
        let mut tasks = vec![OperatorTask::Visit(self.clone(), false)];
        let mut values: Vec<Self> = Vec::new();
        while let Some(task) = tasks.pop() {
            match task {
                OperatorTask::Visit(op, adjoint) => match op.node() {
                    OperatorNode::Primitive(p) => {
                        if adjoint {
                            values.push(Self::from_terms(p.adjoint().map(Self::primitive)));
                        } else {
                            values.push(op);
                        }
                    }
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
                        values.push(Self::from_terms(a.product(b).map(Self::primitive)));
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
}
