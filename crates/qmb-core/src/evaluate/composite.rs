use alloc::{vec, vec::Vec};

use crate::expr::composite::{ExpressionError, OpExpr, OpNode, State, StateExpr, StateNode};
use crate::{Act, Adjoint, Conjugate, HilbertSpace, LocalOperator, coefficient::conjugate};

enum Task<I, S, O, C> {
    Evaluate(StateExpr<I, S, O, C>),
    Apply(StateExpr<I, S, O, C>, OpExpr<I, O, C>, bool),
    ThenApply(OpExpr<I, O, C>, bool),
    Local(StateExpr<I, S, O, C>, HilbertSpace<I>, LocalOperator<O, C>),
    Scale(C),
    Sum(usize),
    Rebuild(StateExpr<I, S, O, C>, usize),
}

impl<I: Eq, S, O, C> State<I, S, O, C> {
    /// Evaluates delayed actions using local rules, without a basis or matrix.
    ///
    /// An operator product applies its right operand first; primitive `Product`
    /// is not required. Ready subtrees, spectator states and spaces are shared.
    pub fn evaluate(&self) -> Result<Self, ExpressionError>
    where
        C: Clone + Conjugate<Output = C>,
        S: Act<C, Operator = O, Term = S>,
        O: Adjoint<C, Term = O>,
    {
        let mut tasks = vec![Task::Evaluate(self.expression.clone())];
        let mut values: Vec<StateExpr<I, S, O, C>> = Vec::new();
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
}
