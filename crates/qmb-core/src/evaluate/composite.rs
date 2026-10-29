use alloc::{vec, vec::Vec};

use crate::expr::composite::{ExpressionError, OpExpr, OpNode, State, StateExpr, StateNode};
use crate::{
    Act, Adjoint, Conjugate, HilbertSpace, HomogeneousParity, LocalOperator, Parity,
    ParityOperator, TensorKind, coefficient::conjugate,
};

enum Task<I, S, O, C> {
    Evaluate(StateExpr<I, S, O, C>),
    Apply(
        StateExpr<I, S, O, C>,
        OpExpr<I, O, C>,
        bool,
        HilbertSpace<I>,
    ),
    ThenApply(OpExpr<I, O, C>, bool, HilbertSpace<I>),
    Local(
        StateExpr<I, S, O, C>,
        HilbertSpace<I>,
        LocalOperator<O, C>,
        Option<HilbertSpace<I>>,
    ),
    Placed(
        StateExpr<I, S, O, C>,
        HilbertSpace<I>,
        LocalOperator<O, C>,
        bool,
        HilbertSpace<I>,
    ),
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
        O: Adjoint<C, Term = O> + HomogeneousParity + ParityOperator,
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
                        StateNode::Sum(s) | StateNode::Tensor(_, s) => {
                            tasks.push(Task::Rebuild(state.clone(), s.len()));
                            for child in s.iter().rev() {
                                tasks.push(Task::Evaluate(child.clone()));
                            }
                        }
                        StateNode::Act(a) => {
                            tasks.push(Task::ThenApply(a.operator.clone(), false, a.space.clone()));
                            tasks.push(Task::Evaluate(a.state.clone()));
                        }
                    }
                }
                Task::ThenApply(op, adjoint, scope) => {
                    let state = values.pop().expect("evaluated action input");
                    tasks.push(Task::Apply(state, op, adjoint, scope));
                }
                Task::Apply(state, op, adjoint, scope) => match op.node() {
                    OpNode::Identity => values.push(state),
                    OpNode::Local(factor, local) => tasks.push(Task::Placed(
                        state,
                        factor.clone(),
                        local.clone(),
                        adjoint,
                        scope,
                    )),
                    OpNode::Scale(s) => {
                        tasks.push(Task::Scale(if adjoint {
                            conjugate(s.coefficient.clone())
                        } else {
                            s.coefficient.clone()
                        }));
                        tasks.push(Task::Apply(state, s.operand.clone(), adjoint, scope));
                    }
                    OpNode::Sum(s) => {
                        tasks.push(Task::Sum(s.len()));
                        for child in s.iter().rev() {
                            tasks.push(Task::Apply(
                                state.clone(),
                                child.clone(),
                                adjoint,
                                scope.clone(),
                            ));
                        }
                    }
                    OpNode::Product(p) => {
                        let (first, second) = if adjoint {
                            (&p.lhs, &p.rhs)
                        } else {
                            (&p.rhs, &p.lhs)
                        };
                        tasks.push(Task::ThenApply(second.clone(), adjoint, scope.clone()));
                        tasks.push(Task::Apply(state, first.clone(), adjoint, scope));
                    }
                    OpNode::Adjoint(a) => {
                        tasks.push(Task::Apply(state, a.operand.clone(), !adjoint, scope))
                    }
                },
                Task::Placed(state, factor, operator, adjoint, scope) => {
                    // Ordinary (including unit/single-factor) action can retain
                    // a whole local expression without distributing its sums.
                    if scope.effective_kind() != Some(TensorKind::Graded) {
                        let operator = if adjoint {
                            operator.adjoint()
                        } else {
                            operator
                        };
                        tasks.push(Task::Local(state, factor, operator, None));
                        continue;
                    }
                    use crate::expr::local::OperatorNode;
                    match operator.node() {
                        OperatorNode::Scale(s) => {
                            tasks.push(Task::Scale(if adjoint {
                                conjugate(s.coefficient.clone())
                            } else {
                                s.coefficient.clone()
                            }));
                            tasks.push(Task::Placed(
                                state,
                                factor,
                                s.operand.clone(),
                                adjoint,
                                scope,
                            ));
                        }
                        OperatorNode::Sum(s) => {
                            tasks.push(Task::Sum(s.len()));
                            for child in s.iter().rev() {
                                tasks.push(Task::Placed(
                                    state.clone(),
                                    factor.clone(),
                                    child.clone(),
                                    adjoint,
                                    scope.clone(),
                                ));
                            }
                        }
                        OperatorNode::Product(p) => {
                            let (first, second) = if adjoint {
                                (&p.lhs, &p.rhs)
                            } else {
                                (&p.rhs, &p.lhs)
                            };
                            // Keep both actions in the original selected scope.
                            let second = OpExpr::new(OpNode::Local(factor.clone(), second.clone()));
                            tasks.push(Task::ThenApply(second, adjoint, scope.clone()));
                            tasks.push(Task::Placed(state, factor, first.clone(), adjoint, scope));
                        }
                        OperatorNode::Adjoint(a) => tasks.push(Task::Placed(
                            state,
                            factor,
                            a.operand.clone(),
                            !adjoint,
                            scope,
                        )),
                        OperatorNode::Primitive(p) if adjoint => {
                            let local = LocalOperator::from_linear_combination(
                                p.adjoint().map(LocalOperator::primitive),
                            );
                            tasks.push(Task::Placed(state, factor, local, false, scope));
                        }
                        OperatorNode::Primitive(p) => {
                            let prefix = if scope.effective_kind() == Some(TensorKind::Graded)
                                && p.parity() == Parity::Odd
                            {
                                Some(scope)
                            } else {
                                None
                            };
                            tasks.push(Task::Local(state, factor, operator, prefix));
                        }
                    }
                }
                Task::Local(state, factor, operator, prefix) => match state.node() {
                    StateNode::Local(target, local) => {
                        if target.factor_index() == factor.factor_index() {
                            values.push(StateExpr::new(StateNode::Local(
                                target.clone(),
                                local.act(&operator).evaluate(),
                            )));
                        } else if prefix.as_ref().is_some_and(|space| {
                            space
                                .indices()
                                .take_while(|i| Some(*i) != factor.factor_index())
                                .any(|i| Some(i) == target.factor_index())
                        }) {
                            values.push(StateExpr::new(StateNode::Local(
                                target.clone(),
                                local
                                    .act(&LocalOperator::primitive(O::parity_operator()))
                                    .evaluate(),
                            )));
                        } else {
                            values.push(state);
                        }
                    }
                    StateNode::Unit => values.push(state),
                    StateNode::Scale(s) => {
                        tasks.push(Task::Rebuild(state.clone(), 1));
                        tasks.push(Task::Local(s.operand.clone(), factor, operator, prefix));
                    }
                    StateNode::Sum(s) | StateNode::Tensor(_, s) => {
                        tasks.push(Task::Rebuild(state.clone(), s.len()));
                        for child in s.iter().rev() {
                            tasks.push(Task::Local(
                                child.clone(),
                                factor.clone(),
                                operator.clone(),
                                prefix.clone(),
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
                        StateNode::Sum(s) | StateNode::Tensor(_, s) => {
                            if s.iter().zip(&parts).all(|(a, b)| a.same_node(b)) {
                                values.push(original);
                            } else {
                                let parts = parts.into_iter().collect();
                                values.push(StateExpr::new(
                                    if matches!(original.node(), StateNode::Tensor(_, _)) {
                                        StateNode::Tensor(
                                            match original.node() {
                                                StateNode::Tensor(k, _) => *k,
                                                _ => unreachable!(),
                                            },
                                            parts,
                                        )
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
