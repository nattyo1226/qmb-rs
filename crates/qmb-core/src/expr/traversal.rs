//! Internal linear-term and operator-term traversal shared by evaluation and bases.

use super::{
    composite::{OpExpr, OpNode},
    local::{LocalState, OperatorNode, StateNode},
    parts::SumExpr,
};
use crate::{
    Conjugate, HilbertSpace, LocalOperator,
    coefficient::{conjugate, multiply},
};
use alloc::{vec, vec::Vec};
use core::ops::Mul;

// Both state and operator scalar reductions use this same linear walker.
// None as coefficient means an implicit one, so generating it needs no One.
pub(crate) enum LinearView<'a, E, P, C> {
    Leaf(&'a P),
    Scale(&'a C, &'a E),
    Sum(&'a SumExpr<E>),
}
pub(crate) struct LinearTerms<'a, E, P, C> {
    pending: Vec<(&'a E, Option<C>)>,
    view: fn(&'a E) -> LinearView<'a, E, P, C>,
}
impl<'a, E, P, C> LinearTerms<'a, E, P, C> {
    pub(crate) fn new(root: &'a E, view: fn(&'a E) -> LinearView<'a, E, P, C>) -> Self {
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
    type Item = (Option<C>, &'a P);
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
pub(crate) fn state_view<S, O, C>(
    state: &LocalState<S, O, C>,
) -> LinearView<'_, LocalState<S, O, C>, S, C> {
    match state.node() {
        StateNode::Primitive(p) => LinearView::Leaf(p),
        StateNode::Scale(s) => LinearView::Scale(&s.coefficient, &s.operand),
        StateNode::Sum(s) => LinearView::Sum(s),
        StateNode::Act(_) => unreachable!("evaluated state"),
    }
}
pub(crate) fn operator_view<O, C>(
    op: &LocalOperator<O, C>,
) -> LinearView<'_, LocalOperator<O, C>, O, C> {
    match op.node() {
        OperatorNode::Primitive(p) => LinearView::Leaf(p),
        OperatorNode::Scale(s) => LinearView::Scale(&s.coefficient, &s.operand),
        OperatorNode::Sum(s) => LinearView::Sum(s),
        _ => unreachable!("evaluated operator"),
    }
}

type OperatorFactors<'a, I, O, C> = Vec<(&'a HilbertSpace<I>, LocalOperator<O, C>)>;
struct OpBranch<'a, I, O, C> {
    pending: Vec<(&'a OpExpr<I, O, C>, bool)>,
    coefficient: Option<C>,
    factors: OperatorFactors<'a, I, O, C>,
}
impl<I, O, C: Clone> Clone for OpBranch<'_, I, O, C> {
    fn clone(&self) -> Self {
        Self {
            pending: self.pending.clone(),
            coefficient: self.coefficient.clone(),
            factors: self.factors.clone(),
        }
    }
}
pub(crate) struct OpTerms<'a, I, O, C> {
    branches: Vec<OpBranch<'a, I, O, C>>,
}
impl<'a, I, O, C> OpTerms<'a, I, O, C> {
    pub(crate) fn new(root: &'a OpExpr<I, O, C>) -> Self {
        Self {
            branches: vec![OpBranch {
                pending: vec![(root, false)],
                coefficient: None,
                factors: Vec::new(),
            }],
        }
    }
}
impl<'a, I: Eq, O, C> Iterator for OpTerms<'a, I, O, C>
where
    C: Clone + Mul<Output = C> + Conjugate<Output = C>,
{
    type Item = (Option<C>, OperatorFactors<'a, I, O, C>);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Operator, test_coefficient::MinimalComplex};
    use core::ops::Neg;

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

    #[test]
    fn operator_term_traversal_needs_no_zero_or_one() {
        // Term traversal only combines stored coefficients and lazy local
        // operators. Generating numeric identities belongs to trace reduction.
        let lhs =
            Operator::<u8, (), MinimalComplex<Ring>>::scalar(MinimalComplex::new(Ring(3), Ring(4)));
        let rhs = Operator::scalar(MinimalComplex::new(Ring(2), Ring(-1)));
        let operator = lhs.product(&rhs).unwrap().adjoint();
        let mut terms = OpTerms::new(&operator.expression);
        let (coefficient, factors) = terms.next().unwrap();
        let coefficient = coefficient.unwrap();
        assert_eq!((coefficient.re.0, coefficient.im.0), (10, -5));
        assert!(factors.is_empty());
        assert!(terms.next().is_none());
    }
}
