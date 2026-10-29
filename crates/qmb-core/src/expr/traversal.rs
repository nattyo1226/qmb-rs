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

type SortedFactor<'a, I, O, C> = (usize, &'a HilbertSpace<I>, LocalOperator<O, C>, bool);
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
    grouped: bool,
    branches: Vec<OpBranch<'a, I, O, C>>,
}
impl<'a, I, O, C> OpTerms<'a, I, O, C> {
    pub(crate) fn new(root: &'a OpExpr<I, O, C>) -> Self {
        Self {
            grouped: true,
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
                        if let Some((_, existing)) = branch.factors.iter_mut().find(|(other, _)| {
                            self.grouped && other.factor_index() == factor.factor_index()
                        }) {
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

// Split only local expressions, never a global state or basis. Handles preserve
// borrowed/non-Clone primitives and scales do not require generating one.
fn parity_parts<O, C>(local: &LocalOperator<O, C>) -> Vec<(crate::Parity, LocalOperator<O, C>)>
where
    O: crate::HomogeneousParity + crate::Product<C, Term = O> + crate::Adjoint<C, Term = O>,
    C: Clone + Mul<Output = C> + Conjugate<Output = C>,
{
    let mut pending = vec![(local.evaluated(), None::<C>)];
    let mut even = Vec::new();
    let mut odd = Vec::new();
    while let Some((op, coefficient)) = pending.pop() {
        match op.node() {
            OperatorNode::Primitive(p) => {
                let parity = p.parity();
                let term = if let Some(c) = coefficient {
                    op.scaled(c)
                } else {
                    op
                };
                if parity == crate::Parity::Odd {
                    odd.push(term);
                } else {
                    even.push(term);
                }
            }
            OperatorNode::Scale(s) => {
                let c = match coefficient {
                    None => s.coefficient.clone(),
                    Some(c) => multiply(c, s.coefficient.clone()),
                };
                pending.push((s.operand.clone(), Some(c)));
            }
            OperatorNode::Sum(s) => pending.extend(
                s.iter()
                    .rev()
                    .map(|child| (child.clone(), coefficient.clone())),
            ),
            _ => unreachable!("evaluated local operator"),
        }
    }
    let mut parts = Vec::new();
    if !even.is_empty() {
        parts.push((crate::Parity::Even, LocalOperator::sum(even)));
    }
    if !odd.is_empty() {
        parts.push((crate::Parity::Odd, LocalOperator::sum(odd)));
    }
    parts
}

// Ordered homogeneous placements, before either Koszul sorting or lowering.
struct HomogeneousBranch<'a, I, O, C> {
    pending: OperatorFactors<'a, I, O, C>,
    factors: Vec<(&'a HilbertSpace<I>, LocalOperator<O, C>, crate::Parity)>,
    coefficient: Option<C>,
}
struct HomogeneousTerms<'a, I, O, C> {
    raw: OpTerms<'a, I, O, C>,
    branches: Vec<HomogeneousBranch<'a, I, O, C>>,
}
impl<'a, I, O, C> HomogeneousTerms<'a, I, O, C> {
    fn new(root: &'a OpExpr<I, O, C>) -> Self {
        let mut raw = OpTerms::new(root);
        raw.grouped = false;
        Self {
            raw,
            branches: Vec::new(),
        }
    }
}
impl<'a, I: Eq, O, C> Iterator for HomogeneousTerms<'a, I, O, C>
where
    O: crate::HomogeneousParity + crate::Product<C, Term = O> + crate::Adjoint<C, Term = O>,
    C: Clone + Mul<Output = C> + Conjugate<Output = C>,
{
    type Item = (
        Option<C>,
        Vec<(&'a HilbertSpace<I>, LocalOperator<O, C>, crate::Parity)>,
    );
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.branches.is_empty() {
                let (coefficient, mut pending) = self.raw.next()?;
                pending.reverse();
                self.branches.push(HomogeneousBranch {
                    pending,
                    factors: Vec::new(),
                    coefficient,
                });
            }
            let mut branch = self.branches.pop().expect("branch");
            if let Some((factor, local)) = branch.pending.pop() {
                for (parity, local) in parity_parts(&local).into_iter().rev() {
                    let mut factors = branch.factors.clone();
                    factors.push((factor, local, parity));
                    self.branches.push(HomogeneousBranch {
                        pending: branch.pending.clone(),
                        factors,
                        coefficient: branch.coefficient.clone(),
                    });
                }
            } else {
                return Some((branch.coefficient, branch.factors));
            }
        }
    }
}

/// Intrinsic graded factors. Never contains inserted parity strings.
pub(crate) struct IntrinsicOpTerms<'a, I, O, C> {
    ordinary: OpTerms<'a, I, O, C>,
    graded: Option<HomogeneousTerms<'a, I, O, C>>,
    ambient: &'a HilbertSpace<I>,
}
impl<'a, I, O, C> IntrinsicOpTerms<'a, I, O, C> {
    pub(crate) fn new(root: &'a OpExpr<I, O, C>, ambient: &'a HilbertSpace<I>) -> Self {
        Self {
            ordinary: OpTerms::new(root),
            graded: (ambient.effective_kind() == Some(crate::TensorKind::Graded))
                .then(|| HomogeneousTerms::new(root)),
            ambient,
        }
    }
}
impl<'a, I: Eq, O, C> Iterator for IntrinsicOpTerms<'a, I, O, C>
where
    O: crate::HomogeneousParity + crate::Product<C, Term = O> + crate::Adjoint<C, Term = O>,
    C: Clone
        + Mul<Output = C>
        + Conjugate<Output = C>
        + num_traits::One
        + core::ops::Neg<Output = C>,
{
    type Item = (Option<C>, OperatorFactors<'a, I, O, C>);
    fn next(&mut self) -> Option<Self::Item> {
        let Some(graded) = &mut self.graded else {
            return self.ordinary.next();
        };
        let (mut coefficient, factors) = graded.next()?;
        let mut sorted: Vec<SortedFactor<'a, I, O, C>> = Vec::new();
        let mut negative = false;
        for (factor, local, parity) in factors {
            let position = self
                .ambient
                .indices()
                .position(|i| Some(i) == factor.factor_index())
                .expect("validated factor");
            let odd = parity == crate::Parity::Odd;
            for (p, _, _, previous_odd) in &sorted {
                if *p > position && *previous_odd && odd {
                    negative = !negative;
                }
            }
            if let Some((_, _, previous, previous_odd)) =
                sorted.iter_mut().find(|(p, _, _, _)| *p == position)
            {
                *previous = previous.product(&local);
                *previous_odd ^= odd;
            } else {
                sorted.push((position, factor, local, odd));
            }
        }
        sorted.sort_by_key(|(position, _, _, _)| *position);
        if negative {
            coefficient = Some(-coefficient.unwrap_or_else(C::one));
        }
        Some((
            coefficient,
            sorted
                .into_iter()
                .map(|(_, f, local, _)| (f, local))
                .collect(),
        ))
    }
}

/// Ordinary representation after lowering placements, for full trace only.
/// These factors must never be interpreted as intrinsic graded coordinates.
pub(crate) struct OrdinaryOpTerms<'a, I, O, C> {
    ordinary: OpTerms<'a, I, O, C>,
    graded: Option<HomogeneousTerms<'a, I, O, C>>,
    ambient: &'a HilbertSpace<I>,
}
impl<'a, I, O, C> OrdinaryOpTerms<'a, I, O, C> {
    pub(crate) fn new(root: &'a OpExpr<I, O, C>, ambient: &'a HilbertSpace<I>) -> Self {
        Self {
            ordinary: OpTerms::new(root),
            graded: (ambient.effective_kind() == Some(crate::TensorKind::Graded))
                .then(|| HomogeneousTerms::new(root)),
            ambient,
        }
    }
}
impl<'a, I: Eq, O, C> Iterator for OrdinaryOpTerms<'a, I, O, C>
where
    O: crate::HomogeneousParity
        + crate::ParityOperator
        + crate::Product<C, Term = O>
        + crate::Adjoint<C, Term = O>,
    C: Clone + Mul<Output = C> + Conjugate<Output = C>,
{
    type Item = (Option<C>, OperatorFactors<'a, I, O, C>);
    fn next(&mut self) -> Option<Self::Item> {
        let Some(graded) = &mut self.graded else {
            return self.ordinary.next();
        };
        let (coefficient, placements) = graded.next()?;
        let mut factors: OperatorFactors<'a, I, O, C> = Vec::new();
        let mut append = |factor: &'a HilbertSpace<I>, local: LocalOperator<O, C>| {
            if let Some((_, previous)) = factors
                .iter_mut()
                .find(|(f, _)| f.factor_index() == factor.factor_index())
            {
                *previous = previous.product(&local);
            } else {
                factors.push((factor, local));
            }
        };
        for (factor, local, parity) in placements {
            if parity == crate::Parity::Odd {
                for prefix in self
                    .ambient
                    .leaves()
                    .take_while(|f| f.factor_index() != factor.factor_index())
                {
                    append(prefix, LocalOperator::primitive(O::parity_operator()));
                }
            }
            append(factor, local);
        }
        Some((coefficient, factors))
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
