use super::composite::{OpNode, StateNode};
use crate::{
    Act, Adjoint, ExpressionError as E, HilbertSpace, InnerProduct, LinearCombination as Lc,
    LocalOperator, LocalState, Operator, State,
};
use alloc::vec::Vec;
use core::cell::Cell;
use num_complex::Complex;

#[derive(PartialEq, Eq)]
struct Index(u8);
#[derive(Default)]
struct Counts {
    acts: Cell<usize>,
    inner_products: Cell<usize>,
    drops: Cell<usize>,
}
struct Primitive<'a>(u8, &'a Counts);
impl Drop for Primitive<'_> {
    fn drop(&mut self) {
        self.1.drops.set(self.1.drops.get() + 1);
    }
}
enum Op {
    Flip,
    Parity,
}
#[allow(non_upper_case_globals)]
const Op: Op = Op::Flip;
impl qmb_core::HomogeneousParity for Op {
    fn parity(&self) -> qmb_core::Parity {
        qmb_core::Parity::Even
    }
}
impl qmb_core::ParityOperator for Op {
    fn parity_operator() -> Self {
        Self::Parity
    }
}

// No Product or Identity capability. Metadata, indices and primitives are not Clone.
impl<'a> Act<Complex<f64>> for Primitive<'a> {
    type Operator = Op;
    type Term = Self;
    fn act(&self, op: &Op) -> Lc<Self, Complex<f64>> {
        self.1.acts.set(self.1.acts.get() + 1);
        Lc::primitive(Primitive(
            if matches!(op, Op::Parity) {
                self.0
            } else {
                1 - self.0
            },
            self.1,
        ))
    }
}
impl Adjoint<Complex<f64>> for Op {
    type Term = Self;
    fn adjoint(&self) -> Lc<Self, Complex<f64>> {
        Lc::primitive(match self {
            Self::Flip => Self::Flip,
            Self::Parity => Self::Parity,
        })
    }
}
impl InnerProduct<Complex<f64>> for Primitive<'_> {
    fn inner_product(&self, rhs: &Self) -> Complex<f64> {
        self.1.inner_products.set(self.1.inner_products.get() + 1);
        Complex::new(if self.0 == rhs.0 { 1.0 } else { 0.0 }, 0.0)
    }
}
type H = HilbertSpace<Index>;
type S<'a, T = f64> = State<Index, Primitive<'a>, Op, Complex<T>>;
type O<T = f64> = Operator<Index, Op, Complex<T>>;

#[test]
fn evaluated_contraction_only_requires_primitive_inner_products() {
    struct InnerOnly;
    struct OpaqueIndex;

    impl InnerProduct<Complex<f64>> for InnerOnly {
        fn inner_product(&self, _: &Self) -> Complex<f64> {
            Complex::new(1.0, 0.0)
        }
    }
    // Neither Act nor Adjoint is implemented. Contraction of already evaluated
    // syntax needs no action capabilities, cloning of primitives, or Index Eq.
    let factor = HilbertSpace::<OpaqueIndex>::factor(OpaqueIndex);
    let local = LocalState::<InnerOnly, (), Complex<f64>>::primitive(InnerOnly);
    let state = local.at(&factor).unwrap();
    assert_eq!(
        local.evaluated_inner_product(&local),
        Complex::new(1.0, 0.0)
    );
    assert_eq!(
        crate::evaluate::inner_product::contract(&state.expression, &state.expression),
        Complex::new(1.0, 0.0)
    );
}

#[test]
fn inner_product_contracts_independent_sums_without_global_expansion() {
    let counts = Counts::default();
    for n in [8, 32] {
        let factors: Vec<_> = (0..n).map(|i| H::factor(Index(255 - i))).collect();
        let space = H::tensor_all(factors.iter().cloned()).unwrap();
        let before = S::tensor_in(
            &space,
            factors.iter().map(|factor| {
                (LocalState::primitive(Primitive(0, &counts))
                    + LocalState::primitive(Primitive(1, &counts)))
                .at(factor)
                .unwrap()
            }),
        )
        .unwrap();
        let after = S::tensor_in(
            &space,
            factors.iter().map(|factor| {
                let up = LocalState::primitive(Primitive(0, &counts))
                    .at(factor)
                    .unwrap();
                let down = LocalState::primitive(Primitive(1, &counts))
                    .at(factor)
                    .unwrap();
                (up.clone() + &down).unwrap()
            }),
        )
        .unwrap();
        let expected = Complex::new(2_f64.powi(i32::from(n)), 0.0);
        for (lhs, rhs) in [
            (&before, &before),
            (&before, &after),
            (&after, &before),
            (&after, &after),
        ] {
            counts.inner_products.set(0);
            assert_eq!(lhs.inner_product(rhs).unwrap(), expected);
            assert_eq!(counts.inner_products.get(), 4 * usize::from(n));
        }
        let op = LocalOperator::primitive(Op).at(&factors[0]).unwrap();
        counts.inner_products.set(0);
        assert_eq!(after.expectation(&op).unwrap(), expected);
        assert_eq!(counts.inner_products.get(), 4 * usize::from(n));
    }
}

#[test]
fn inner_product_aligns_overlapping_subsystems_and_scalar_factors() {
    let counts = Counts::default();
    let factors = [5, 1, 9, 3].map(|i| H::factor(Index(i)));
    let ket = |site: usize, value| {
        LocalState::primitive(Primitive(value, &counts))
            .at(&factors[site])
            .unwrap()
    };
    let ab = (S::tensor_in(
        &H::tensor_all([factors[0].clone(), factors[1].clone()]).unwrap(),
        [ket(0, 0), ket(1, 0)],
    )
    .unwrap()
    .clone()
        + &S::tensor_in(
            &H::tensor_all([factors[0].clone(), factors[1].clone()]).unwrap(),
            [ket(0, 1), ket(1, 1)],
        )
        .unwrap()
        .scaled(Complex::new(0.0, 1.0)))
        .unwrap();
    let cd = (S::tensor_in(
        &H::tensor_all([factors[2].clone(), factors[3].clone()]).unwrap(),
        [ket(2, 0), ket(3, 0)],
    )
    .unwrap()
    .clone()
        + &S::tensor_in(
            &H::tensor_all([factors[2].clone(), factors[3].clone()]).unwrap(),
            [ket(2, 1), ket(3, 1)],
        )
        .unwrap()
        .scaled(Complex::new(1.0, 1.0)))
        .unwrap();
    let lhs = S::tensor_in(
        &H::tensor_all(factors.iter().cloned()).unwrap(),
        [ab.clone(), cd.clone()],
    )
    .unwrap();
    let a = (ket(0, 0).clone() + &ket(0, 1).scaled(Complex::new(2.0, 1.0))).unwrap();
    let bc = (S::tensor_in(
        &H::tensor_all([factors[1].clone(), factors[2].clone()]).unwrap(),
        [ket(1, 0), ket(2, 0)],
    )
    .unwrap()
    .clone()
        + &S::tensor_in(
            &H::tensor_all([factors[1].clone(), factors[2].clone()]).unwrap(),
            [ket(1, 1), ket(2, 1)],
        )
        .unwrap()
        .scaled(Complex::new(2.0, -1.0)))
        .unwrap();
    let d = (ket(3, 0).clone() + &ket(3, 1).scaled(Complex::new(0.0, 1.0))).unwrap();
    let rhs = S::tensor_in(
        &H::tensor_all(factors.iter().cloned()).unwrap(),
        [
            a,
            S::tensor_in(
                &H::tensor_all(factors[1..].iter().cloned()).unwrap(),
                [bc, d],
            )
            .unwrap(),
        ],
    )
    .unwrap();
    // Independent amplitudes: the only shared nonzero coordinates are 0000
    // and 1111, with amplitudes (1, -1+i) and (1, 5i), respectively.
    assert_eq!(lhs.inner_product(&rhs).unwrap(), Complex::new(6.0, -5.0));
    assert_eq!(rhs.inner_product(&lhs).unwrap(), Complex::new(6.0, 5.0));
    assert_eq!(lhs.inner_product(&lhs).unwrap(), Complex::new(6.0, 0.0));
    assert_eq!(rhs.inner_product(&rhs).unwrap(), Complex::new(72.0, 0.0));

    let unit = S::unit();
    let scalar = (unit.scaled(Complex::new(2.0, 1.0)).clone()
        + &unit.scaled(Complex::new(-1.0, 2.0)))
        .unwrap();
    let with_scalar = S::tensor_in(
        lhs.space(),
        [
            unit.clone(),
            ab.clone(),
            scalar.clone(),
            cd.clone(),
            unit.clone(),
        ],
    )
    .unwrap();
    assert_eq!(
        with_scalar.inner_product(&rhs).unwrap(),
        Complex::new(-9.0, -23.0)
    );
    assert_eq!(
        rhs.inner_product(&with_scalar).unwrap(),
        Complex::new(-9.0, 23.0)
    );
    assert_eq!(
        scalar.inner_product(&scalar).unwrap(),
        Complex::new(10.0, 0.0)
    );
    let zero = S::zero(lhs.space());
    let with_zero = (zero.clone() + &with_scalar).unwrap();
    assert_eq!(
        with_zero.inner_product(&rhs).unwrap(),
        Complex::new(-9.0, -23.0)
    );
    let zero_tensor = S::tensor_in(lhs.space(), [S::zero(unit.space()), lhs.clone()]).unwrap();
    assert_eq!(
        zero_tensor.inner_product(&rhs).unwrap(),
        Complex::new(0.0, 0.0)
    );
    assert_eq!(
        lhs.inner_product(&zero_tensor).unwrap(),
        Complex::new(0.0, 0.0)
    );
}

#[test]
fn construction_uses_opaque_scalars_and_shared_nonclone_borrowed_leaves() {
    struct Scalar;
    let counts = Counts::default();
    let h = H::factor(Index(0));
    let local = LocalState::<Primitive<'_>, Op, Complex<Scalar>>::primitive(Primitive(0, &counts));
    let state = local.at(&h).unwrap();
    let op = LocalOperator::primitive(Op).at(&h).unwrap();
    let a = (state.clone() + &state)
        .unwrap()
        .scaled(Complex::new(Scalar, Scalar));
    let action = a.act(&op).unwrap();
    let clone = action.clone();
    assert!(clone.space().shares_node_with(&h));
    let _ = S::<Scalar>::tensor_in(&H::unit(), []).unwrap();
    let _ = S::<Scalar>::zero(&h);
    let _ = (op.clone() * &op).unwrap().adjoint();
    assert_eq!(counts.acts.get(), 0);
    drop(local);
    drop(state);
    drop(a);
    drop(action);
    assert_eq!(counts.drops.get(), 0);
    drop(clone);
    assert_eq!(counts.drops.get(), 1);
}

#[test]
fn action_uses_no_product_shares_spectators_and_preserves_original_inputs() {
    let counts = Counts::default();
    let a = H::factor(Index(1));
    let b = H::factor(Index(2));
    let left = LocalState::primitive(Primitive(0, &counts)).at(&a).unwrap();
    let right = LocalState::primitive(Primitive(1, &counts)).at(&b).unwrap();
    let state = S::tensor_in(&a.tensor(&b).unwrap(), [left.clone(), right]).unwrap();
    let op = LocalOperator::primitive(Op).at(&a).unwrap();
    let action = state.act(&op).unwrap();
    assert_eq!(counts.acts.get(), 0);
    let result = action.evaluate().unwrap();
    assert_eq!(counts.acts.get(), 1);
    let StateNode::Tensor(before) = state.expression.node() else {
        panic!("tensor");
    };
    let StateNode::Tensor(after) = result.expression.node() else {
        panic!("tensor");
    };
    assert!(before.as_slice()[1].shares_node_with(&after.as_slice()[1]));
    assert!(!state.expression.pending_action && action.expression.pending_action);
    assert!(!result.expression.pending_action);
    assert!(
        result
            .evaluate()
            .unwrap()
            .expression
            .shares_node_with(&result.expression)
    );
    let identity = O::identity();
    assert!(
        state
            .act(&identity)
            .unwrap()
            .evaluate()
            .unwrap()
            .expression
            .shares_node_with(&state.expression)
    );
    let x2 = (op.clone() * &op).unwrap();
    let p = state.act(&x2).unwrap().evaluate().unwrap();
    let q = state
        .act(&op)
        .unwrap()
        .act(&op)
        .unwrap()
        .evaluate()
        .unwrap();
    assert_eq!(p.inner_product(&q).unwrap(), Complex::new(1.0, 0.0));
    assert_eq!(state.inner_product(&p).unwrap(), Complex::new(1.0, 0.0));
}

#[test]
fn binding_merges_and_explicit_contexts_follow_the_selected_space() {
    let counts = Counts::default();
    let a = H::factor(Index(1));
    let b = H::factor(Index(2));
    let ab = a.tensor(&b).unwrap();
    let op = LocalOperator::<Op, Complex<f64>>::primitive(Op)
        .at(&a)
        .unwrap();
    let small = op.bind(&a).unwrap();
    let large = op.bind(&ab).unwrap();
    assert!(matches!(small.clone() + &large, Err(E::IncompatibleSpaces)));
    assert!(matches!(small.clone() * &large, Err(E::IncompatibleSpaces)));
    // The third, unbound operand suppresses comparisons of the first two.
    let mixed = O::sum([small.clone(), large.clone(), O::identity()]).unwrap();
    assert!(mixed.space().is_none());
    assert!(
        O::product([small.clone(), large.clone(), O::identity()])
            .unwrap()
            .space()
            .is_none()
    );
    let rebound = mixed.bind(&ab).unwrap();
    assert!(rebound.domain().unwrap().is_compatible_with(&ab));
    assert!(rebound.codomain().unwrap().is_compatible_with(&ab));
    assert!(small.space().unwrap().is_compatible_with(&a));
    let left = LocalState::primitive(Primitive(0, &counts)).at(&a).unwrap();
    let right = LocalState::primitive(Primitive(0, &counts)).at(&b).unwrap();
    let state = S::tensor_in(&a.tensor(&b).unwrap(), [left.clone(), right]).unwrap();
    assert!(matches!(state.act(&small), Err(E::IncompatibleSpaces)));
    assert!(state.act_in(&small, &ab).is_ok());
    assert!(state.act(&op).is_ok());
    assert!(matches!(left.act(&large), Err(E::IncompatibleSpaces)));
    assert!(left.act_in(&large, &a).is_ok());
    assert!(matches!(op.bind(&b), Err(E::MissingIndex)));
    assert!(matches!(
        op.scaled(Complex::new(0.0, 0.0)).bind(&b),
        Err(E::MissingIndex)
    ));
    assert!(matches!(
        LocalOperator::<Op, Complex<f64>>::zero()
            .at(&a)
            .unwrap()
            .bind(&b),
        Err(E::MissingIndex)
    ));
    assert!(matches!(
        LocalOperator::<Op, Complex<f64>>::primitive(Op)
            .at(&a)
            .unwrap()
            .bind(&b),
        Err(E::MissingIndex)
    ));
    assert!(matches!(
        LocalState::<Primitive<'_>, Op, Complex<f64>>::zero().at(&ab),
        Err(E::NotAFactor)
    ));
    assert!(matches!(
        S::<f64>::zero(&a).clone() + &S::zero(&b),
        Err(E::IncompatibleSpaces)
    ));
    assert_eq!(counts.acts.get(), 0);
    let OpNode::Sum(children) = rebound.expression.node() else {
        panic!("sum");
    };
    assert_eq!(children.len(), 3);
}

#[test]
fn evaluate_does_not_require_coefficient_arithmetic() {
    #[derive(Clone)]
    struct Scalar;
    impl crate::Conjugate for Scalar {
        type Output = Self;
        fn conjugate(&self) -> Self {
            Self
        }
    }
    impl<'a> Act<Scalar> for Primitive<'a> {
        type Operator = Op;
        type Term = Self;
        fn act(&self, _: &Op) -> Lc<Self, Scalar> {
            Lc::term(Scalar, Primitive(self.0, self.1))
        }
    }
    impl Adjoint<Scalar> for Op {
        type Term = Self;
        fn adjoint(&self) -> Lc<Self, Scalar> {
            Lc::term(
                Scalar,
                match self {
                    Self::Flip => Self::Flip,
                    Self::Parity => Self::Parity,
                },
            )
        }
    }
    let counts = Counts::default();
    let factor = H::factor(Index(1));
    let state = LocalState::<Primitive<'_>, Op, Scalar>::primitive(Primitive(0, &counts))
        .at(&factor)
        .unwrap();
    let op = LocalOperator::primitive(Op).at(&factor).unwrap();
    assert!(state.act(&op).unwrap().evaluate().is_ok());
}

#[test]
fn scalar_syntax_needs_no_index_equality_or_numeric_capabilities() {
    struct OpaqueIndex;
    struct Scalar(u8);
    let coefficient = || Complex::new(Scalar(2), Scalar(3));
    let counts = Counts::default();
    let factor = HilbertSpace::<OpaqueIndex>::factor(OpaqueIndex);
    let state = LocalState::<Primitive<'_>, Op, Complex<Scalar>>::primitive(Primitive(0, &counts))
        .at(&factor)
        .unwrap();
    let operator = LocalOperator::<Op, Complex<Scalar>>::primitive(Op)
        .at(&factor)
        .unwrap();
    for scaled in [
        state.clone() * coefficient(),
        &state * coefficient(),
        state.scaled(coefficient()),
        state.scaled(coefficient()),
    ] {
        let StateNode::Scale(s) = scaled.expression.node() else {
            panic!("scale")
        };
        assert_eq!((s.coefficient.re.0, s.coefficient.im.0), (2, 3));
        assert!(s.operand.shares_node_with(&state.expression));
        assert!(scaled.space().shares_node_with(&factor));
    }
    for scaled in [
        operator.clone() * coefficient(),
        &operator * coefficient(),
        operator.scaled(coefficient()),
        operator.scaled(coefficient()),
    ] {
        let OpNode::Scale(s) = scaled.expression.node() else {
            panic!("scale")
        };
        assert_eq!((s.coefficient.re.0, s.coefficient.im.0), (2, 3));
        assert!(core::ptr::eq(s.operand.node(), operator.expression.node()));
        assert!(scaled.space().is_none());
    }

    // Borrowing the coefficient adds only Clone, with no arithmetic or Index Eq.
    #[derive(Clone)]
    struct CloneOnly(u8);
    let coefficient = Complex::new(CloneOnly(4), CloneOnly(5));
    let state =
        LocalState::<Primitive<'_>, Op, Complex<CloneOnly>>::primitive(Primitive(1, &counts))
            .at(&factor)
            .unwrap();
    let operator = LocalOperator::<Op, Complex<CloneOnly>>::primitive(Op)
        .at(&factor)
        .unwrap();
    for scaled in [
        state.clone() * &coefficient,
        &state * &coefficient,
        state.scaled(coefficient.clone()),
        state.scaled(coefficient.clone()),
    ] {
        let StateNode::Scale(s) = scaled.expression.node() else {
            panic!("scale")
        };
        assert_eq!((s.coefficient.re.0, s.coefficient.im.0), (4, 5));
        assert!(s.operand.shares_node_with(&state.expression));
    }
    for scaled in [
        operator.clone() * &coefficient,
        &operator * &coefficient,
        operator.scaled(coefficient.clone()),
        operator.scaled(coefficient.clone()),
    ] {
        let OpNode::Scale(s) = scaled.expression.node() else {
            panic!("scale")
        };
        assert_eq!((s.coefficient.re.0, s.coefficient.im.0), (4, 5));
        assert!(core::ptr::eq(s.operand.node(), operator.expression.node()));
    }
    assert_eq!(counts.acts.get(), 0);
}

#[test]
#[allow(clippy::op_ref)] // Exercise the public borrowed-scalar implementations.
fn all_binary_ownership_forms_preserve_results_and_shared_operands() {
    let counts = Counts::default();
    let factor = H::factor(Index(0));
    let state = LocalState::primitive(Primitive(0, &counts))
        .at(&factor)
        .unwrap();
    let op = LocalOperator::primitive(Op)
        .at(&factor)
        .unwrap()
        .bind(&factor)
        .unwrap();
    let original = state.expression.clone();
    macro_rules! forms {
        ($a:expr, $symbol:tt, $b:expr) => {
            [$a.clone() $symbol $b.clone(), $a.clone() $symbol &$b,
             &$a $symbol $b.clone(), &$a $symbol &$b]
        };
    }
    for sum in forms!(state, +, state) {
        let sum = sum.unwrap();
        let StateNode::Sum(parts) = sum.expression.node() else {
            panic!("sum")
        };
        assert!(parts.iter().all(|p| p.shares_node_with(&original)));
        assert_eq!(state.inner_product(&sum).unwrap(), Complex::new(2.0, 0.0));
    }
    for difference in forms!(state, -, state) {
        assert_eq!(
            difference.unwrap().inner_product(&state).unwrap(),
            Complex::new(0.0, 0.0)
        );
    }
    for sum in forms!(op, +, op) {
        let sum = sum.unwrap();
        assert!(sum.space().unwrap().shares_node_with(&factor));
        let OpNode::Sum(parts) = sum.expression.node() else {
            panic!("sum")
        };
        assert!(
            parts
                .iter()
                .all(|p| core::ptr::eq(p.node(), op.expression.node()))
        );
        let result = state.act(&sum).unwrap();
        assert_eq!(
            result.inner_product(&result).unwrap(),
            Complex::new(4.0, 0.0)
        );
    }
    for difference in forms!(op, -, op) {
        let result = state.act(&difference.unwrap()).unwrap();
        assert_eq!(
            result.inner_product(&result).unwrap(),
            Complex::new(0.0, 0.0)
        );
    }
    for product in forms!(op, *, op) {
        assert_eq!(
            state.expectation(&product.unwrap()).unwrap(),
            Complex::new(1.0, 0.0)
        );
    }
    for action in forms!(op, *, state) {
        let action = action.unwrap();
        assert!(action.space().shares_node_with(&factor));
        assert_eq!(
            action.inner_product(&action).unwrap(),
            Complex::new(1.0, 0.0)
        );
        assert_eq!(
            state.inner_product(&action).unwrap(),
            Complex::new(0.0, 0.0)
        );
    }
    let coefficient = Complex::new(2.0, -1.0);
    for scaled in [
        state.clone() * coefficient,
        &state * coefficient,
        coefficient * state.clone(),
        coefficient * &state,
        state.clone() * &coefficient,
        &state * &coefficient,
        &coefficient * state.clone(),
        &coefficient * &state,
    ] {
        assert_eq!(state.inner_product(&scaled).unwrap(), coefficient);
        assert_eq!(
            scaled.inner_product(&state).unwrap(),
            Complex::new(2.0, 1.0)
        );
    }
    for scaled in [
        op.clone() * coefficient,
        &op * coefficient,
        coefficient * op.clone(),
        coefficient * &op,
        op.clone() * &coefficient,
        &op * &coefficient,
        &coefficient * op.clone(),
        &coefficient * &op,
    ] {
        let result = state.act(&scaled).unwrap();
        assert_eq!(
            result.inner_product(&result).unwrap(),
            Complex::new(5.0, 0.0)
        );
    }
    for negated in [-state.clone(), -&state] {
        assert_eq!(
            state.inner_product(&negated).unwrap(),
            Complex::new(-1.0, 0.0)
        );
    }
    for negated in [-op.clone(), -&op] {
        let expected = state.act(&op).unwrap();
        assert_eq!(
            expected
                .inner_product(&state.act(&negated).unwrap())
                .unwrap(),
            Complex::new(-1.0, 0.0)
        );
    }

    let other_factor = H::factor(Index(1));
    let other_state = LocalState::primitive(Primitive(0, &counts))
        .at(&other_factor)
        .unwrap();
    let other_op = LocalOperator::primitive(Op)
        .at(&other_factor)
        .unwrap()
        .bind(&other_factor)
        .unwrap();
    for result in forms!(state, +, other_state)
        .into_iter()
        .chain(forms!(state, -, other_state))
        .chain(forms!(op, *, other_state))
    {
        assert!(matches!(result, Err(E::IncompatibleSpaces)));
    }
    for result in forms!(op, +, other_op)
        .into_iter()
        .chain(forms!(op, -, other_op))
        .chain(forms!(op, *, other_op))
    {
        assert!(matches!(result, Err(E::IncompatibleSpaces)));
    }
    assert!(state.expression.shares_node_with(&original));
}

#[test]
fn unit_zero_singleton_and_placement_preserve_spaces_and_lazy_local_expressions() {
    let counts = Counts::default();
    let factor = H::factor(Index(1));
    let unit = S::<f64>::unit();
    let leaf = LocalState::primitive(Primitive(0, &counts))
        .at(&factor)
        .unwrap();
    let nested = S::tensor_in(
        &H::tensor_all([unit.space().clone(), factor.clone()]).unwrap(),
        [unit.clone(), S::tensor_in(&factor, [leaf.clone()]).unwrap()],
    )
    .unwrap();
    assert!(nested.space().is_compatible_with(&factor));
    assert_eq!(nested.inner_product(&leaf).unwrap(), Complex::new(1.0, 0.0));
    assert!(leaf.domain().is_compatible_with(unit.space()));
    assert!(leaf.codomain().shares_node_with(&factor));
    assert_eq!(
        S::sum(&factor, []).unwrap().inner_product(&leaf).unwrap(),
        Complex::new(0.0, 0.0)
    );
    assert!(
        S::sum(&factor, [leaf.clone()])
            .unwrap()
            .expression
            .shares_node_with(&leaf.expression)
    );
    assert!(matches!(
        S::tensor_in(&factor, [leaf.clone(), leaf.clone()]),
        Err(E::IncompatibleSpaces)
    ));
    assert!(matches!(
        O::<f64>::place(&LocalOperator::primitive(Op), unit.space()),
        Err(E::NotAFactor)
    ));
    assert!(matches!(
        unit.act(&LocalOperator::primitive(Op).at(&factor).unwrap()),
        Err(E::MissingIndex)
    ));
    // Placing an entire local action keeps the local expression lazy until evaluate.
    let local = LocalState::primitive(Primitive(0, &counts)).act(&LocalOperator::primitive(Op));
    let delayed = local.at(&factor).unwrap();
    assert_eq!(counts.acts.get(), 0);
    let result = delayed.evaluate().unwrap();
    assert_eq!(counts.acts.get(), 1);
    assert_eq!(leaf.inner_product(&result).unwrap(), Complex::new(0.0, 0.0));
}

#[test]
fn shared_composite_types_follow_stored_primitive_data() {
    fn send_sync<T: Send + Sync>() {}
    send_sync::<HilbertSpace<u8>>();
    send_sync::<State<u8, (), (), Complex<f64>>>();
    send_sync::<Operator<u8, (), Complex<f64>>>();
}

#[test]
fn composite_scalar_reductions_require_neither_num_float_copy_nor_division() {
    use crate::test_coefficient::MinimalComplex;
    use crate::{Identity, Product, Trace};
    use core::ops::{Add, Mul, Neg, Sub};
    use num_traits::{One, Zero};

    #[derive(Clone, Debug, PartialEq)]
    struct Ring(i64);
    impl Add for Ring {
        type Output = Self;
        fn add(self, rhs: Self) -> Self {
            Self(self.0 + rhs.0)
        }
    }
    impl Sub for Ring {
        type Output = Self;
        fn sub(self, rhs: Self) -> Self {
            Self(self.0 - rhs.0)
        }
    }
    impl Mul for Ring {
        type Output = Self;
        fn mul(self, rhs: Self) -> Self {
            Self(self.0 * rhs.0)
        }
    }
    impl Neg for Ring {
        type Output = Self;
        fn neg(self) -> Self {
            Self(-self.0)
        }
    }
    impl Zero for Ring {
        fn zero() -> Self {
            Self(0)
        }
        fn is_zero(&self) -> bool {
            self.0 == 0
        }
    }
    impl One for Ring {
        fn one() -> Self {
            Self(1)
        }
    }
    struct Primitive;
    impl qmb_core::HomogeneousParity for Primitive {
        fn parity(&self) -> qmb_core::Parity {
            qmb_core::Parity::Even
        }
    }
    impl qmb_core::ParityOperator for Primitive {
        fn parity_operator() -> Self {
            Self
        }
    }

    impl Identity for Primitive {
        fn identity() -> Self {
            Self
        }
    }
    impl Act<MinimalComplex<Ring>> for Primitive {
        type Operator = Self;
        type Term = Self;
        fn act(&self, _: &Self) -> Lc<Self, MinimalComplex<Ring>> {
            Lc::primitive(Primitive)
        }
    }
    impl Product<MinimalComplex<Ring>> for Primitive {
        type Term = Self;
        fn product(&self, _: &Self) -> Lc<Self, MinimalComplex<Ring>> {
            Lc::primitive(Primitive)
        }
    }
    impl Adjoint<MinimalComplex<Ring>> for Primitive {
        type Term = Self;
        fn adjoint(&self) -> Lc<Self, MinimalComplex<Ring>> {
            Lc::primitive(Primitive)
        }
    }
    impl InnerProduct<MinimalComplex<Ring>> for Primitive {
        fn inner_product(&self, _: &Self) -> MinimalComplex<Ring> {
            MinimalComplex::new(Ring(1), Ring(0))
        }
    }
    impl Trace<MinimalComplex<Ring>> for Primitive {
        fn trace(&self) -> MinimalComplex<Ring> {
            MinimalComplex::new(Ring(1), Ring(0))
        }
    }
    let a = HilbertSpace::<_>::factor(0u8);
    let b = HilbertSpace::factor(1);
    let state = qmb_core::State::tensor_in(
        &a.tensor(&b).unwrap(),
        [
            LocalState::<Primitive, Primitive, MinimalComplex<Ring>>::primitive(Primitive)
                .at(&a)
                .unwrap(),
            LocalState::primitive(Primitive).at(&b).unwrap(),
        ],
    )
    .unwrap()
    .scaled(MinimalComplex::new(Ring(3), Ring(4)));
    let operator = LocalOperator::primitive(Primitive)
        .at(&a)
        .unwrap()
        .scaled(MinimalComplex::new(Ring(2), Ring(3)));
    let norm_squared = state.inner_product(&state).unwrap();
    assert_eq!(norm_squared, MinimalComplex::new(Ring(25), Ring(0)));
    assert_eq!(
        state.expectation(&operator).unwrap(),
        MinimalComplex::new(Ring(50), Ring(75))
    );
    assert_eq!(
        operator.adjoint().trace_in(state.space()).unwrap(),
        MinimalComplex::new(Ring(2), Ring(-3))
    );
    assert_eq!(
        operator
            .hs_inner_product_in(&operator, state.space())
            .unwrap(),
        MinimalComplex::new(Ring(13), Ring(0))
    );
}

#[test]
fn deep_composite_actions_and_operators_evaluate_and_drop_iteratively() {
    std::thread::Builder::new()
        .stack_size(128 * 1024)
        .spawn(|| {
            let counts = Counts::default();
            let h = H::factor(Index(0));
            let leaf = LocalState::primitive(Primitive(0, &counts)).at(&h).unwrap();
            let mut state = leaf.clone();
            let mut op = O::identity();
            for _ in 0..20_000 {
                state = state.act(&O::identity()).unwrap();
                op = (op.clone() * &O::identity()).unwrap().adjoint();
            }
            let result = state.evaluate().unwrap();
            assert!(result.expression.shares_node_with(&leaf.expression));
            assert!(
                leaf.act(&op)
                    .unwrap()
                    .evaluate()
                    .unwrap()
                    .expression
                    .shares_node_with(&leaf.expression)
            );
            drop(state);
            drop(op);
            drop(result);
            assert_eq!(counts.drops.get(), 0);
            assert_eq!(leaf.inner_product(&leaf).unwrap(), Complex::new(1.0, 0.0));

            let local = LocalOperator::primitive(Op).at(&h).unwrap();
            let mut nested = leaf.clone();
            let mut scaled = leaf.clone();
            let mut sum = leaf.clone();
            for _ in 0..20_000 {
                nested = S::tensor_in(&h, [S::unit(), nested]).unwrap();
                scaled = scaled.scaled(Complex::new(1.0, 0.0));
                sum = (sum.clone() + &S::zero(&h)).unwrap();
            }
            for state in [nested, scaled, sum] {
                let result = state.act(&local).unwrap().evaluate().unwrap();
                assert_eq!(
                    result.inner_product(&result).unwrap(),
                    Complex::new(1.0, 0.0)
                );
                drop(state);
            }
            assert_eq!(leaf.inner_product(&leaf).unwrap(), Complex::new(1.0, 0.0));
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn fallible_sums_retain_order_duplicates_and_explicit_zero_scaling() {
    let counts = Counts::default();
    let factor = H::factor(Index(0));
    let local = LocalState::primitive(Primitive(0, &counts));
    let first: S<'_> = local.at(&factor).unwrap();
    let sum = S::try_sum(
        &factor,
        [Ok(first.clone()), Ok(first.scaled(0.0)), Ok(first)],
    )
    .unwrap();
    let StateNode::Sum(parts) = sum.expression.node() else {
        panic!("expected sum")
    };
    assert_eq!(parts.len(), 3);
    assert!(matches!(parts.as_slice()[1].node(), StateNode::Scale(_)));
    assert_eq!(counts.acts.get(), 0);
    let first = O::<f64>::place(&LocalOperator::primitive(Op), &factor).unwrap();
    let sum = O::try_sum([Ok(first.clone()), Ok(first.scaled(0.0)), Ok(first)]).unwrap();
    let OpNode::Sum(parts) = sum.expression.node() else {
        panic!("expected sum")
    };
    assert_eq!(parts.len(), 3);
    assert!(matches!(parts.as_slice()[1].node(), OpNode::Scale(_)));
}
