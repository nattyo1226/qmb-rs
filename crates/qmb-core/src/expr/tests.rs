use alloc::{vec, vec::Vec};
use core::{cell::Cell, num::NonZeroUsize};
use std::{sync::Arc, thread};

use super::{
    LocalOperator, LocalState,
    local::{OperatorNode, StateNode},
    parts::*,
};
use crate::{Act, Adjoint, Complex, LinearCombination, LocalStructure};

#[derive(Default)]
struct Counts {
    actions: Cell<usize>,
    adjoints: Cell<usize>,
    drops: Cell<usize>,
}
// Neither primitives nor metadata are Clone. Primitives borrow non-Send data.
struct State<'a>(u8, &'a Counts);
impl Drop for State<'_> {
    fn drop(&mut self) {
        self.1.drops.set(self.1.drops.get() + 1)
    }
}
struct Operator<'a>(u8, &'a Counts);
struct Structure<'a>(core::marker::PhantomData<&'a Counts>);
impl<'a> LocalStructure for Structure<'a> {
    type State = State<'a>;
    type Operator = Operator<'a>;
    const DIMENSION: NonZeroUsize = NonZeroUsize::new(2).unwrap();
}
// This fixture tests traversal, not algebraic laws. There is deliberately no
// Product or Identity implementation. Matrix-oracle integration tests check laws.
impl<'a> Act<f64> for State<'a> {
    type Operator = Operator<'a>;
    type Term = Self;
    fn act(&self, op: &Operator<'a>) -> LinearCombination<Self, f64> {
        self.1.actions.set(self.1.actions.get() + 1);
        let label = match op.0 {
            0 => self.0,
            1 => self.0 + 1,
            _ => self.0 * 2,
        };
        LinearCombination::primitive(State(label, self.1))
    }
}
impl<'a> Adjoint<f64> for Operator<'a> {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, f64> {
        self.1.adjoints.set(self.1.adjoints.get() + 1);
        LinearCombination::primitive(Operator(self.0, self.1))
    }
}
type S<'a, T = f64> = LocalState<Structure<'a>, T>;
type O<'a, T = f64> = LocalOperator<Structure<'a>, T>;
fn coefficients(state: &S<'_>) -> Vec<(u8, Complex<f64>)> {
    let mut stack = vec![(state, Complex::new(1.0, 0.0))];
    let mut result = Vec::new();
    while let Some((state, c)) = stack.pop() {
        match state.node() {
            StateNode::Primitive(p) => result.push((p.0, c)),
            StateNode::Scale(s) => stack.push((&s.operand, c * s.coefficient)),
            StateNode::Sum(s) => {
                for child in s.iter().rev() {
                    stack.push((child, c));
                }
            }
            StateNode::Act(_) => panic!("unevaluated action in result"),
        }
    }
    result
}

#[test]
fn construction_shares_without_numeric_or_clone_bounds() {
    struct Scalar;
    let _: SumExpr<Scalar> = SumExpr::default();
    let counts = Counts::default();
    let state = S::<Scalar>::primitive(State(1, &counts));
    let op = O::<Scalar>::primitive(Operator(0, &counts));
    let sum = &state + &state;
    let StateNode::Sum(children) = sum.node() else {
        panic!("sum")
    };
    for child in children {
        assert!(Arc::ptr_eq(
            state.root.as_ref().unwrap(),
            child.root.as_ref().unwrap()
        ));
    }
    let scaled = state.scaled(Complex::new(Scalar, Scalar));
    let StateNode::Scale(s) = scaled.node() else {
        panic!("scale")
    };
    assert!(Arc::ptr_eq(
        state.root.as_ref().unwrap(),
        s.operand.root.as_ref().unwrap()
    ));
    let pending = state.act(&op);
    assert!(pending.has_action);
    assert_eq!(counts.actions.get(), 0);
    let op_clone = op.clone();
    let product = &op * &op_clone;
    let OperatorNode::Product(p) = product.node() else {
        panic!("product")
    };
    assert!(Arc::ptr_eq(
        p.lhs.root.as_ref().unwrap(),
        p.rhs.root.as_ref().unwrap()
    ));
    let _ = S::<Scalar>::from_linear_combination(LinearCombination::term(
        Complex::new(Scalar, Scalar),
        state.clone(),
    ));
    let _ = O::<Scalar>::identity().adjoint();
    let unrelated = S::<Scalar>::primitive(State(1, &counts));
    assert!(!Arc::ptr_eq(
        state.root.as_ref().unwrap(),
        unrelated.root.as_ref().unwrap()
    ));
}

#[test]
fn scalar_syntax_moves_opaque_coefficients_and_shares_expression_nodes() {
    struct Scalar(u8);
    let coefficient = || Complex::new(Scalar(2), Scalar(3));
    let counts = Counts::default();
    let state = S::<Scalar>::primitive(State(1, &counts));
    let operator = O::<Scalar>::primitive(Operator(0, &counts));
    for scaled in [
        state.clone() * coefficient(),
        &state * coefficient(),
        coefficient() * state.clone(),
        coefficient() * &state,
    ] {
        let StateNode::Scale(s) = scaled.node() else {
            panic!("scale")
        };
        assert_eq!((s.coefficient.re.0, s.coefficient.im.0), (2, 3));
        assert!(Arc::ptr_eq(
            s.operand.root.as_ref().unwrap(),
            state.root.as_ref().unwrap()
        ));
    }
    for scaled in [
        operator.clone() * coefficient(),
        &operator * coefficient(),
        coefficient() * operator.clone(),
        coefficient() * &operator,
    ] {
        let OperatorNode::Scale(s) = scaled.node() else {
            panic!("scale")
        };
        assert_eq!((s.coefficient.re.0, s.coefficient.im.0), (2, 3));
        assert!(Arc::ptr_eq(
            s.operand.root.as_ref().unwrap(),
            operator.root.as_ref().unwrap()
        ));
    }
    assert_eq!(counts.actions.get(), 0);
    assert_eq!(counts.drops.get(), 0);
}

#[test]
fn constructor_cardinality_sharing_and_single_final_drop() {
    let counts = Counts::default();
    let leaf = S::<f64>::primitive(State(1, &counts));
    let single_state = S::sum([leaf.clone()]);
    assert!(Arc::ptr_eq(
        leaf.root.as_ref().unwrap(),
        single_state.root.as_ref().unwrap()
    ));
    let ready = single_state.evaluate();
    assert!(Arc::ptr_eq(
        leaf.root.as_ref().unwrap(),
        ready.root.as_ref().unwrap()
    ));
    assert!(matches!(
        S::<f64>::zero().node(),
        StateNode::Sum(SumExpr::Empty)
    ));
    assert!(matches!(
        O::<f64>::product_all([]).node(),
        OperatorNode::Identity
    ));
    let op = O::<f64>::primitive(Operator(0, &counts));
    let single = O::product_all([op.clone()]);
    assert!(Arc::ptr_eq(
        op.root.as_ref().unwrap(),
        single.root.as_ref().unwrap()
    ));
    let shared = &leaf + &leaf;
    drop(leaf);
    drop(ready);
    drop(single);
    drop(single_state);
    assert_eq!(counts.drops.get(), 0);
    drop(shared);
    assert_eq!(counts.drops.get(), 1);
}

#[test]
fn action_order_no_product_and_immutable_reusable_results() {
    let counts = Counts::default();
    let state = S::primitive(State(1, &counts));
    let a = O::primitive(Operator(2, &counts));
    let b = O::primitive(Operator(1, &counts));
    let chain = state.act(&b).act(&a);
    assert_eq!(counts.actions.get(), 0);
    let result = chain.evaluate();
    assert_eq!(coefficients(&result), vec![(4, Complex::new(1.0, 0.0))]);
    assert_eq!(counts.actions.get(), 2);
    assert!(chain.has_action);
    assert!(!result.has_action);
    assert_eq!(
        coefficients(&state.act(&a.product(&b)).evaluate()),
        coefficients(&result)
    );
    let again = result.evaluate();
    assert!(Arc::ptr_eq(
        result.root.as_ref().unwrap(),
        again.root.as_ref().unwrap()
    ));
    assert_eq!(
        coefficients(&again.act(&b).evaluate()),
        vec![(5, Complex::new(1.0, 0.0))]
    );
    assert_eq!(coefficients(&state), vec![(1, Complex::new(1.0, 0.0))]);
}

#[test]
fn ordered_zero_and_duplicate_terms_survive_mapping_and_evaluation() {
    let counts = Counts::default();
    let terms = LinearCombination::term(Complex::new(0.0, 0.0), State(2, &counts))
        + LinearCombination::term(Complex::new(3.0, 4.0), State(1, &counts))
        + LinearCombination::term(Complex::new(1.0, 0.0), State(1, &counts));
    let state = S::from_linear_combination(terms.map(S::primitive));
    let result = state.act(&O::primitive(Operator(0, &counts))).evaluate();
    assert_eq!(
        coefficients(&result),
        vec![
            (2, Complex::new(0.0, 0.0)),
            (1, Complex::new(3.0, 4.0)),
            (1, Complex::new(1.0, 0.0))
        ]
    );
    assert_eq!(counts.actions.get(), 3);
    let zero = S::zero()
        .act(&O::primitive(Operator(0, &counts)))
        .evaluate();
    assert!(coefficients(&zero).is_empty());
    assert_eq!(counts.actions.get(), 3);
    assert!(!zero.has_action);
}

#[test]
fn operator_sum_and_state_sum_preserve_left_to_right_term_order() {
    let counts = Counts::default();
    let state = S::sum([
        S::primitive(State(1, &counts)),
        S::primitive(State(2, &counts)),
    ]);
    let operator = O::sum([
        O::primitive(Operator(1, &counts)),
        O::primitive(Operator(2, &counts)),
    ]);
    let result = state.act(&operator).evaluate();
    assert_eq!(
        coefficients(&result)
            .into_iter()
            .map(|(label, _)| label)
            .collect::<Vec<_>>(),
        vec![2, 3, 2, 4]
    );
    assert_eq!(counts.actions.get(), 4);
}

#[test]
fn deep_action_product_scale_sum_and_adjoint_chains_use_work_stacks() {
    thread::Builder::new()
        .stack_size(128 * 1024)
        .spawn(|| {
            let counts = Counts::default();
            let leaf = S::primitive(State(1, &counts));
            let mut chain = leaf.clone();
            let mut product = O::identity();
            let mut adjoint = O::identity();
            let mut scaled = leaf.clone();
            let mut sum = leaf.clone();
            for _ in 0..20_000 {
                chain = chain.act(&O::identity());
                product = product.product(&O::identity());
                adjoint = adjoint.adjoint();
                scaled = scaled.scaled(Complex::new(1.0, 0.0));
                sum = &sum + S::zero();
            }
            for state in [
                chain.evaluate(),
                leaf.act(&product).evaluate(),
                leaf.act(&adjoint).evaluate(),
                scaled.act(&O::primitive(Operator(0, &counts))).evaluate(),
                sum.act(&O::identity()).evaluate(),
            ] {
                assert_eq!(coefficients(&state), vec![(1, Complex::new(1.0, 0.0))]);
            }
            drop(chain);
            drop(product);
            drop(adjoint);
            drop(scaled);
            drop(sum);
            assert_eq!(coefficients(&leaf), vec![(1, Complex::new(1.0, 0.0))]);
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn send_and_sync_follow_stored_types_without_metadata_bounds() {
    struct Metadata {
        _not_send: std::rc::Rc<()>,
    }
    impl LocalStructure for Metadata {
        type State = u8;
        type Operator = u8;
        const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
    }
    fn send_sync<T: Send + Sync>() {}
    send_sync::<LocalState<Metadata, f64>>();
    send_sync::<LocalOperator<Metadata, f64>>();
}

#[test]
fn delayed_evaluation_needs_no_coefficient_arithmetic() {
    #[derive(Clone)]
    struct Scalar(i32);
    impl core::ops::Neg for Scalar {
        type Output = Self;
        fn neg(self) -> Self {
            Self(-self.0)
        }
    }
    impl<'a> Act<Scalar> for State<'a> {
        type Operator = Operator<'a>;
        type Term = Self;
        fn act(&self, _: &Self::Operator) -> LinearCombination<Self, Scalar> {
            LinearCombination::term(Complex::new(Scalar(2), Scalar(0)), State(self.0, self.1))
        }
    }
    impl Adjoint<Scalar> for Operator<'_> {
        type Term = Self;
        fn adjoint(&self) -> LinearCombination<Self, Scalar> {
            LinearCombination::term(Complex::new(Scalar(1), Scalar(0)), Operator(self.0, self.1))
        }
    }
    let counts = Counts::default();
    let state = S::<Scalar>::primitive(State(1, &counts));
    let op = O::<Scalar>::primitive(Operator(0, &counts));
    let result = state.act(&op).evaluate();
    let StateNode::Scale(s) = result.node() else {
        panic!("rule scale")
    };
    assert_eq!(s.coefficient.re.0, 2);
    assert!(matches!(s.operand.node(), StateNode::Primitive(_)));
    assert!(!result.has_action);
}

#[test]
fn scalar_reductions_need_no_num_division_remainder_one_or_copy() {
    use crate::{Identity, InnerProduct, Product, Trace};
    use core::ops::{Add, Mul, Neg, Sub};
    use num_traits::Zero;
    #[derive(Clone)]
    struct Ring(i32);
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
    struct Primitive;
    struct RingStructure;
    impl LocalStructure for RingStructure {
        type State = Primitive;
        type Operator = Primitive;
        const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
    }
    fn one() -> LinearCombination<Primitive, Ring> {
        LinearCombination::term(Complex::new(Ring(1), Ring(0)), Primitive)
    }
    impl Identity<Ring> for RingStructure {
        fn identity() -> LinearCombination<Primitive, Ring> {
            one()
        }
    }
    impl Adjoint<Ring> for Primitive {
        type Term = Self;
        fn adjoint(&self) -> LinearCombination<Self, Ring> {
            one()
        }
    }
    impl Product<Ring> for Primitive {
        type Term = Self;
        fn product(&self, _: &Self) -> LinearCombination<Self, Ring> {
            one()
        }
    }
    impl Act<Ring> for Primitive {
        type Operator = Self;
        type Term = Self;
        fn act(&self, _: &Self) -> LinearCombination<Self, Ring> {
            one()
        }
    }
    impl Trace<Ring> for Primitive {
        fn trace(&self) -> Complex<Ring> {
            Complex::new(Ring(1), Ring(0))
        }
    }
    impl InnerProduct<Ring> for Primitive {
        fn inner_product(&self, _: &Self) -> Complex<Ring> {
            Complex::new(Ring(1), Ring(0))
        }
    }
    let op =
        LocalOperator::<RingStructure, Ring>::identity().scaled(Complex::new(Ring(3), Ring(4)));
    let trace = op
        .adjoint()
        .product(&LocalOperator::primitive(Primitive))
        .trace();
    assert_eq!((trace.re.0, trace.im.0), (3, -4));
    let state = LocalState::<RingStructure, Ring>::primitive(Primitive);
    assert_eq!(state.norm_squared().0, 1);
    let expectation = state.expectation(&op);
    assert_eq!((expectation.re.0, expectation.im.0), (3, 4));
}
