use super::composite::{OpNode, StateNode};
use crate::{
    Act, Adjoint, Complex, ExpressionError as E, HilbertSpace, InnerProduct,
    LinearCombination as Lc, LocalOperator, LocalState, LocalStructure, Operator, State,
};
use core::{cell::Cell, marker::PhantomData, num::NonZeroUsize};

#[derive(PartialEq, Eq)]
struct Index(u8);
#[derive(Default)]
struct Counts {
    acts: Cell<usize>,
    drops: Cell<usize>,
}
struct Primitive<'a>(u8, &'a Counts);
impl Drop for Primitive<'_> {
    fn drop(&mut self) {
        self.1.drops.set(self.1.drops.get() + 1);
    }
}
struct Op;
struct L<'a>(PhantomData<&'a Counts>);
impl<'a> LocalStructure for L<'a> {
    type State = Primitive<'a>;
    type Operator = Op;
    const DIMENSION: NonZeroUsize = NonZeroUsize::new(2).unwrap();
}
// No Product or Identity capability. Metadata, indices and primitives are not Clone.
impl<'a> Act<f64> for Primitive<'a> {
    type Operator = Op;
    type Term = Self;
    fn act(&self, _: &Op) -> Lc<Self, f64> {
        self.1.acts.set(self.1.acts.get() + 1);
        Lc::primitive(Primitive(1 - self.0, self.1))
    }
}
impl Adjoint<f64> for Op {
    type Term = Self;
    fn adjoint(&self) -> Lc<Self, f64> {
        Lc::primitive(Op)
    }
}
impl InnerProduct<f64> for Primitive<'_> {
    fn inner_product(&self, rhs: &Self) -> Complex<f64> {
        Complex::new(if self.0 == rhs.0 { 1.0 } else { 0.0 }, 0.0)
    }
}
type H<'a> = HilbertSpace<Index, L<'a>>;
type S<'a, T = f64> = State<Index, L<'a>, T>;
type O<'a, T = f64> = Operator<Index, L<'a>, T>;

#[test]
fn construction_uses_opaque_scalars_and_shared_nonclone_borrowed_leaves() {
    struct Scalar;
    let counts = Counts::default();
    let h = H::factor(Index(0));
    let local = LocalState::<L<'_>, Scalar>::primitive(Primitive(0, &counts));
    let state = local.at(&h).unwrap();
    let op = LocalOperator::primitive(Op).at(&h).unwrap();
    let a = state
        .add(&state)
        .unwrap()
        .scaled(Complex::new(Scalar, Scalar));
    let action = a.act(&op).unwrap();
    let clone = action.clone();
    assert!(clone.space().shares_node_with(&h));
    let _ = S::<Scalar>::tensor([]).unwrap();
    let _ = S::<Scalar>::zero(&h);
    let _ = op.product(&op).unwrap().adjoint();
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
    let state = left.tensor_product(&right).unwrap();
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
    assert!(before.as_slice()[1].same_node(&after.as_slice()[1]));
    assert!(!state.expression.pending_action && action.expression.pending_action);
    assert!(!result.expression.pending_action);
    assert!(
        result
            .evaluate()
            .unwrap()
            .expression
            .same_node(&result.expression)
    );
    let identity = O::identity();
    assert!(
        state
            .act(&identity)
            .unwrap()
            .evaluate()
            .unwrap()
            .expression
            .same_node(&state.expression)
    );
    let x2 = op.product(&op).unwrap();
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
    let ab = a.tensor_product(&b).unwrap();
    let op = LocalOperator::<L<'_>, f64>::primitive(Op).at(&a).unwrap();
    let small = op.bind(&a).unwrap();
    let large = op.bind(&ab).unwrap();
    assert!(matches!(small.add(&large), Err(E::IncompatibleSpaces)));
    assert!(matches!(small.product(&large), Err(E::IncompatibleSpaces)));
    // The third, unbound operand suppresses comparisons of the first two.
    let mixed = O::sum([small.clone(), large.clone(), O::identity()]).unwrap();
    assert!(mixed.space().is_none());
    assert!(
        O::product_all([small.clone(), large.clone(), O::identity()])
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
    let state = left.tensor_product(&right).unwrap();
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
        LocalOperator::<L<'_>, f64>::zero().at(&a).unwrap().bind(&b),
        Err(E::MissingIndex)
    ));
    assert!(matches!(
        LocalOperator::<L<'_>, f64>::identity()
            .at(&a)
            .unwrap()
            .bind(&b),
        Err(E::MissingIndex)
    ));
    assert!(matches!(
        LocalState::<L<'_>, f64>::zero().at(&ab),
        Err(E::NotAFactor)
    ));
    assert!(matches!(
        S::<f64>::zero(&a).add(&S::zero(&b)),
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
    impl core::ops::Neg for Scalar {
        type Output = Self;
        fn neg(self) -> Self {
            self
        }
    }
    impl<'a> Act<Scalar> for Primitive<'a> {
        type Operator = Op;
        type Term = Self;
        fn act(&self, _: &Op) -> Lc<Self, Scalar> {
            Lc::term(Complex::new(Scalar, Scalar), Primitive(self.0, self.1))
        }
    }
    impl Adjoint<Scalar> for Op {
        type Term = Self;
        fn adjoint(&self) -> Lc<Self, Scalar> {
            Lc::term(Complex::new(Scalar, Scalar), Op)
        }
    }
    let counts = Counts::default();
    let factor = H::factor(Index(1));
    let state = LocalState::<L<'_>, Scalar>::primitive(Primitive(0, &counts))
        .at(&factor)
        .unwrap();
    let op = LocalOperator::primitive(Op).at(&factor).unwrap();
    assert!(state.act(&op).unwrap().evaluate().is_ok());
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
        assert!(parts.iter().all(|p| p.same_node(&original)));
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
    assert!(state.expression.same_node(&original));
}

#[test]
fn unit_zero_singleton_and_placement_preserve_spaces_and_lazy_local_expressions() {
    let counts = Counts::default();
    let factor = H::factor(Index(1));
    let unit = S::<f64>::unit();
    let leaf = LocalState::primitive(Primitive(0, &counts))
        .at(&factor)
        .unwrap();
    let nested = unit
        .tensor_product(&S::tensor([leaf.clone()]).unwrap())
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
            .same_node(&leaf.expression)
    );
    assert!(matches!(
        leaf.tensor_product(&leaf),
        Err(E::Space(crate::SpaceError::DuplicateIndex))
    ));
    assert!(matches!(
        O::<f64>::place(&LocalOperator::identity(), unit.space()),
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
fn shared_composite_types_follow_stored_data_without_metadata_bounds() {
    struct Metadata(std::rc::Rc<()>);
    impl LocalStructure for Metadata {
        type State = ();
        type Operator = ();
        const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
    }
    fn send_sync<T: Send + Sync>() {}
    send_sync::<HilbertSpace<u8, Metadata>>();
    send_sync::<State<u8, Metadata, f64>>();
    send_sync::<Operator<u8, Metadata, f64>>();
    // The metadata value is not stored by any of these types.
    let metadata = Metadata(std::rc::Rc::new(()));
    assert_eq!(std::rc::Rc::strong_count(&metadata.0), 1);
}

#[test]
fn composite_scalar_reductions_require_neither_num_float_copy_nor_division() {
    use crate::{
        Identity, Product, Trace,
        num_traits::{One, Zero},
    };
    use core::ops::{Add, Mul, Neg, Sub};

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
    struct Point;
    struct Primitive;
    impl LocalStructure for Point {
        type State = Primitive;
        type Operator = Primitive;
        const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
    }
    impl Identity<Ring> for Point {
        fn identity() -> Lc<Primitive, Ring> {
            Lc::primitive(Primitive)
        }
    }
    impl Act<Ring> for Primitive {
        type Operator = Self;
        type Term = Self;
        fn act(&self, _: &Self) -> Lc<Self, Ring> {
            Lc::primitive(Primitive)
        }
    }
    impl Product<Ring> for Primitive {
        type Term = Self;
        fn product(&self, _: &Self) -> Lc<Self, Ring> {
            Lc::primitive(Primitive)
        }
    }
    impl Adjoint<Ring> for Primitive {
        type Term = Self;
        fn adjoint(&self) -> Lc<Self, Ring> {
            Lc::primitive(Primitive)
        }
    }
    impl InnerProduct<Ring> for Primitive {
        fn inner_product(&self, _: &Self) -> Complex<Ring> {
            Complex::new(Ring(1), Ring(0))
        }
    }
    impl Trace<Ring> for Primitive {
        fn trace(&self) -> Complex<Ring> {
            Complex::new(Ring(1), Ring(0))
        }
    }
    let a = HilbertSpace::<_, Point>::factor(0u8);
    let b = HilbertSpace::factor(1);
    let state = LocalState::<Point, Ring>::primitive(Primitive)
        .at(&a)
        .unwrap()
        .tensor_product(&LocalState::primitive(Primitive).at(&b).unwrap())
        .unwrap()
        .scaled(Complex::new(Ring(3), Ring(4)));
    let operator = LocalOperator::primitive(Primitive)
        .at(&a)
        .unwrap()
        .scaled(Complex::new(Ring(2), Ring(3)));
    let norm_squared = state.inner_product(&state).unwrap();
    assert_eq!(norm_squared, Complex::new(Ring(25), Ring(0)));
    assert_eq!(
        state.expectation(&operator).unwrap(),
        Complex::new(Ring(50), Ring(75))
    );
    assert_eq!(
        operator.adjoint().trace_in(state.space()).unwrap(),
        Complex::new(Ring(2), Ring(-3))
    );
    assert_eq!(
        operator
            .hilbert_schmidt_inner_product_in(&operator, state.space())
            .unwrap(),
        Complex::new(Ring(13), Ring(0))
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
                op = op.product(&O::identity()).unwrap().adjoint();
            }
            let result = state.evaluate().unwrap();
            assert!(result.expression.same_node(&leaf.expression));
            assert!(
                leaf.act(&op)
                    .unwrap()
                    .evaluate()
                    .unwrap()
                    .expression
                    .same_node(&leaf.expression)
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
                nested = S::tensor([nested]).unwrap();
                scaled = scaled.scaled(Complex::new(1.0, 0.0));
                sum = sum.add(&S::zero(&h)).unwrap();
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
