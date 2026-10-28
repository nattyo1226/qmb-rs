//! Primitive capabilities and index structure are independent.

use core::num::NonZeroUsize;
use qmb_core::{
    Dimension, DimensionOverflow, HilbertSpace, Identity, LocalOperator, LocalState,
    operator_basis_len, state_basis_len,
};

struct Ket;
struct Op;
struct OtherKet;
struct OtherOp;
struct Opaque;

impl Identity for Op {
    fn identity() -> Self {
        Self
    }
}
impl Dimension for Ket {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(3).unwrap())
    }
}
impl Dimension for Op {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(9).unwrap())
    }
}

#[test]
fn unrelated_primitive_types_share_and_bind_the_same_index_structure() {
    let factor = HilbertSpace::factor(7u8);
    // Neither coefficient nor primitives implement numeric rules or Clone.
    let first = LocalState::<Ket, Op, Opaque>::primitive(Ket)
        .at(&factor)
        .unwrap();
    let second = LocalState::<OtherKet, OtherOp, Opaque>::primitive(OtherKet)
        .at(&factor)
        .unwrap();
    assert!(first.space().shares_node_with(second.space()));
    let op = LocalOperator::<Op, Opaque>::primitive(Op::identity())
        .at(&factor)
        .unwrap()
        .bind(&factor)
        .unwrap();
    let other = LocalOperator::<OtherOp, Opaque>::primitive(OtherOp)
        .at(&factor)
        .unwrap()
        .bind(&factor)
        .unwrap();
    assert!(op.space().unwrap().shares_node_with(other.space().unwrap()));
    assert!(first.act(&op).unwrap().space().shares_node_with(&factor));
    assert!(
        second
            .act(&other)
            .unwrap()
            .space()
            .shares_node_with(&factor)
    );
}

#[test]
fn span_queries_distinguish_state_and_operator_spaces_without_squaring_twice() {
    assert_eq!(state_basis_len::<Ket>(), Ok(3));
    assert_eq!(operator_basis_len::<Op>(), Ok(9));
    let space = HilbertSpace::tensor([HilbertSpace::factor(1u8), HilbertSpace::factor(2)]).unwrap();
    assert_eq!(space.dimension::<Ket>(), Ok(9));
    assert_eq!(space.dimension::<Op>(), Ok(81));
    // Querying a dimension does not attach its primitive type to the space.
    assert_eq!(space.dimension_with(NonZeroUsize::new(2).unwrap()), Ok(4));
    assert_eq!(space.indices().copied().collect::<Vec<_>>(), [1, 2]);
}

struct Overflow;
impl Dimension for Overflow {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Err(DimensionOverflow)
    }
}

#[test]
fn local_dimension_errors_are_requested_only_for_queries_on_actual_factors() {
    let space = HilbertSpace::factor(0u8);
    let state = LocalState::<Overflow, Op, Opaque>::zero()
        .at(&space)
        .unwrap();
    assert!(state.space().shares_node_with(&space));
    assert_eq!(state_basis_len::<Overflow>(), Err(DimensionOverflow));
    assert_eq!(operator_basis_len::<Overflow>(), Err(DimensionOverflow));
    assert_eq!(space.dimension::<Overflow>(), Err(DimensionOverflow));
    assert_eq!(HilbertSpace::<u8>::unit().dimension::<Overflow>(), Ok(1));
}

#[test]
fn identities_share_a_trait_without_numeric_or_unrelated_primitive_bounds() {
    fn identity<T: qmb_core::ops::Identity>() -> T {
        T::identity()
    }

    let _: Op = identity();
    let local: LocalOperator<Op, Opaque> = identity();
    let factor = HilbertSpace::factor(0u8);
    assert!(local.at(&factor).is_ok());
    let _: LocalOperator<Op, Opaque> = LocalOperator::product_all([]);

    // Neither index, operator primitive nor coefficient has any capabilities.
    let unbound: qmb_core::Operator<Opaque, OtherOp, Opaque> = identity();
    assert!(unbound.space().is_none());
    let composite: qmb_core::Operator<u8, OtherOp, Opaque> = identity();
    assert!(
        composite
            .bind(&factor)
            .unwrap()
            .space()
            .unwrap()
            .shares_node_with(&factor)
    );
}

#[test]
fn moved_common_parts_and_dimension_capabilities_keep_their_public_paths() {
    fn dimension<P: qmb_core::ops::Dimension>() -> usize {
        P::dim().unwrap().get()
    }
    assert_eq!(dimension::<Ket>(), 3);
    assert_eq!(dimension::<Op>(), 9);
    let terms = qmb_core::expr::parts::LinearCombination::term(Opaque, OtherOp);
    let _: qmb_core::LinearCombination<OtherOp, Opaque> = terms;
}

#[test]
fn local_identity_uses_primitive_rules_without_dimension_queries() {
    use qmb_core::{Act, Adjoint, InnerProduct, LinearCombination as Lc, Product, Trace};
    use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

    static ACTIONS: AtomicUsize = AtomicUsize::new(0);
    static PRODUCTS: AtomicUsize = AtomicUsize::new(0);
    static ADJOINTS: AtomicUsize = AtomicUsize::new(0);
    struct UnitState;
    struct UnitOperator;
    impl Identity for UnitOperator {
        fn identity() -> Self {
            Self
        }
    }
    impl Act<f64> for UnitState {
        type Operator = UnitOperator;
        type Term = Self;
        fn act(&self, _: &UnitOperator) -> Lc<Self, f64> {
            ACTIONS.fetch_add(1, Relaxed);
            Lc::primitive(Self)
        }
    }
    impl InnerProduct<f64> for UnitState {
        fn inner_product(&self, _: &Self) -> f64 {
            1.0
        }
    }
    impl Product<f64> for UnitOperator {
        type Term = Self;
        fn product(&self, _: &Self) -> Lc<Self, f64> {
            PRODUCTS.fetch_add(1, Relaxed);
            Lc::primitive(Self)
        }
    }
    impl Adjoint<f64> for UnitOperator {
        type Term = Self;
        fn adjoint(&self) -> Lc<Self, f64> {
            ADJOINTS.fetch_add(1, Relaxed);
            Lc::primitive(Self)
        }
    }
    impl Trace<f64> for UnitOperator {
        fn trace(&self) -> f64 {
            1.0
        }
    }

    let identity = LocalOperator::<UnitOperator, f64>::identity();
    let state = LocalState::primitive(UnitState);
    assert_eq!(state.act(&identity).inner_product(&state), 1.0);
    assert_eq!(ACTIONS.load(Relaxed), 1);
    assert_eq!(identity.product(&identity).trace(), 1.0);
    assert_eq!(PRODUCTS.load(Relaxed), 1);
    assert_eq!(identity.adjoint().trace(), 1.0);
    assert_eq!(ADJOINTS.load(Relaxed), 1);
}
