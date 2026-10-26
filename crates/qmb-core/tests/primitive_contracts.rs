//! Primitive capabilities and index structure are independent.

use core::num::NonZeroUsize;
use qmb_core::{
    DimensionOverflow, HilbertSpace, Identity, LocalOperator, LocalState, SpanDimension,
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
impl SpanDimension for Ket {
    fn span_dimension() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(3).unwrap())
    }
}
impl SpanDimension for Op {
    fn span_dimension() -> Result<NonZeroUsize, DimensionOverflow> {
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
impl SpanDimension for Overflow {
    fn span_dimension() -> Result<NonZeroUsize, DimensionOverflow> {
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
