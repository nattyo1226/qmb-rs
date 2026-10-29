use core::num::NonZeroUsize;
use qmb_core::{Dimension, DimensionOverflow, HilbertSpace, SpaceError};

// No Clone, Hash, Ord, Debug or generation implementation.
#[derive(PartialEq, Eq)]
struct Index(u32);
struct Structure<const D: usize>;
impl<const D: usize> Dimension for Structure<D> {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(D).unwrap())
    }
}

type Space = HilbertSpace<Index>;

#[test]
fn grouping_and_noncontiguous_subsystems_use_index_identity_and_order() {
    let a = Space::factor(Index(30));
    let b = Space::factor(Index(10));
    let c = Space::factor(Index(20));
    let ab = a.tensor(&b).unwrap();
    let left = ab.tensor(&c).unwrap();
    let right = a.tensor(&b.tensor(&c).unwrap()).unwrap();
    assert!(left.is_compatible_with(&right));
    assert_eq!(left.dim::<Structure<2>>(), Ok(8));
    assert_eq!(
        left.indices().map(|i| i.0).collect::<Vec<_>>(),
        vec![30, 10, 20]
    );
    let children = left.children().collect::<Vec<_>>();
    assert_eq!(children.len(), 2);
    assert!(children[0].shares_node_with(&ab));
    assert!(children[1].shares_node_with(&c));
    assert!(a.tensor(&c).unwrap().is_subsystem_of(&left));
    assert!(!c.tensor(&a).unwrap().is_subsystem_of(&left));
    assert!(!Space::factor(Index(99)).is_subsystem_of(&left));
    assert!(left.is_subsystem_of(&left));
    assert!(Space::unit().is_subsystem_of(&left));
    assert!(!a.is_subsystem_of(&Space::unit()));
    assert!(!ab.is_compatible_with(&b.tensor(&a).unwrap()));
    assert!(a.is_compatible_with(&Space::factor(Index(30))));
}

#[test]
fn duplicates_across_nested_and_shared_subspaces_are_rejected() {
    let a = Space::factor(Index(1));
    let b = Space::factor(Index(2));
    let ab = a.tensor(&b).unwrap();
    assert!(matches!(ab.tensor(&a), Err(SpaceError::DuplicateIndex)));
    assert!(matches!(ab.tensor(&ab), Err(SpaceError::DuplicateIndex)));
    assert!(matches!(
        Space::tensor_all([a.clone(), Space::factor(Index(1))]),
        Err(SpaceError::DuplicateIndex)
    ));
}

#[test]
fn empty_singleton_and_indexed_one_dimensional_factor_keep_their_meanings() {
    let empty = Space::tensor_all([]).unwrap();
    assert_eq!(empty.dim::<Structure<2>>(), Ok(1));
    assert_eq!(empty.indices().count(), 0);
    assert_eq!(empty.children().len(), 0);
    let a = Space::factor(Index(1));
    let singleton = Space::tensor_all([a.clone()]).unwrap();
    assert_eq!(a.children().len(), 0);
    assert_eq!(singleton.children().len(), 1);
    assert!(a.is_compatible_with(&singleton));
    assert!(a.is_compatible_with(&Space::tensor_all([empty, a.clone()]).unwrap()));
    let unit_factor = HilbertSpace::<_>::factor(Index(1));
    assert_eq!(unit_factor.dim::<Structure<1>>(), Ok(1));
    assert!(!unit_factor.is_compatible_with(&HilbertSpace::unit()));
}

#[test]
fn large_dimension_is_an_error_only_when_queried() {
    let h = Space::tensor_all((0..usize::BITS + 1).map(|n| Space::factor(Index(n)))).unwrap();
    assert_eq!(h.dim::<Structure<2>>(), Err(DimensionOverflow));
    assert_eq!(h.indices().count(), usize::BITS as usize + 1);
    assert!(Space::factor(Index(0)).is_subsystem_of(&h));
}

#[test]
fn deep_spaces_walk_and_drop_on_a_small_stack_with_shared_leaves() {
    std::thread::Builder::new()
        .stack_size(128 * 1024)
        .spawn(|| {
            let leaf = Space::factor(Index(8));
            let mut h = leaf.clone();
            for _ in 0..20_000 {
                h = Space::tensor_all([h]).unwrap();
            }
            assert_eq!(h.dim::<Structure<2>>(), Ok(2));
            assert!(h.is_compatible_with(&leaf));
            assert_eq!(h.indices().next().unwrap().0, 8);
            drop(h);
            assert_eq!(leaf.factor_index().unwrap().0, 8);
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn leaves_borrow_actual_handles_in_order_and_remain_exhausted() {
    // Reading leaves needs no Eq, Clone, Default or numeric capabilities.
    struct Opaque;
    let leaf = HilbertSpace::factor(Opaque);
    let mut leaves = leaf.leaves();
    assert!(std::ptr::eq(leaves.next().unwrap(), &leaf));
    assert!(leaves.next().is_none());
    assert!(leaves.next().is_none());
    let unit = HilbertSpace::<Opaque>::unit();
    let mut empty = unit.leaves();
    assert!(empty.next().is_none());
    assert!(empty.next().is_none());
}
