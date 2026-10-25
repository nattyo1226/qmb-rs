use core::num::NonZeroUsize;
use qmb_core::{DimensionOverflow, HilbertSpace, LocalStructure, SpaceError};

// No Clone, Hash, Ord, Debug or generation implementation.
#[derive(PartialEq, Eq)]
struct Index(u32);
struct Structure<const D: usize>;
impl<const D: usize> LocalStructure for Structure<D> {
    type State = ();
    type Operator = ();
    const DIMENSION: NonZeroUsize = NonZeroUsize::new(D).unwrap();
}
type Space = HilbertSpace<Index, Structure<2>>;

#[test]
fn grouping_and_noncontiguous_subsystems_use_index_identity_and_order() {
    let a = Space::factor(Index(30));
    let b = Space::factor(Index(10));
    let c = Space::factor(Index(20));
    let ab = a.tensor_product(&b).unwrap();
    let left = ab.tensor_product(&c).unwrap();
    let right = a.tensor_product(&b.tensor_product(&c).unwrap()).unwrap();
    assert!(left.is_compatible_with(&right));
    assert_eq!(left.dimension(), Ok(8));
    assert_eq!(
        left.indices().map(|i| i.0).collect::<Vec<_>>(),
        vec![30, 10, 20]
    );
    let children = left.children().collect::<Vec<_>>();
    assert_eq!(children.len(), 2);
    assert!(children[0].shares_node_with(&ab));
    assert!(children[1].shares_node_with(&c));
    assert!(a.tensor_product(&c).unwrap().is_subsystem_of(&left));
    assert!(!c.tensor_product(&a).unwrap().is_subsystem_of(&left));
    assert!(!Space::factor(Index(99)).is_subsystem_of(&left));
    assert!(left.is_subsystem_of(&left));
    assert!(Space::unit().is_subsystem_of(&left));
    assert!(!a.is_subsystem_of(&Space::unit()));
    assert!(!ab.is_compatible_with(&b.tensor_product(&a).unwrap()));
    assert!(a.is_compatible_with(&Space::factor(Index(30))));
}

#[test]
fn duplicates_across_nested_and_shared_subspaces_are_rejected() {
    let a = Space::factor(Index(1));
    let b = Space::factor(Index(2));
    let ab = a.tensor_product(&b).unwrap();
    assert!(matches!(
        ab.tensor_product(&a),
        Err(SpaceError::DuplicateIndex)
    ));
    assert!(matches!(
        ab.tensor_product(&ab),
        Err(SpaceError::DuplicateIndex)
    ));
    assert!(matches!(
        Space::tensor([a.clone(), Space::factor(Index(1))]),
        Err(SpaceError::DuplicateIndex)
    ));
}

#[test]
fn empty_singleton_and_indexed_one_dimensional_factor_keep_their_meanings() {
    let empty = Space::tensor([]).unwrap();
    assert_eq!(empty.dimension(), Ok(1));
    assert_eq!(empty.indices().count(), 0);
    assert_eq!(empty.children().len(), 0);
    let a = Space::factor(Index(1));
    let singleton = Space::tensor([a.clone()]).unwrap();
    assert_eq!(a.children().len(), 0);
    assert_eq!(singleton.children().len(), 1);
    assert!(a.is_compatible_with(&singleton));
    assert!(a.is_compatible_with(&Space::tensor([empty, a.clone()]).unwrap()));
    let unit_factor = HilbertSpace::<_, Structure<1>>::factor(Index(1));
    assert_eq!(unit_factor.dimension(), Ok(1));
    assert!(!unit_factor.is_compatible_with(&HilbertSpace::unit()));
}

#[test]
fn large_dimension_is_an_error_only_when_queried() {
    let h = Space::tensor((0..usize::BITS + 1).map(|n| Space::factor(Index(n)))).unwrap();
    assert_eq!(h.dimension(), Err(DimensionOverflow));
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
                h = Space::tensor([h]).unwrap();
            }
            assert_eq!(h.dimension(), Ok(2));
            assert!(h.is_compatible_with(&leaf));
            assert_eq!(h.indices().next().unwrap().0, 8);
            drop(h);
            assert_eq!(leaf.factor_index().unwrap().0, 8);
        })
        .unwrap()
        .join()
        .unwrap();
}
