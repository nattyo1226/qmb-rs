use qmb_core::{HilbertSpace, SpaceError};
use qmb_geometry::{
    Boundary::*, Chain, Geometry, GeometryError, NearestNeighborRelation, index_pairs,
    index_pairs_by_pos,
};
use std::{
    num::NonZeroUsize,
    sync::atomic::{AtomicUsize, Ordering},
};

// Neither Clone nor Default; holding and enumerating existing spaces need neither.
#[derive(Debug, PartialEq, Eq)]
struct Id(usize);
type H = HilbertSpace<Id>;
fn factor(n: usize) -> H {
    H::factor(Id(n))
}
fn pairs(relation: &impl NearestNeighborRelation<Index = Id>) -> Vec<(usize, usize)> {
    relation.nn_pairs().map(|(a, b)| (a.0, b.0)).collect()
}

#[test]
fn boundaries_are_normalized_sorted_unique_and_loop_free() {
    for n in 0..7 {
        for boundary in [Open, Periodic] {
            let chain = Chain::from_sites((0..n).map(factor), boundary).unwrap();
            let edges: Vec<_> = chain.nn_site_pairs().collect();
            let mut expected: Vec<_> = (0..n.saturating_sub(1)).map(|i| (i, i + 1)).collect();
            if boundary == Periodic && n > 2 {
                expected.push((0, n - 1));
                expected.sort();
            }
            assert_eq!(edges, expected);
            assert!(edges.iter().all(|&(i, j)| i < j && j < n));
            assert!(edges.windows(2).all(|w| w[0] < w[1]));
            assert_eq!(pairs(&chain), expected);
            assert_eq!(chain.boundary(), boundary);
            assert!(chain.space().shares_node_with(chain.geometry().space()));
        }
    }
}

#[test]
fn nested_sites_pair_orders_and_empty_sides() {
    let a = factor(9);
    let b = factor(3);
    let left = H::tensor([H::unit(), H::tensor([a.clone(), b.clone()]).unwrap()]).unwrap();
    let right = H::tensor([factor(8), H::tensor([factor(2), factor(7)]).unwrap()]).unwrap();
    let all = vec![(9, 8), (9, 2), (9, 7), (3, 8), (3, 2), (3, 7)];
    assert_eq!(
        index_pairs(&left, &right)
            .map(|(a, b)| (a.0, b.0))
            .collect::<Vec<_>>(),
        all
    );
    assert_eq!(
        index_pairs_by_pos(&left, &right)
            .map(|(a, b)| (a.0, b.0))
            .collect::<Vec<_>>(),
        [(9, 8), (3, 2)]
    );
    let chain = Chain::from_sites(
        [left.clone(), right.clone(), H::unit(), factor(6)],
        Periodic,
    )
    .unwrap();
    assert_eq!(pairs(&chain), [all, vec![(9, 6), (3, 6)]].concat());
    assert_eq!(
        chain
            .nn_pairs_by_pos()
            .map(|(a, b)| (a.0, b.0))
            .collect::<Vec<_>>(),
        [(9, 8), (3, 2), (9, 6)]
    );
    let copy = chain.clone();
    assert!(copy.space().shares_node_with(chain.space()));
    assert!(copy.sites()[0].shares_node_with(&left));
    assert!(
        copy.sites()[0]
            .children()
            .nth(1)
            .unwrap()
            .children()
            .next()
            .unwrap()
            .shares_node_with(&a)
    );
    assert_eq!(
        chain.geometry().site(99).err(),
        Some(GeometryError::SiteOutOfBounds)
    );
    assert!(chain.sites().get(99).is_none());
    assert_eq!(index_pairs(&left, &H::unit()).count(), 0);
    assert_eq!(index_pairs(&H::unit(), &right).count(), 0);
}

#[test]
fn existing_full_order_grouping_and_errors() {
    let a = factor(5);
    let b = factor(1);
    let c = factor(8);
    let full = H::tensor([c.clone(), H::tensor([a.clone(), b.clone()]).unwrap()]).unwrap();
    let geometry = Geometry::with_space([a.clone(), b.clone(), c.clone()], &full).unwrap();
    assert!(geometry.space().shares_node_with(&full));
    assert_eq!(
        geometry.space().indices().map(|i| i.0).collect::<Vec<_>>(),
        [8, 5, 1]
    );
    assert!(geometry.sites()[0].shares_node_with(&a));
    assert_eq!(
        Geometry::from_sites([a.clone(), a.clone()]).err(),
        Some(GeometryError::Space(SpaceError::DuplicateIndex))
    );
    assert_eq!(
        Geometry::with_space([a.clone(), b.clone()], &full).err(),
        Some(GeometryError::IncompatibleSpace)
    );
    assert_eq!(
        Geometry::with_space([a.clone(), factor(99), c], &full).err(),
        Some(GeometryError::IncompatibleSpace)
    );
    assert_eq!(
        Chain::with_sites(2, [a], Open).err(),
        Some(GeometryError::SiteCountMismatch)
    );
}

static ISSUED: AtomicUsize = AtomicUsize::new(0);
#[derive(Debug, PartialEq, Eq)]
struct Generated(usize);
impl Default for Generated {
    fn default() -> Self {
        Self(ISSUED.fetch_add(1, Ordering::SeqCst))
    }
}
#[test]
fn generation_is_exact_ordered_and_never_repeated_by_access_or_clone() {
    let start = ISSUED.load(Ordering::SeqCst);
    let chain = Chain::<Generated>::new(3, 2, Open).unwrap();
    assert_eq!(ISSUED.load(Ordering::SeqCst), start + 6);
    assert_eq!(
        chain.space().indices().map(|i| i.0).collect::<Vec<_>>(),
        (start..start + 6).collect::<Vec<_>>()
    );
    let copy = chain.clone();
    assert!(copy.space().shares_node_with(chain.space()));
    assert_eq!(chain.nn_pairs().count(), 8);
    let empty = Chain::<Generated>::new(0, 3, Open).unwrap();
    assert!(empty.sites().is_empty());
    assert_eq!(
        empty.space().dimension_with(NonZeroUsize::new(2).unwrap()),
        Ok(1)
    );
    let units = Chain::<Generated>::new(4, 0, Periodic).unwrap();
    assert_eq!(units.sites().len(), 4);
    assert_eq!(
        units.space().dimension_with(NonZeroUsize::new(2).unwrap()),
        Ok(1)
    );
    assert!(units.sites().iter().all(|s| s.indices().next().is_none()));
    assert_eq!(units.nn_pairs().count(), 0);
    assert_eq!(ISSUED.load(Ordering::SeqCst), start + 6);
    // A Default implementation that reuses identity violates the issuance
    // contract and is caught by ordinary space validation.
    assert_eq!(
        Chain::<usize>::new(2, 1, Open).err(),
        Some(GeometryError::Space(SpaceError::DuplicateIndex))
    );
}

// A second, private geometry implements the same public capability. Shared
// conversion functions require neither custom index generation nor Clone.
struct TwoSites([H; 2]);
impl NearestNeighborRelation for TwoSites {
    type Index = Id;
    fn nn_pairs(&self) -> impl Iterator<Item = (&Id, &Id)> + '_ {
        index_pairs(&self.0[0], &self.0[1])
    }
    fn nn_pairs_by_pos(&self) -> impl Iterator<Item = (&Id, &Id)> + '_ {
        index_pairs_by_pos(&self.0[0], &self.0[1])
    }
}
#[test]
fn relationship_capability_is_shared_by_distinct_geometries() {
    let other = TwoSites([factor(4), H::tensor([factor(1), factor(9)]).unwrap()]);
    let chain = Chain::from_sites(other.0.iter().cloned(), Open).unwrap();
    assert_eq!(pairs(&other), pairs(&chain));
    assert_eq!(other.nn_pairs_by_pos().count(), 1);
}

#[test]
fn geometry_error_conversion_preserves_message_and_source() {
    let original = SpaceError::DuplicateIndex;
    let error = GeometryError::from(original);
    assert_eq!(error, GeometryError::Space(original));
    assert_eq!(error.to_string(), original.to_string());
    assert_eq!(
        std::error::Error::source(&error).unwrap().to_string(),
        original.to_string()
    );
    for (error, message) in [
        (
            GeometryError::SiteOutOfBounds,
            "site number is out of bounds",
        ),
        (GeometryError::SiteCountMismatch, "site count mismatch"),
        (
            GeometryError::IncompatibleSpace,
            "sites do not partition the full space",
        ),
    ] {
        assert_eq!(error.to_string(), message);
        assert!(std::error::Error::source(&error).is_none());
    }
}
