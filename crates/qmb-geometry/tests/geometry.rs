use qmb_geometry::{
    Boundary::*, Chain, ChainError, NearestNeighborRelation, index_pairs, index_pairs_by_pos,
};
use std::sync::atomic::{AtomicUsize, Ordering};

// Neither Clone, Default, Hash nor Ord: existing index inputs move into a chain.
#[derive(Debug, PartialEq, Eq)]
struct Id(usize);
fn pairs(relation: &impl NearestNeighborRelation<Index = Id>) -> Vec<(usize, usize)> {
    relation.nn_pairs().map(|(a, b)| (a.0, b.0)).collect()
}

#[test]
fn boundaries_are_normalized_sorted_unique_and_loop_free() {
    for n in 0..7 {
        for boundary in [Open, Periodic] {
            let chain = Chain::from_sites((0..n).map(|i| [Id(i)]), boundary).unwrap();
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
        }
    }
}

#[test]
fn variable_length_sites_pair_orders_and_empty_sides() {
    let left = [Id(9), Id(3)];
    let right = [Id(8), Id(2), Id(7)];
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
        [Vec::from(left), Vec::from(right), vec![], vec![Id(6)]],
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
    assert_eq!(chain.site(0).unwrap(), [Id(9), Id(3)]);
    assert!(chain.site(2).unwrap().is_empty());
    assert_eq!(chain.site(99).err(), Some(ChainError::SiteOutOfBounds));
    assert!(chain.sites().get(99).is_none());
    assert_eq!(index_pairs(chain.site(0).unwrap(), &[]).count(), 0);
    assert_eq!(index_pairs(&[], chain.site(1).unwrap()).count(), 0);
}

#[test]
fn duplicates_within_and_across_sites_and_count_mismatches_are_rejected() {
    assert_eq!(
        Chain::from_sites([[Id(1), Id(1)]], Open).err(),
        Some(ChainError::DuplicateIndex)
    );
    assert_eq!(
        Chain::from_sites([[Id(1)], [Id(1)]], Open).err(),
        Some(ChainError::DuplicateIndex)
    );
    assert_eq!(
        Chain::with_sites(2, [[Id(1)]], Open).err(),
        Some(ChainError::SiteCountMismatch)
    );
    assert_eq!(
        Chain::with_sites(0, [[Id(1)]], Open).err(),
        Some(ChainError::SiteCountMismatch)
    );
    assert_eq!(
        Chain::with_sites(1, [[Id(1)]], Open).unwrap().sites().len(),
        1
    );
    // Equal identities may be reused in separately constructed chains.
    assert_eq!(
        Chain::from_sites([[Id(1)]], Open).unwrap().site(0).unwrap(),
        [Id(1)]
    );
}

static ISSUED: AtomicUsize = AtomicUsize::new(0);
#[derive(Debug, Clone, PartialEq, Eq)]
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
        chain
            .sites()
            .iter()
            .flatten()
            .map(|i| i.0)
            .collect::<Vec<_>>(),
        (start..start + 6).collect::<Vec<_>>()
    );
    let copy = chain.clone();
    assert_eq!(copy.sites(), chain.sites());
    assert_eq!(copy.boundary(), chain.boundary());
    assert_eq!(
        copy.site(0).unwrap(),
        [Generated(start), Generated(start + 1)]
    );
    assert_eq!(chain.nn_pairs().count(), 8);
    let empty = Chain::<Generated>::new(0, 3, Open).unwrap();
    assert!(empty.sites().is_empty());
    let units = Chain::<Generated>::new(4, 0, Periodic).unwrap();
    assert_eq!(units.sites().len(), 4);
    assert!(units.sites().iter().all(Vec::is_empty));
    assert_eq!(units.nn_site_pairs().count(), 4);
    assert_eq!(units.nn_pairs().count(), 0);
    assert_eq!(ISSUED.load(Ordering::SeqCst), start + 6);
    // Default alone does not guarantee freshness; validate ordinary equality.
    assert_eq!(
        Chain::<usize>::new(2, 1, Open).err(),
        Some(ChainError::DuplicateIndex)
    );
}

#[test]
fn borrowed_indices_need_no_static_lifetime_or_clone() {
    let names = [String::from("left"), String::from("right")];
    let chain = Chain::from_sites(names.iter().map(|s| [s.as_str()]), Open).unwrap();
    assert_eq!(
        chain.nn_pairs().map(|(&a, &b)| (a, b)).collect::<Vec<_>>(),
        [("left", "right")]
    );
}

// Relationship implementations need no common Geometry container.
struct TwoSites([Vec<Id>; 2]);
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
    let other = TwoSites([vec![Id(4)], vec![Id(1), Id(9)]]);
    let chain = Chain::from_sites([vec![Id(4)], vec![Id(1), Id(9)]], Open).unwrap();
    assert_eq!(pairs(&other), pairs(&chain));
    assert_eq!(other.nn_pairs_by_pos().count(), 1);
}

#[test]
fn chain_errors_have_messages_and_no_sources() {
    for (error, message) in [
        (ChainError::SiteOutOfBounds, "site number is out of bounds"),
        (ChainError::SiteCountMismatch, "site count mismatch"),
        (ChainError::DuplicateIndex, "duplicate factor index"),
    ] {
        assert_eq!(error.to_string(), message);
        assert!(std::error::Error::source(&error).is_none());
    }
}
