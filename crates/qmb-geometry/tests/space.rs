use qmb_core::HilbertSpace;
use qmb_core::HilbertSpace as H;
use qmb_geometry::{
    Boundary::*, Chain, Consecutive, MappingError, NearestNeighbors, Relation, SiteMapper,
};

fn pairs<'s, I>(
    mapper: &SiteMapper<'_, 's, Chain, H<I>>,
    zip: bool,
) -> std::vec::IntoIter<(&'s H<I>, &'s H<I>)> {
    mapper
        .groups(&NearestNeighbors)
        .unwrap()
        .flat_map(|group| {
            let group = group.unwrap();
            let rows: Vec<_> = if zip {
                group.leaves_zip().collect()
            } else {
                group.leaves_cartesian().collect()
            };
            rows.into_iter().map(|row| (row[0], row[1]))
        })
        .collect::<Vec<_>>()
        .into_iter()
}

fn ids<'a>(pairs: impl Iterator<Item = (&'a H<u8>, &'a H<u8>)>) -> Vec<(u8, u8)> {
    pairs
        .map(|(a, b)| (*a.factor_index().unwrap(), *b.factor_index().unwrap()))
        .collect()
}

#[test]
fn nested_variable_length_sites_cartesian_and_positional_pairs_keep_order() {
    let a = H::factor(9);
    let b = H::factor(3);
    let factors = [
        H::tensor_all([a.clone(), H::tensor_all([H::unit(), b.clone()]).unwrap()]).unwrap(),
        H::tensor_all([H::factor(8), H::factor(2), H::factor(7)]).unwrap(),
        H::unit(),
        H::factor(6),
    ];
    let chain = Chain::new(4, Periodic);
    assert_eq!(
        ids(pairs(&SiteMapper::new(&chain, &factors).unwrap(), false)),
        [
            (9, 8),
            (9, 2),
            (9, 7),
            (3, 8),
            (3, 2),
            (3, 7),
            (9, 6),
            (3, 6)
        ]
    );
    assert_eq!(
        ids(pairs(&SiteMapper::new(&chain, &factors).unwrap(), true)),
        [(9, 8), (3, 2), (9, 6)]
    );
    let (left, right) = pairs(&SiteMapper::new(&chain, &factors).unwrap(), false)
        .next()
        .unwrap();
    assert!(std::ptr::eq(left, factors[0].children().next().unwrap()));
    assert!(left.shares_node_with(&a));
    assert!(std::ptr::eq(right, factors[1].children().next().unwrap()));
    let full = H::tensor_all(factors.iter().cloned()).unwrap();
    assert_eq!(
        ids(pairs(
            &SiteMapper::from_space(&chain, &full).unwrap(),
            false
        )),
        ids(pairs(&SiteMapper::new(&chain, &factors).unwrap(), false))
    );
    assert_eq!(
        ids(pairs(&SiteMapper::from_space(&chain, &full).unwrap(), true)),
        ids(pairs(&SiteMapper::new(&chain, &factors).unwrap(), true))
    );
    let mapper = SiteMapper::from_space(&chain, &full).unwrap();
    for (group, [i, j]) in mapper
        .groups(&NearestNeighbors)
        .unwrap()
        .zip(NearestNeighbors.groups(&chain).unwrap())
    {
        let group = group.unwrap();
        let (a, b) = (group.get(0).unwrap(), group.get(1).unwrap());
        assert!(std::ptr::eq(a, &full.children().as_slice()[i]));
        assert!(std::ptr::eq(b, &full.children().as_slice()[j]));
        assert!(a.shares_node_with(&factors[i]));
    }
    assert_eq!(full.children().count(), 4);
    assert_eq!(full.children().next().unwrap().children().count(), 2);
}

#[test]
fn on_reads_top_level_partition_and_distinguishes_empty_and_singleton_tensors() {
    let leaf = H::<u8>::factor(0);
    let unit = H::unit();
    let single_empty = H::tensor_all([unit.clone()]).unwrap();
    let single_leaf = H::tensor_all([leaf.clone()]).unwrap();
    let nested = H::tensor_all([single_leaf.clone()]).unwrap();
    for (space, count) in [
        (&leaf, 1),
        (&unit, 0),
        (&single_empty, 1),
        (&single_leaf, 1),
        (&nested, 1),
    ] {
        let chain = Chain::new(count, Open);
        assert_eq!(
            SiteMapper::from_space(&chain, space)
                .unwrap()
                .groups(&NearestNeighbors)
                .unwrap()
                .count(),
            0
        );
        assert_eq!(
            pairs(&SiteMapper::from_space(&chain, space).unwrap(), false).count(),
            0
        );
        assert_eq!(
            pairs(&SiteMapper::from_space(&chain, space).unwrap(), true).count(),
            0
        );
        let wrong = Chain::new(count + 1, Open);
        assert_eq!(
            SiteMapper::from_space(&wrong, space).err(),
            Some(MappingError::SiteCountMismatch {
                expected: count + 1,
                actual: count
            })
        );
    }
    let full = H::tensor_all([
        H::tensor_all([H::factor(1), H::factor(2)]).unwrap(),
        H::factor(3),
    ])
    .unwrap();
    assert_eq!(
        SiteMapper::from_space(&Chain::new(3, Open), &full).err(),
        Some(MappingError::SiteCountMismatch {
            expected: 3,
            actual: 2
        })
    );
    assert_eq!(
        ids(pairs(
            &SiteMapper::from_space(&Chain::new(2, Open), &full).unwrap(),
            false
        )),
        [(1, 3), (2, 3)]
    );
}

#[test]
fn leaf_connections_need_neither_eq_nor_clone_and_validate_empty_cases() {
    struct Opaque;
    let factors = [H::factor(Opaque), H::factor(Opaque)];
    let chain = Chain::new(2, Open);
    let (a, b) = pairs(&SiteMapper::new(&chain, &factors).unwrap(), false)
        .next()
        .unwrap();
    assert!(std::ptr::eq(a, &factors[0]));
    assert!(std::ptr::eq(b, &factors[1]));
    assert_eq!(
        pairs(&SiteMapper::new(&chain, &factors).unwrap(), true).count(),
        1
    );
    for n in 0..=4 {
        let chain = Chain::new(n, Open);
        for count in 0..=5 {
            let sites: Vec<_> = (0..count).map(|_| H::<Opaque>::unit()).collect();
            let expected = (n != count).then_some(MappingError::SiteCountMismatch {
                expected: n,
                actual: count,
            });
            assert_eq!(SiteMapper::new(&chain, &sites).err(), expected);
        }
    }
}

// Arbitrary saved groups are test fixtures; subset enumeration is not a Relation API.
struct TestGroups {
    groups: Vec<Vec<usize>>,
}
impl<G: qmb_geometry::Geometry> Relation<G> for TestGroups {
    type Error = core::convert::Infallible;
    type Group<'a>
        = &'a [usize]
    where
        Self: 'a,
        G: 'a;
    fn groups<'a>(
        &'a self,
        _: &'a G,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
        Ok(self.groups.iter().map(Vec::as_slice))
    }
}
fn leaf_ids<'a>(rows: impl Iterator<Item = Vec<&'a HilbertSpace<u8>>>) -> Vec<Vec<u8>> {
    rows.map(|row| {
        row.into_iter()
            .map(|leaf| *leaf.factor_index().unwrap())
            .collect()
    })
    .collect()
}

#[test]
fn leaf_combiners_keep_tensor_order_for_uneven_nested_sites() {
    let sites = [
        HilbertSpace::tensor_all([
            HilbertSpace::factor(9u8),
            HilbertSpace::tensor_all([HilbertSpace::unit(), HilbertSpace::factor(3)]).unwrap(),
        ])
        .unwrap(),
        HilbertSpace::tensor_all([
            HilbertSpace::factor(8),
            HilbertSpace::factor(2),
            HilbertSpace::factor(7),
        ])
        .unwrap(),
        HilbertSpace::tensor_all([HilbertSpace::factor(6), HilbertSpace::factor(5)]).unwrap(),
    ];
    let chain = Chain::new(3, Open);
    let mapped = SiteMapper::new(&chain, &sites).unwrap();
    let rule = Consecutive::<3>;
    let group = mapped.groups(&rule).unwrap().next().unwrap().unwrap();
    assert_eq!(
        leaf_ids(group.leaves_cartesian()),
        [
            [9, 8, 6],
            [9, 8, 5],
            [9, 2, 6],
            [9, 2, 5],
            [9, 7, 6],
            [9, 7, 5],
            [3, 8, 6],
            [3, 8, 5],
            [3, 2, 6],
            [3, 2, 5],
            [3, 7, 6],
            [3, 7, 5],
        ]
    );
    assert_eq!(leaf_ids(group.leaves_zip()), [[9, 8, 6], [3, 2, 5]]);
    let custom = group.combine_leaves(|mut streams| {
        streams.reverse();
        let mut done = false;
        std::iter::from_fn(move || {
            if done {
                return None;
            }
            let row: Option<Vec<_>> = streams.iter_mut().map(Iterator::next).collect();
            if row.is_none() {
                done = true;
            }
            row
        })
    });
    assert_eq!(leaf_ids(custom), [[6, 8, 9], [5, 2, 3]]);
    let expected: Vec<_> = sites
        .iter()
        .map(|site| site.leaves().collect::<Vec<_>>())
        .collect();
    for row in group.leaves_cartesian() {
        for (slot, leaf) in row.into_iter().enumerate() {
            assert!(
                expected[slot]
                    .iter()
                    .any(|original| std::ptr::eq(*original, leaf))
            );
        }
    }
    assert_eq!(*group.source(), [0, 1, 2]);
}

#[test]
fn leaf_combiners_handle_no_sites_empty_site_single_site_and_exhaustion() {
    let chain = Chain::new(1, Open);
    let sites = [HilbertSpace::factor(4u8)];
    let mapped = SiteMapper::new(&chain, &sites).unwrap();
    let rule = TestGroups {
        groups: vec![vec![], vec![0]],
    };
    let groups: Vec<_> = mapped.groups(&rule).unwrap().map(Result::unwrap).collect();
    assert_eq!(leaf_ids(groups[0].leaves_cartesian()), vec![vec![]]);
    assert_eq!(leaf_ids(groups[0].leaves_zip()), vec![vec![]]);
    assert_eq!(leaf_ids(groups[1].leaves_cartesian()), vec![vec![4]]);
    assert_eq!(leaf_ids(groups[1].leaves_zip()), vec![vec![4]]);
    for group in &groups {
        let mut product = group.leaves_cartesian();
        assert!(product.next().is_some());
        assert!(product.next().is_none());
        assert!(product.next().is_none());
        let mut zipped = group.leaves_zip();
        assert!(zipped.next().is_some());
        assert!(zipped.next().is_none());
        assert!(zipped.next().is_none());
    }
    let sites = [HilbertSpace::<u8>::unit(), HilbertSpace::factor(7)];
    let chain = Chain::new(2, Open);
    let mapped = SiteMapper::new(&chain, &sites).unwrap();
    let rule = NearestNeighbors;
    let group = mapped.groups(&rule).unwrap().next().unwrap().unwrap();
    let mut product = group.leaves_cartesian();
    assert!(product.next().is_none());
    assert!(product.next().is_none());
    let mut zipped = group.leaves_zip();
    assert!(zipped.next().is_none());
    assert!(zipped.next().is_none());
}

#[test]
fn leaf_combination_needs_no_index_traits_and_accepts_fn_once() {
    struct Opaque;
    let sites = [HilbertSpace::factor(Opaque), HilbertSpace::factor(Opaque)];
    let chain = Chain::new(2, Open);
    let mapped = SiteMapper::new(&chain, &sites).unwrap();
    let rule = NearestNeighbors;
    let group = mapped.groups(&rule).unwrap().next().unwrap().unwrap();
    assert_eq!(group.leaves_cartesian().count(), 1);
    assert_eq!(group.leaves_zip().count(), 1);
    let token = String::from("consumed by the combiner");
    let mut custom = group.combine_leaves(move |mut streams| {
        drop(token);
        streams.remove(1).map(|leaf| vec![leaf])
    });
    let row = custom.next().unwrap();
    assert!(std::ptr::eq(row[0], &sites[1]));
    assert!(custom.next().is_none());
}
