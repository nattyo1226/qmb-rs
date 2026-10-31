use qmb_geometry::{
    Boundary::*, Chain, Consecutive, ConsecutiveError, Geometry, MappingError, NearestNeighbors,
    Relation, SiteGroup, SiteGroupError, SiteMapper,
};

#[test]
fn boundaries_are_normalized_sorted_unique_and_loop_free() {
    for n in 0..7 {
        for boundary in [Open, Periodic] {
            let chain = Chain::new(n, boundary);
            let edges: Vec<_> = NearestNeighbors.groups(&chain).unwrap().collect();
            let mut expected: Vec<_> = (0..n.saturating_sub(1)).map(|i| [i, i + 1]).collect();
            if boundary == Periodic && n > 2 {
                expected.push([0, n - 1]);
                expected.sort();
            }
            assert_eq!(edges, expected);
            assert!(edges.iter().all(|&[i, j]| i < j && j < n));
            assert!(edges.windows(2).all(|w| w[0] < w[1]));
            assert_eq!(chain.num_sites(), n);
            assert_eq!(chain.boundary(), boundary);
            // No Clone, Eq, Default, Hash or Ord requirement on site values.
            struct Opaque(usize);
            let values: Vec<_> = (0..n).map(Opaque).collect();
            let mapper = SiteMapper::new(&chain, &values).unwrap();
            let mapped: Vec<_> = mapper
                .groups(&NearestNeighbors)
                .unwrap()
                .map(Result::unwrap)
                .collect();
            assert_eq!(mapped.len(), expected.len());
            for (group, [i, j]) in mapped.into_iter().zip(expected) {
                let (a, b) = (group.get(0).unwrap(), group.get(1).unwrap());
                assert_eq!((a.0, b.0), (i, j));
                assert!(std::ptr::eq(a, &values[i]));
                assert!(std::ptr::eq(b, &values[j]));
            }
        }
    }
}

#[test]
fn exact_length_is_validated_before_returning_even_without_pairs() {
    for n in 0..=4 {
        let chain = Chain::new(n, Periodic);
        for count in 0..=5 {
            let values = vec![(); count];
            assert_eq!(
                SiteMapper::new(&chain, &values).err(),
                (n != count).then_some(MappingError::SiteCountMismatch {
                    expected: n,
                    actual: count
                })
            );
        }
    }
    let error = MappingError::SiteCountMismatch {
        expected: 2,
        actual: 1,
    };
    assert_eq!(error.to_string(), "site count mismatch: expected 2, got 1");
    assert!(std::error::Error::source(&error).is_none());
}

// Small independently written expectations cover direction and periodic rotations.
fn consecutive<const K: usize>(n: usize, boundary: qmb_geometry::Boundary) -> Vec<[usize; K]> {
    Consecutive::<K>
        .groups(&Chain::new(n, boundary))
        .unwrap()
        .collect()
}

#[test]
fn consecutive_placements_preserve_roles_and_small_chain_boundaries() {
    for boundary in [Open, Periodic] {
        assert!(consecutive::<1>(0, boundary).is_empty());
        assert!(consecutive::<3>(2, boundary).is_empty());
        assert!(consecutive::<4>(3, boundary).is_empty());
        assert_eq!(consecutive::<1>(3, boundary), [[0], [1], [2]]);
        for n in [0, 1, 4] {
            assert_eq!(
                Consecutive::<0>.groups(&Chain::new(n, boundary)).err(),
                Some(ConsecutiveError::ZeroLength)
            );
        }
    }
    assert_eq!(consecutive::<3>(3, Open), [[0, 1, 2]]);
    assert_eq!(
        consecutive::<3>(3, Periodic),
        [[0, 1, 2], [1, 2, 0], [2, 0, 1]]
    );
    assert_eq!(consecutive::<3>(5, Open), [[0, 1, 2], [1, 2, 3], [2, 3, 4]]);
    assert_eq!(
        consecutive::<3>(5, Periodic),
        [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 0], [4, 0, 1]]
    );
    assert_eq!(consecutive::<4>(4, Open), [[0, 1, 2, 3]]);
    assert_eq!(
        consecutive::<4>(4, Periodic),
        [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]]
    );
    assert_eq!(consecutive::<4>(5, Open), [[0, 1, 2, 3], [1, 2, 3, 4]]);
    assert_eq!(
        consecutive::<4>(5, Periodic),
        [
            [0, 1, 2, 3],
            [1, 2, 3, 4],
            [2, 3, 4, 0],
            [3, 4, 0, 1],
            [4, 0, 1, 2]
        ]
    );
    assert_eq!(consecutive::<2>(2, Periodic), [[0, 1], [1, 0]]);
    // The first few groups at usize::MAX must not overflow or preallocate.
    let chain = Chain::new(usize::MAX, Periodic);
    assert_eq!(
        Consecutive::<3>
            .groups(&chain)
            .unwrap()
            .take(2)
            .collect::<Vec<_>>(),
        [[0, 1, 2], [1, 2, 3]]
    );
}

#[test]
fn bounded_groups_expose_only_valid_sites_without_imposing_set_semantics() {
    let group = SiteGroup::<5>::new([4, 1, 4]).unwrap();
    assert_eq!(group.as_ref(), [4, 1, 4]);
    assert_eq!(group.len(), 3);
    assert!(!group.is_empty());
    assert_eq!(SiteGroup::<0>::new([]).unwrap().as_ref(), []);
    assert_eq!(SiteGroup::<2>::default().as_ref(), []);
    assert_eq!(
        SiteGroup::<2>::new([0, 1, 2]),
        Err(SiteGroupError::CapacityExceeded {
            capacity: 2,
            len: 3
        })
    );
}

#[test]
fn custom_errors_and_borrowed_groups_are_preserved_with_lazy_validation() {
    use std::cell::Cell;
    struct UserError;
    struct Stored {
        fail: bool,
        groups: Vec<SiteGroup<3>>,
        consumed: Cell<usize>,
    }
    impl<G: Geometry> Relation<G> for Stored {
        type Error = UserError;
        type Group<'a>
            = &'a [usize]
        where
            Self: 'a,
            G: 'a;
        fn groups<'a>(
            &'a self,
            _: &'a G,
        ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
            if self.fail {
                return Err(UserError);
            }
            Ok(self.groups.iter().map(|group| {
                self.consumed.set(self.consumed.get() + 1);
                group.as_ref()
            }))
        }
    }
    struct Opaque;
    let chain = Chain::new(2, Open);
    let values = [Opaque, Opaque];
    let mapper = SiteMapper::new(&chain, &values).unwrap();
    assert!(std::ptr::eq(mapper.geometry(), &chain));
    assert!(std::ptr::eq(mapper.values(), &values));
    let mut rule = Stored {
        fail: true,
        groups: [vec![], vec![1], vec![1, 0, 1], vec![2], vec![0]]
            .map(|g| SiteGroup::new(g).unwrap())
            .into(),
        consumed: Cell::new(0),
    };
    assert!(matches!(mapper.groups(&rule), Err(UserError)));
    rule.fail = false;
    let mut groups = mapper
        .groups(&rule)
        .unwrap_or_else(|_| panic!("unexpected error"));
    assert_eq!(rule.consumed.get(), 0);
    assert!(groups.next().unwrap().unwrap().is_empty());
    let single = groups.next().unwrap().unwrap();
    assert!(std::ptr::eq(single.get(0).unwrap(), &values[1]));
    let repeated = groups.next().unwrap().unwrap();
    assert_eq!(repeated.len(), 3);
    assert!(repeated.get(3).is_none());
    assert!(std::ptr::eq(*repeated.source(), rule.groups[2].as_ref()));
    assert_eq!(
        groups.next().unwrap().err(),
        Some(MappingError::SiteOutOfRange {
            site: 2,
            num_sites: 2
        })
    );
    assert!(groups.next().unwrap().is_ok());
    assert!(groups.next().is_none());
    assert_eq!(rule.consumed.get(), 5);
}

#[test]
fn groups_can_borrow_stored_geometry_data() {
    struct Mesh {
        cells: Vec<[usize; 3]>,
    }
    impl Geometry for Mesh {
        fn num_sites(&self) -> usize {
            4
        }
    }
    struct Cells;
    impl Relation<Mesh> for Cells {
        type Error = core::convert::Infallible;
        type Group<'a> = &'a [usize];
        fn groups<'a>(
            &'a self,
            mesh: &'a Mesh,
        ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
            Ok(mesh.cells.iter().map(|cell| cell.as_slice()))
        }
    }
    let mesh = Mesh {
        cells: vec![[0, 1, 2], [2, 1, 3]],
    };
    let values = ['a', 'b', 'c', 'd'];
    let mapper = SiteMapper::new(&mesh, &values).unwrap();
    for (mapped, cell) in mapper.groups(&Cells).unwrap().zip(&mesh.cells) {
        let mapped = mapped.unwrap();
        assert!(std::ptr::eq(*mapped.source(), cell.as_slice()));
        assert_eq!(
            mapped.iter().copied().collect::<Vec<_>>(),
            cell.map(|i| values[i])
        );
    }
}
