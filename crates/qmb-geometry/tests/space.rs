use qmb_core::HilbertSpace as H;
use qmb_geometry::{Boundary::*, Chain, ChainError};

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
        H::tensor([a.clone(), H::tensor([H::unit(), b.clone()]).unwrap()]).unwrap(),
        H::tensor([H::factor(8), H::factor(2), H::factor(7)]).unwrap(),
        H::unit(),
        H::factor(6),
    ];
    let chain = Chain::new(4, Periodic);
    assert_eq!(
        ids(chain.nn_pairs_with(&factors).unwrap()),
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
        ids(chain.nn_pairs_by_pos_with(&factors).unwrap()),
        [(9, 8), (3, 2), (9, 6)]
    );
    let (left, right) = chain.nn_pairs_with(&factors).unwrap().next().unwrap();
    assert!(std::ptr::eq(left, factors[0].children().next().unwrap()));
    assert!(left.shares_node_with(&a));
    assert!(std::ptr::eq(right, factors[1].children().next().unwrap()));
    let full = H::tensor(factors.iter().cloned()).unwrap();
    assert_eq!(
        ids(chain.nn_pairs_on(&full).unwrap()),
        ids(chain.nn_pairs_with(&factors).unwrap())
    );
    assert_eq!(
        ids(chain.nn_pairs_by_pos_on(&full).unwrap()),
        ids(chain.nn_pairs_by_pos_with(&factors).unwrap())
    );
    for ((a, b), (i, j)) in chain
        .nn_site_pairs_on(&full)
        .unwrap()
        .zip(chain.nn_site_pairs())
    {
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
    let single_empty = H::tensor([unit.clone()]).unwrap();
    let single_leaf = H::tensor([leaf.clone()]).unwrap();
    let nested = H::tensor([single_leaf.clone()]).unwrap();
    for (space, count) in [
        (&leaf, 1),
        (&unit, 0),
        (&single_empty, 1),
        (&single_leaf, 1),
        (&nested, 1),
    ] {
        let chain = Chain::new(count, Open);
        assert_eq!(chain.nn_site_pairs_on(space).unwrap().count(), 0);
        assert_eq!(chain.nn_pairs_on(space).unwrap().count(), 0);
        assert_eq!(chain.nn_pairs_by_pos_on(space).unwrap().count(), 0);
        let wrong = Chain::new(count + 1, Open);
        assert_eq!(
            wrong.nn_site_pairs_on(space).err(),
            Some(ChainError::SiteCountMismatch)
        );
        assert_eq!(
            wrong.nn_pairs_on(space).err(),
            Some(ChainError::SiteCountMismatch)
        );
        assert_eq!(
            wrong.nn_pairs_by_pos_on(space).err(),
            Some(ChainError::SiteCountMismatch)
        );
    }
    let full = H::tensor([
        H::tensor([H::factor(1), H::factor(2)]).unwrap(),
        H::factor(3),
    ])
    .unwrap();
    assert_eq!(
        Chain::new(3, Open).nn_pairs_on(&full).err(),
        Some(ChainError::SiteCountMismatch)
    );
    assert_eq!(
        ids(Chain::new(2, Open).nn_pairs_on(&full).unwrap()),
        [(1, 3), (2, 3)]
    );
}

#[test]
fn leaf_connections_need_neither_eq_nor_clone_and_validate_empty_cases() {
    struct Opaque;
    let factors = [H::factor(Opaque), H::factor(Opaque)];
    let chain = Chain::new(2, Open);
    let (a, b) = chain.nn_pairs_with(&factors).unwrap().next().unwrap();
    assert!(std::ptr::eq(a, &factors[0]));
    assert!(std::ptr::eq(b, &factors[1]));
    assert_eq!(chain.nn_pairs_by_pos_with(&factors).unwrap().count(), 1);
    for n in 0..=4 {
        let chain = Chain::new(n, Open);
        for count in 0..=5 {
            let sites: Vec<_> = (0..count).map(|_| H::<Opaque>::unit()).collect();
            let expected = (n != count).then_some(ChainError::SiteCountMismatch);
            assert_eq!(chain.nn_pairs_with(&sites).err(), expected);
            assert_eq!(chain.nn_pairs_by_pos_with(&sites).err(), expected);
        }
    }
}
