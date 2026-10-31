use num_complex::Complex64 as C;
use qmb_core::HilbertSpace as H;
use qmb_fermion::{Fermion, FermionOperator as F, FermionState as Occupation};
use qmb_geometry::{Boundary, Chain, Consecutive, NearestNeighbors, SiteMapper};
use qmb_spin_half::{SpinHalf, SpinHalfOperator as S, SpinHalfState as Spin};

fn close(actual: C, expected: C) {
    assert!(
        (actual - expected).norm_sqr() < 1e-24,
        "{actual} != {expected}"
    );
}

// Independent tensor product of literal Pauli matrices, I=0, X=1, Y=2, Z=3.
fn pauli_element(labels: &[usize], row: usize, col: usize) -> C {
    let mut value = C::new(1., 0.);
    for (site, &axis) in labels.iter().enumerate() {
        let r = (row >> (labels.len() - 1 - site)) & 1;
        let c = (col >> (labels.len() - 1 - site)) & 1;
        value *= match (axis, r, c) {
            (0, a, b) if a == b => C::new(1., 0.),
            (1, a, b) if a != b => C::new(1., 0.),
            (2, 0, 1) => C::new(0., -1.),
            (2, 1, 0) => C::new(0., 1.),
            (3, 0, 0) => C::new(1., 0.),
            (3, 1, 1) => C::new(-1., 0.),
            _ => C::new(0., 0.),
        };
    }
    value
}

fn check_spin_placements<const K: usize>(boundary: Boundary, expected: &[[usize; K]]) {
    let sites: Vec<_> = (0..5).map(|i| H::factor(10 - i)).collect();
    let ambient = H::tensor_all(sites.iter().cloned()).unwrap();
    let chain = Chain::new(5, boundary);
    let mapper = SiteMapper::new(&chain, &sites).unwrap();
    let rule = Consecutive::<K>;
    let groups: Vec<_> = mapper.groups(&rule).unwrap().map(Result::unwrap).collect();
    assert_eq!(groups.len(), expected.len());
    let states: Vec<_> = (0..32)
        .map(|bits| {
            SpinHalf::<C>::product_state(
                &ambient,
                (0..5).map(|i| {
                    if bits & (1 << (4 - i)) == 0 {
                        Spin::Up
                    } else {
                        Spin::Down
                    }
                }),
            )
            .unwrap()
        })
        .collect();
    for (group, placement) in groups.iter().zip(expected) {
        let axes = [S::Sx, S::Sy, S::Sz, S::Sx];
        let leaves = group.leaves_cartesian().next().unwrap();
        let mapped = SpinHalf::<C>::try_operator_product(
            leaves
                .iter()
                .zip(axes)
                .map(|(leaf, axis)| SpinHalf::<C>::operator(axis).scaled(2.).at(leaf)),
        )
        .unwrap()
        .bind(&ambient)
        .unwrap();
        let manual = SpinHalf::<C>::try_operator_product(
            placement
                .iter()
                .zip(axes)
                .map(|(&i, axis)| SpinHalf::<C>::operator(axis).scaled(2.).at(&sites[i])),
        )
        .unwrap()
        .bind(&ambient)
        .unwrap();
        let mut labels = [0; 5];
        for (slot, &site) in placement.iter().enumerate() {
            labels[site] = [1, 2, 3, 1][slot];
        }
        for (row, bra) in states.iter().enumerate() {
            for (col, ket) in states.iter().enumerate() {
                let reference = pauli_element(&labels, row, col);
                close(
                    bra.inner_product(&ket.act(&mapped).unwrap()).unwrap(),
                    reference,
                );
                close(
                    bra.inner_product(&ket.act(&manual).unwrap()).unwrap(),
                    reference,
                );
            }
        }
    }
}

#[test]
fn three_and_four_site_ordered_pauli_terms_match_manual_and_independent_matrices() {
    check_spin_placements::<3>(Boundary::Open, &[[0, 1, 2], [1, 2, 3], [2, 3, 4]]);
    check_spin_placements::<3>(
        Boundary::Periodic,
        &[[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 0], [4, 0, 1]],
    );
    check_spin_placements::<4>(Boundary::Open, &[[0, 1, 2, 3], [1, 2, 3, 4]]);
    check_spin_placements::<4>(
        Boundary::Periodic,
        &[
            [0, 1, 2, 3],
            [1, 2, 3, 4],
            [2, 3, 4, 0],
            [3, 4, 0, 1],
            [4, 0, 1, 2],
        ],
    );
}

#[test]
fn mapper_constructs_all_37_pauli_strings_with_support_at_most_two_once() {
    let sites: Vec<_> = (0..3).map(H::factor).collect();
    let ambient = H::tensor_all(sites.iter().cloned()).unwrap();
    // Three periodic sites have all three undirected pairs. Singles use another
    // geometric rule on the same mapper; identity is constructed by core.
    let chain = Chain::new(3, Boundary::Periodic);
    let mapper = SiteMapper::new(&chain, &sites).unwrap();
    let mut supports = vec![vec![]];
    for group in mapper.groups(&Consecutive::<1>).unwrap() {
        supports.push(group.unwrap().iter().collect());
    }
    for group in mapper.groups(&NearestNeighbors).unwrap() {
        supports.push(group.unwrap().iter().collect());
    }
    let states: Vec<_> = (0..8)
        .map(|bits| {
            SpinHalf::<C>::product_state(
                &ambient,
                (0..3).map(|i| {
                    if bits & (1 << (2 - i)) == 0 {
                        Spin::Up
                    } else {
                        Spin::Down
                    }
                }),
            )
            .unwrap()
        })
        .collect();
    let mut seen = std::collections::BTreeSet::new();
    for support in supports {
        for word in 0..3usize.pow(support.len() as u32) {
            let mut digits = word;
            let mut labels = [0; 3];
            let mut factors = Vec::new();
            for site in &support {
                let axis = digits % 3;
                digits /= 3;
                labels[*site.factor_index().unwrap()] = axis + 1;
                factors.push(
                    SpinHalf::<C>::operator([S::Sx, S::Sy, S::Sz][axis])
                        .scaled(2.)
                        .at(site),
                );
            }
            assert!(seen.insert(labels));
            let operator = SpinHalf::<C>::try_operator_product(factors)
                .unwrap()
                .bind(&ambient)
                .unwrap();
            for (row, bra) in states.iter().enumerate() {
                for (col, ket) in states.iter().enumerate() {
                    close(
                        bra.inner_product(&ket.act(&operator).unwrap()).unwrap(),
                        pauli_element(&labels, row, col),
                    );
                }
            }
        }
    }
    let expected: std::collections::BTreeSet<_> = (0..4)
        .flat_map(|a| (0..4).flat_map(move |b| (0..4).map(move |c| [a, b, c])))
        .filter(|labels| labels.iter().filter(|&&axis| axis != 0).count() <= 2)
        .collect();
    assert_eq!(seen, expected);
    assert_eq!(seen.len(), 37);
}

#[test]
fn graded_three_and_four_site_products_preserve_periodic_role_and_odd_order() {
    let sites: Vec<_> = (0..4).map(|i| H::factor(20 - i)).collect();
    let ambient = H::tensor_all_g(sites.iter().cloned()).unwrap();
    let states: Vec<_> = (0..16)
        .map(|bits| {
            Fermion::<C>::product_state(
                &ambient,
                (0..4).map(|i| {
                    if bits & (1 << (3 - i)) == 0 {
                        Occupation::Empty
                    } else {
                        Occupation::Occupied
                    }
                }),
            )
            .unwrap()
        })
        .collect();
    for boundary in [Boundary::Open, Boundary::Periodic] {
        let chain = Chain::new(4, boundary);
        let mapper = SiteMapper::new(&chain, &sites).unwrap();
        for k in [3, 4] {
            let supports: Vec<Vec<_>> = if k == 3 {
                mapper
                    .groups(&Consecutive::<3>)
                    .unwrap()
                    .map(|group| group.unwrap().iter().collect())
                    .collect()
            } else {
                mapper
                    .groups(&Consecutive::<4>)
                    .unwrap()
                    .map(|group| group.unwrap().iter().collect())
                    .collect()
            };
            let expected = if boundary == Boundary::Open { 5 - k } else { 4 };
            assert_eq!(supports.len(), expected);
            for (start, support) in supports.iter().enumerate() {
                let placement: Vec<_> = (0..k).map(|offset| (start + offset) % 4).collect();
                let mapped = Fermion::<C>::operator_product_at(
                    support.iter().map(|&leaf| (leaf, F::Create)),
                )
                .unwrap()
                .bind(&ambient)
                .unwrap();
                let manual = Fermion::<C>::operator_product_at(
                    placement.iter().map(|&i| (&sites[i], F::Create)),
                )
                .unwrap()
                .bind(&ambient)
                .unwrap();
                // Independent Fock action: rightmost operator acts first. Earlier
                // occupied tensor positions supply each Jordan-Wigner sign.
                for (col, ket) in states.iter().enumerate() {
                    let mut result = Some((col, 1.));
                    for &mode in placement.iter().rev() {
                        result = result.and_then(|(bits, sign)| {
                            let mask = 1 << (3 - mode);
                            if bits & mask != 0 {
                                return None;
                            }
                            let earlier = (bits >> (4 - mode)).count_ones();
                            Some((bits | mask, if earlier % 2 == 0 { sign } else { -sign }))
                        });
                    }
                    for (row, bra) in states.iter().enumerate() {
                        let reference =
                            C::new(result.filter(|&(r, _)| r == row).map_or(0., |(_, s)| s), 0.);
                        close(
                            bra.inner_product(&ket.act(&mapped).unwrap()).unwrap(),
                            reference,
                        );
                        close(
                            bra.inner_product(&ket.act(&manual).unwrap()).unwrap(),
                            reference,
                        );
                    }
                }
            }
        }
    }
}
