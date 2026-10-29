mod support;
use num_complex::Complex;
use qmb_core::{HilbertSpace as H, QuantumSystem, SpaceError, TensorKind};
use qmb_fermion::{
    Fermion, FermionOperator as O, FermionState as S, LadderBasis, MajoranaBasis, OccupationBasis,
    UtilityError,
};
use qmb_geometry::{Boundary, Chain, ChainError};
use qmb_normalize::{Norm, Normalize, NormalizedExpectation};
use support::*;

fn sites(n: usize) -> Vec<H<usize>> {
    // Deliberately descending indices: ordering comes from the tensor.
    (0..n)
        .map(|j| H::tensor_g([H::factor(2 * (n - j)), H::factor(2 * (n - j) - 1)]).unwrap())
        .collect()
}
fn reference(n: usize, boundary: Boundary, t: f64, u: f64, symmetric: bool) -> Matrix {
    let dim: usize = 1 << (2 * n);
    let mut m = vec![Complex::new(0., 0.); dim * dim];
    let occupation = |bits: usize, mode: usize| ((bits >> (2 * n - 1 - mode)) & 1) as f64;
    // Independent edges; do not use Chain or production geometry.
    let mut edges: Vec<_> = (1..n).map(|j| (j - 1, j)).collect();
    if boundary == Boundary::Periodic && n > 2 {
        edges.push((0, n - 1));
    }
    for col in 0..dim {
        for j in 0..n {
            let (up, down) = (occupation(col, 2 * j), occupation(col, 2 * j + 1));
            let interaction = if symmetric {
                (up - 0.5) * (down - 0.5)
            } else {
                up * down
            };
            m[col * dim + col] += Complex::new(u * interaction, 0.);
        }
        for &(left, right) in &edges {
            for spin in 0..2 {
                for (from, to) in [
                    (2 * left + spin, 2 * right + spin),
                    (2 * right + spin, 2 * left + spin),
                ] {
                    if let Some((middle, s1)) = fock_action(col, from, 2 * n, false) {
                        if let Some((row, s2)) = fock_action(middle, to, 2 * n, true) {
                            m[row * dim + col] += Complex::new(-t * s1 * s2, 0.);
                        }
                    }
                }
            }
        }
    }
    m
}
fn check<T: Real>() {
    for n in 0..=3 {
        let sites = sites(n);
        for boundary in [Boundary::Open, Boundary::Periodic] {
            let chain = Chain::new(n, boundary);
            for (t, u) in [(0.75, 1.25), (0., 2.), (1., 0.)] {
                let standard = Fermion::<Complex<T>>::hubbard(
                    &chain,
                    &sites,
                    T::from_f64(t).unwrap(),
                    T::from_f64(u).unwrap(),
                )
                .unwrap();
                let symmetric = Fermion::<Complex<T>>::symmetric_hubbard(
                    &chain,
                    &sites,
                    T::from_f64(t).unwrap(),
                    T::from_f64(u).unwrap(),
                )
                .unwrap();
                let std = exported_matrix(&standard);
                let sym = exported_matrix(&symmetric);
                values::<T>(&std, &reference(n, boundary, t, u, false));
                values::<T>(&sym, &reference(n, boundary, t, u, true));
                values::<T>(&std, &adjoint(&std));
                values::<T>(&sym, &adjoint(&sym));
                close::<T>(to64(standard.trace().unwrap()), trace(&std));
                close::<T>(to64(symmetric.trace().unwrap()), trace(&sym));
                assert_eq!(standard.space().unwrap().kind(), Some(TensorKind::Graded));
                for (given, retained) in sites.iter().zip(standard.space().unwrap().children()) {
                    assert!(given.shares_node_with(retained));
                }
                let dim: usize = 1 << (2 * n);
                for row in 0..dim {
                    for col in 0..dim {
                        let shift = if row == col {
                            -0.5 * u * row.count_ones() as f64 + u * n as f64 / 4.
                        } else {
                            0.
                        };
                        close::<T>(sym[row * dim + col], std[row * dim + col] + shift);
                        if std[row * dim + col] != Complex::new(0., 0.) {
                            for spin in 0..2 {
                                let count = |bits: usize| {
                                    (0..n)
                                        .map(|j| (bits >> (2 * n - 1 - 2 * j - spin)) & 1)
                                        .sum::<usize>()
                                };
                                assert_eq!(count(row), count(col));
                            }
                        }
                    }
                }
            }
        }
    }
    let sites = sites(1);
    let chain = Chain::new(1, Boundary::Open);
    let sym = Fermion::<Complex<T>>::symmetric_hubbard(
        &chain,
        &sites,
        T::one(),
        T::from_f64(4.).unwrap(),
    )
    .unwrap();
    let leaves: Vec<_> = sites[0].leaves().collect();
    let double = Fermion::<Complex<T>>::state_tensor_g_at([
        (leaves[0], S::Occupied),
        (leaves[1], S::Occupied),
    ])
    .unwrap()
    .scaled(c::<T>(0., 2.));
    let normalized = double.normalized().unwrap();
    close::<T>(
        normalized.norm().unwrap().to_f64().unwrap().into(),
        1_f64.into(),
    );
    close::<T>(
        to64(normalized.normalized_expectation(&sym).unwrap()),
        1_f64.into(),
    );
    assert_eq!(
        normalized
            .coordinates_in(&OccupationBasis)
            .unwrap()
            .dimension()
            .unwrap(),
        4
    );
}
#[test]
fn models_reference_limits_boundaries_and_exports_f32() {
    check::<f32>();
}
#[test]
fn models_reference_limits_boundaries_and_exports_f64() {
    check::<f64>();
}

#[test]
fn invalid_sites_lengths_kinds_indices_and_error_sources() {
    type F = Fermion<Complex<f64>>;
    for symmetric in [false, true] {
        let build = if symmetric {
            F::symmetric_hubbard
        } else {
            F::hubbard
        };
        let chain = Chain::new(1, Boundary::Open);
        for bad in [
            H::unit(),
            H::factor(0),
            H::tensor_g([H::factor(0)]).unwrap(),
            H::tensor_g([H::factor(0), H::factor(1), H::factor(2)]).unwrap(),
            H::tensor([H::factor(0), H::factor(1)]).unwrap(),
        ] {
            assert!(matches!(
                build(&chain, &[bad], 0., 0.),
                Err(UtilityError::InvalidSite { site: 0 })
            ));
        }
        let mismatch = build(&chain, &[], 0., 0.).err().unwrap();
        assert_eq!(mismatch, UtilityError::Chain(ChainError::SiteCountMismatch));
        assert!(std::error::Error::source(&mismatch).is_some());
        let pair = sites(1).pop().unwrap();
        let duplicate = build(
            &Chain::new(2, Boundary::Open),
            &[pair.clone(), pair],
            0.,
            0.,
        )
        .err()
        .unwrap();
        assert_eq!(duplicate, UtilityError::Space(SpaceError::DuplicateIndex));
        assert!(std::error::Error::source(&duplicate).is_some());
    }
    let nested = H::tensor([H::unit(), sites(1).pop().unwrap()]).unwrap();
    assert!(F::hubbard(&Chain::new(1, Boundary::Open), &[nested], 1., 2.).is_ok());
}

fn graded<T: Real>() {
    let modes = [H::factor(9), H::factor(1)];
    let space = H::tensor_g(modes.clone()).unwrap();
    let annihilate: Vec<_> = modes
        .iter()
        .map(|h| {
            Fermion::<Complex<T>>::operator(O::Annihilate)
                .at(h)
                .unwrap()
                .bind(&space)
                .unwrap()
        })
        .collect();
    let dim = 4;
    for i in 0..2 {
        values::<T>(&exported_matrix(&annihilate[i]), &fock_matrix(i, 2, false));
        for j in 0..2 {
            let aa = annihilate[i].anticommutator(&annihilate[j]).unwrap();
            let ad = annihilate[i]
                .anticommutator(&annihilate[j].adjoint())
                .unwrap();
            values::<T>(&exported_matrix(&aa), &vec![Complex::new(0., 0.); 16]);
            let expected: Vec<_> = (0..16)
                .map(|k| Complex::new(if i == j && k / dim == k % dim { 1. } else { 0. }, 0.))
                .collect();
            values::<T>(&exported_matrix(&ad), &expected);
        }
    }
    let double = Fermion::<Complex<T>>::state_tensor_g_at([
        (&modes[0], S::Occupied),
        (&modes[1], S::Occupied),
    ])
    .unwrap();
    for (mode, labels, sign) in [
        (0, [S::Empty, S::Occupied], 1.),
        (1, [S::Occupied, S::Empty], -1.),
    ] {
        let target = Fermion::<Complex<T>>::state_tensor_g_at([
            (&modes[0], labels[0]),
            (&modes[1], labels[1]),
        ])
        .unwrap();
        close::<T>(
            to64(
                target
                    .inner_product(&double.act(&annihilate[mode]).unwrap())
                    .unwrap(),
            ),
            sign.into(),
        );
    }
    let parity = Fermion::<Complex<T>>::operator_product_at([
        (&modes[0], O::Parity),
        (&modes[1], O::Parity),
    ])
    .unwrap()
    .bind(&space)
    .unwrap();
    let m = exported_matrix(&parity);
    assert_eq!([m[0].re, m[5].re, m[10].re, m[15].re], [1., -1., -1., 1.]);
    let reversed = annihilate[1]
        .product(&annihilate[0])
        .unwrap()
        .scaled(c::<T>(0.3, -0.7));
    let reference = scale(
        &product(&fock_matrix(1, 2, false), &fock_matrix(0, 2, false)),
        Complex::new(0.3, -0.7),
    );
    let round = reversed
        .expand_in(&MajoranaBasis)
        .unwrap()
        .expand_in(&LadderBasis)
        .unwrap();
    values::<T>(&exported_matrix(&round), &reference);
    values::<T>(&exported_matrix(&round.adjoint()), &adjoint(&reference));
    assert_eq!(
        round
            .coordinates_in(&MajoranaBasis)
            .unwrap()
            .dimension()
            .unwrap(),
        16
    );
}
#[test]
fn car_spinful_signs_and_parity_f32() {
    graded::<f32>();
}
#[test]
fn car_spinful_signs_and_parity_f64() {
    graded::<f64>();
}
