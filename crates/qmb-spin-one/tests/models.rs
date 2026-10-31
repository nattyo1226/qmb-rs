mod support;

use num_complex::Complex;
use qmb_core::{ExpressionError, HilbertSpace as H, OperatorOf, SpaceError};
use qmb_export::{Export, Exported, MatrixArgs};
use qmb_geometry::{Boundary, Chain, MappingError};
use qmb_spin_one::{SpinOne, SpinOneOperator as O, SpinOneState as S, SzBasis, UtilityError};
use support::*;

type Spin<T> = SpinOne<Complex<T>>;

#[derive(Clone, Copy)]
enum Model {
    Xyz([f64; 3]),
    Heisenberg(f64),
    Blbq(f64, f64),
}

impl Model {
    fn build<T: Real>(self, chain: &Chain, sites: &[H<u8>]) -> OperatorOf<Spin<T>, u8> {
        let real = |x| T::from_f64(x).unwrap();
        match self {
            Self::Xyz([x, y, z]) => qmb_spin_one::models::SpinOneModels::<Spin<T>>::xyz(
                chain,
                sites,
                real(x),
                real(y),
                real(z),
            ),
            Self::Heisenberg(j) => {
                qmb_spin_one::models::SpinOneModels::<Spin<T>>::heisenberg(chain, sites, real(j))
            }
            Self::Blbq(j, k) => {
                qmb_spin_one::models::SpinOneModels::<Spin<T>>::blbq(chain, sites, real(j), real(k))
            }
        }
        .unwrap()
    }
}

// Independent dense construction from explicit spin matrices and leaf pairs.
// In particular, it uses neither model helpers nor Chain's edge enumeration.
fn reference(n: usize, pairs: &[(usize, usize)], model: Model) -> Matrix {
    let dimension = 3_usize.pow(n as u32);
    let mut result = vec![0.0.into(); dimension * dimension];
    for &(left, right) in pairs {
        let axes = [O::Sx, O::Sy, O::Sz].map(|axis| {
            (0..n).fold(vec![1.0.into()], |previous, position| {
                kron(
                    &previous,
                    &physical_matrix(if position == left || position == right {
                        axis
                    } else {
                        O::I
                    }),
                )
            })
        });
        let contribution = match model {
            Model::Xyz(couplings) => axes
                .iter()
                .zip(couplings)
                .fold(vec![0.0.into(); result.len()], |sum, (axis, coupling)| {
                    add(&sum, &scale(axis, (-coupling).into()))
                }),
            Model::Heisenberg(j) => scale(&add(&add(&axes[0], &axes[1]), &axes[2]), (-j).into()),
            Model::Blbq(j, k) => {
                let a = add(&add(&axes[0], &axes[1]), &axes[2]);
                add(
                    &scale(&a, (-j).into()),
                    &scale(&product(&a, &a), (-k).into()),
                )
            }
        };
        result = add(&result, &contribution);
    }
    result
}

fn exported<T: Real>(h: &OperatorOf<Spin<T>, u8>) -> Matrix {
    let dense: Exported<Vec<Complex<T>>, _> = h.export(MatrixArgs::<S, _>::new(&SzBasis)).unwrap();
    let dimension = h.space().unwrap().dim::<S>().unwrap();
    assert_eq!(dense.metadata().shape(), [dimension, dimension]);
    assert!(
        dense
            .metadata()
            .space()
            .shares_node_with(h.space().unwrap())
    );
    dense.data().iter().map(|z| to64(*z)).collect()
}

fn matches_reference<T: Real>(
    model: Model,
    chain: &Chain,
    sites: &[H<u8>],
    pairs: &[(usize, usize)],
) {
    let h = model.build::<T>(chain, sites);
    let n = h.space().unwrap().indices().count();
    let expected = reference(n, pairs, model);
    let actual = exported(&h);
    values::<T>(&actual, &expected);
    values::<T>(&actual, &adjoint(&actual));
    close::<T>(to64(h.trace().unwrap()), trace(&expected));
}

fn small_chains<T: Real>() {
    for boundary in [Boundary::Open, Boundary::Periodic] {
        for n in 0..=3 {
            let sites: Vec<_> = [11, 7, 9].into_iter().take(n).map(H::factor).collect();
            let chain = Chain::new(n, boundary);
            let pairs = match (n, boundary) {
                (0 | 1, _) => vec![],
                (2, _) => vec![(0, 1)],
                (3, Boundary::Open) => vec![(0, 1), (1, 2)],
                (3, Boundary::Periodic) => vec![(0, 1), (0, 2), (1, 2)],
                _ => unreachable!(),
            };
            for model in [
                Model::Xyz([0.7, -0.4, 1.2]),
                Model::Xyz([0.0, 0.0, -0.8]),
                Model::Xyz([0.0; 3]),
                Model::Heisenberg(1.0),
                Model::Heisenberg(-0.6),
                Model::Blbq(1.0, 0.5),
                Model::Blbq(-0.7, -0.3),
                Model::Blbq(0.0, 0.8),
                Model::Blbq(0.0, 0.0),
            ] {
                matches_reference::<T>(model, &chain, &sites, &pairs);
            }

            let isotropic = exported(&Model::Heisenberg(-0.7).build::<T>(&chain, &sites));
            values::<T>(
                &isotropic,
                &exported(&Model::Xyz([-0.7; 3]).build::<T>(&chain, &sites)),
            );
            values::<T>(
                &isotropic,
                &exported(&Model::Blbq(-0.7, 0.0).build::<T>(&chain, &sites)),
            );
        }
    }
}

#[test]
fn matrices_match_independent_small_chains_f32() {
    small_chains::<f32>();
}

#[test]
fn matrices_match_independent_small_chains_f64() {
    small_chains::<f64>();
}

fn total_spin<T: Real>() {
    let sites = [H::factor(11), H::factor(7)];
    let chain = Chain::new(2, Boundary::Open);
    let space = H::tensor_all(sites.iter().cloned()).unwrap();
    let singlet = vec![
        (1.0 / 3_f64.sqrt(), S::Plus, S::Minus),
        (-1.0 / 3_f64.sqrt(), S::Zero, S::Zero),
        (1.0 / 3_f64.sqrt(), S::Minus, S::Plus),
    ];
    let triplet = vec![
        (1.0 / 2_f64.sqrt(), S::Plus, S::Zero),
        (-1.0 / 2_f64.sqrt(), S::Zero, S::Plus),
    ];
    let quintet = vec![(1.0, S::Plus, S::Plus)];
    // A = S_1 dot S_2 has eigenvalues -2, -1, 1 at total spin 0, 1, 2.
    for (a, terms) in [(-2.0, singlet), (-1.0, triplet), (1.0, quintet)] {
        let psi = Spin::<T>::state_sum(
            &space,
            terms.into_iter().map(|(weight, left, right)| {
                Spin::<T>::product_state(
                    &qmb_core::HilbertSpace::tensor_all([(sites[0]).clone(), (sites[1]).clone()])
                        .unwrap(),
                    [left, right],
                )
                .unwrap()
                .scaled(c::<T>(weight, 0.0))
            }),
        )
        .unwrap();
        close::<T>(to64(psi.inner_product(&psi).unwrap()), 1.0.into());
        let vector: Exported<Vec<Complex<T>>, _> = psi.export(&SzBasis).unwrap();
        let vector: Vec<_> = vector.data().iter().map(|z| to64(*z)).collect();
        for (j, k) in [(-1.0, 0.0), (1.0, 0.5), (-0.7, -0.3), (-1.0, -1.0 / 3.0)] {
            let h = Model::Blbq(j, k).build::<T>(&chain, &sites);
            let energy = -j * a - k * a * a;
            close::<T>(to64(psi.expectation(&h).unwrap()), energy.into());
            let acted: Exported<Vec<Complex<T>>, _> =
                psi.act(&h).unwrap().export(&SzBasis).unwrap();
            values::<T>(
                &acted.data().iter().map(|z| to64(*z)).collect::<Vec<_>>(),
                &scale(&vector, energy.into()),
            );
        }
    }
}

#[test]
fn blbq_total_spin_energies_and_eigenstates_f32() {
    total_spin::<f32>();
}

#[test]
fn blbq_total_spin_energies_and_eigenstates_f64() {
    total_spin::<f64>();
}

fn site_spaces<T: Real>() {
    let models = [
        Model::Xyz([0.7, -0.4, 1.2]),
        Model::Heisenberg(-0.6),
        Model::Blbq(0.4, -0.7),
    ];
    let a = H::factor(11);
    let b = H::factor(7);
    let c = H::factor(9);
    let first = H::tensor_all([H::unit(), H::tensor_all([a.clone(), b.clone()]).unwrap()]).unwrap();
    let sites = [first, c];
    let chain = Chain::new(2, Boundary::Periodic);
    for model in models {
        // Each of a and b couples to c. There is no a-b interaction, and BLBQ
        // squares the two pair interactions separately, omitting cross-pair terms.
        matches_reference::<T>(model, &chain, &sites, &[(0, 2), (1, 2)]);
        let h = model.build::<T>(&chain, &sites);
        let full = h.space().unwrap();
        assert_eq!(full.indices().copied().collect::<Vec<_>>(), [11, 7, 9]);
        assert_eq!(full.children().count(), 2);
        for (actual, input) in full.children().zip(&sites) {
            assert!(actual.shares_node_with(input));
        }
    }

    let sites = [a.clone(), H::unit(), b.clone()];
    for boundary in [Boundary::Open, Boundary::Periodic] {
        let chain = Chain::new(3, boundary);
        let pairs = if boundary == Boundary::Periodic {
            vec![(0, 1)]
        } else {
            vec![]
        };
        for model in models {
            matches_reference::<T>(model, &chain, &sites, &pairs);
            let h = model.build::<T>(&chain, &sites);
            assert_eq!(h.space().unwrap().children().count(), 3);
            for (actual, input) in h.space().unwrap().children().zip(&sites) {
                assert!(actual.shares_node_with(input));
            }
        }
    }

    for n in 0..=3 {
        let sites = vec![H::unit(); n];
        for boundary in [Boundary::Open, Boundary::Periodic] {
            for model in models {
                let h = model.build::<T>(&Chain::new(n, boundary), &sites);
                values::<T>(&exported(&h), &[0.0.into()]);
                assert_eq!(h.space().unwrap().children().count(), n);
            }
        }
    }
    // A single composite site is a spectator; no within-site terms are added.
    let sites = [H::tensor_all([a, b]).unwrap()];
    for model in models {
        matches_reference::<T>(model, &Chain::new(1, Boundary::Periodic), &sites, &[]);
        let h = model.build::<T>(&Chain::new(1, Boundary::Open), &sites);
        assert_eq!(h.space().unwrap().dim::<S>().unwrap(), 9);
        assert!(
            h.space()
                .unwrap()
                .children()
                .next()
                .unwrap()
                .shares_node_with(&sites[0])
        );
    }
}

#[test]
fn leaf_pair_squares_empty_sites_and_shared_grouping_f32() {
    site_spaces::<f32>();
}

#[test]
fn leaf_pair_squares_empty_sites_and_shared_grouping_f64() {
    site_spaces::<f64>();
}

#[test]
fn every_model_validates_even_with_no_edges_and_zero_couplings() {
    type Q = Spin<f64>;
    for boundary in [Boundary::Open, Boundary::Periodic] {
        for n in 0..=3 {
            for count in 0..=4 {
                let sites: Vec<_> = (0..count).map(H::factor).collect();
                let chain = Chain::new(n, boundary);
                let expected = (n != count).then_some(UtilityError::Mapping(
                    MappingError::SiteCountMismatch {
                        expected: n,
                        actual: count,
                    },
                ));
                assert_eq!(
                    qmb_spin_one::models::SpinOneModels::<Q>::xyz(&chain, &sites, 0.0, 0.0, 0.0)
                        .err(),
                    expected
                );
                assert_eq!(
                    qmb_spin_one::models::SpinOneModels::<Q>::heisenberg(&chain, &sites, 0.0).err(),
                    expected
                );
                assert_eq!(
                    qmb_spin_one::models::SpinOneModels::<Q>::blbq(&chain, &sites, 0.0, 0.0).err(),
                    expected
                );
            }
        }
        let sites = [
            H::factor(0),
            H::unit(),
            H::tensor_all([H::factor(0)]).unwrap(),
        ];
        let chain = Chain::new(3, boundary);
        let expected = Some(UtilityError::Space(SpaceError::DuplicateIndex));
        assert_eq!(
            qmb_spin_one::models::SpinOneModels::<Q>::xyz(&chain, &sites, 0.0, 0.0, 0.0).err(),
            expected
        );
        assert_eq!(
            qmb_spin_one::models::SpinOneModels::<Q>::heisenberg(&chain, &sites, 0.0).err(),
            expected
        );
        assert_eq!(
            qmb_spin_one::models::SpinOneModels::<Q>::blbq(&chain, &sites, 0.0, 0.0).err(),
            expected
        );
    }
}

#[test]
fn construction_needs_rational_arithmetic_but_no_float_or_evaluation() {
    #[derive(PartialEq, Eq)]
    struct Index(u8);
    struct Input(i32);
    #[derive(Clone)]
    struct Coefficient(f64);
    impl core::ops::Add for Coefficient {
        type Output = Self;
        fn add(self, rhs: Self) -> Self {
            Self(self.0 + rhs.0)
        }
    }

    impl core::ops::Mul for Coefficient {
        type Output = Self;
        fn mul(self, rhs: Self) -> Self {
            Self(self.0 * rhs.0)
        }
    }

    impl core::ops::Div for Coefficient {
        type Output = Self;
        fn div(self, rhs: Self) -> Self {
            Self(self.0 / rhs.0)
        }
    }
    impl num_traits::One for Coefficient {
        fn one() -> Self {
            Self(1.0)
        }
    }

    impl From<Input> for Coefficient {
        fn from(value: Input) -> Self {
            Self(f64::from(value.0))
        }
    }
    impl From<i32> for Coefficient {
        fn from(value: i32) -> Self {
            Self(f64::from(value))
        }
    }
    impl core::ops::Neg for Coefficient {
        type Output = Self;
        fn neg(self) -> Self {
            Self(-self.0)
        }
    }
    // No Zero, Float, primitive rules or index Clone/Hash/Ord.
    let sites = [H::factor(Index(11)), H::factor(Index(7))];
    let chain = Chain::new(2, Boundary::Open);
    type Q = SpinOne<Coefficient>;
    for h in [
        qmb_spin_one::models::SpinOneModels::<Q>::xyz(&chain, &sites, Input(1), Input(2), Input(3))
            .unwrap(),
        qmb_spin_one::models::SpinOneModels::<Q>::heisenberg(&chain, &sites, 1).unwrap(),
        qmb_spin_one::models::SpinOneModels::<Q>::blbq(&chain, &sites, Input(1), Input(2)).unwrap(),
    ] {
        assert_eq!(
            h.space()
                .unwrap()
                .indices()
                .map(|i| i.0)
                .collect::<Vec<_>>(),
            [11, 7]
        );
        assert!(
            h.space()
                .unwrap()
                .children()
                .next()
                .unwrap()
                .shares_node_with(&sites[0])
        );
    }
}

#[test]
fn model_errors_preserve_the_original_source_and_message() {
    for (error, source) in [
        (
            UtilityError::from(MappingError::SiteCountMismatch {
                expected: 2,
                actual: 1,
            }),
            MappingError::SiteCountMismatch {
                expected: 2,
                actual: 1,
            }
            .to_string(),
        ),
        (
            UtilityError::from(SpaceError::DuplicateIndex),
            SpaceError::DuplicateIndex.to_string(),
        ),
        (
            UtilityError::from(ExpressionError::NotAFactor),
            ExpressionError::NotAFactor.to_string(),
        ),
    ] {
        assert_eq!(error.to_string(), source);
        assert_eq!(
            std::error::Error::source(&error).unwrap().to_string(),
            source
        );
    }
}
