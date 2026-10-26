use num_complex::Complex;
#[cfg(feature = "normalize")]
use num_traits::Float;
use num_traits::{FromPrimitive, Num, ToPrimitive};
use qmb::spin_half::utility::{UtilityError, heisenberg, neel, xyz};
use qmb::{
    Boundary, Chain, Divide, HilbertSpace, Id128, LocalState, Operator, SpinHalfState, State,
    ToSpace,
};
#[cfg(feature = "normalize")]
use qmb::{Norm, Normalize, NormalizedExpectation};
use std::{fmt::Debug, ops::Neg};

trait Real: Num + Copy + Neg<Output = Self> + FromPrimitive + ToPrimitive + Debug {
    const ATOL: f64;
    const RTOL: f64;
}
impl Real for f32 {
    const ATOL: f64 = 2e-5;
    const RTOL: f64 = 2e-5;
}
impl Real for f64 {
    const ATOL: f64 = 1e-12;
    const RTOL: f64 = 1e-12;
}
fn real<T: Real>(value: f64) -> T {
    T::from_f64(value).unwrap()
}
fn c<T: Real>(re: f64, im: f64) -> Complex<T> {
    Complex::new(real(re), real(im))
}
fn close<T: Real>(actual: Complex<T>, expected: Complex<T>) {
    for (a, b) in [(actual.re, expected.re), (actual.im, expected.im)] {
        let (a, b) = (a.to_f64().unwrap(), b.to_f64().unwrap());
        assert!(
            (a - b).abs() <= T::ATOL + T::RTOL * b.abs(),
            "{actual:?} != {expected:?}"
        );
    }
}
fn basis<T: Real>(
    chain: &Chain<Id128>,
    bits: usize,
) -> State<Id128, qmb_spin_half::SpinHalfState, qmb_spin_half::SpinHalfOperator, Complex<T>> {
    State::tensor(chain.sites().iter().enumerate().map(|(i, site)| {
        // Build the same indexed factors without sharing Chain storage.
        LocalState::primitive(if (bits >> (chain.sites().len() - 1 - i)) & 1 == 0 {
            SpinHalfState::Up
        } else {
            SpinHalfState::Down
        })
        .at(&HilbertSpace::factor(site[0]))
        .unwrap()
    }))
    .unwrap()
}
fn explicit_chain(n: usize, boundary: Boundary) -> Chain<Id128> {
    Chain::from_sites((0..n).map(|_| [Id128::new()]), boundary).unwrap()
}

// Independent reference: Kronecker products of the three Pauli/2 matrices.
// This deliberately does not use chain edges or primitive algebra rules.
fn matrix_element<T: Real>(
    n: usize,
    boundary: Boundary,
    row: usize,
    col: usize,
    couplings: [T; 3],
) -> Complex<T> {
    let z = c::<T>(0.0, 0.0);
    let half = real::<T>(0.5);
    let spins = [
        [[z, c(0.5, 0.0)], [c(0.5, 0.0), z]],
        [[z, c(0.0, -0.5)], [c(0.0, 0.5), z]],
        [
            [Complex::new(half, T::zero()), z],
            [z, Complex::new(-half, T::zero())],
        ],
    ];
    let mut edges: Vec<_> = (1..n).map(|i| (i - 1, i)).collect();
    if boundary == Boundary::Periodic && n > 2 {
        edges.push((0, n - 1));
    }
    let mut value = z;
    for (a, b) in edges {
        for (spin, coupling) in spins.into_iter().zip(couplings) {
            let mut term = c::<T>(1.0, 0.0);
            for i in 0..n {
                let r = (row >> (n - 1 - i)) & 1;
                let s = (col >> (n - 1 - i)) & 1;
                term = term
                    * if i == a || i == b {
                        spin[r][s]
                    } else if r == s {
                        c(1.0, 0.0)
                    } else {
                        z
                    };
            }
            value = value - term * Complex::new(coupling, T::zero());
        }
    }
    value
}

fn check<T: Real>() {
    for n in 0..=6 {
        for boundary in [Boundary::Open, Boundary::Periodic] {
            let chain = Chain::<Id128>::new(n, 1, boundary).unwrap();
            let space = chain.to_space().unwrap();
            assert!(chain.sites().iter().all(|s| s.len() == 1));
            let j = real::<T>(1.3);
            let h = heisenberg(&chain, &space, j).unwrap();
            assert!(h.space().unwrap().shares_node_with(&space));
            let before: Vec<_> = space.indices().copied().collect();
            for first in [SpinHalfState::Up, SpinHalfState::Down] {
                let psi = neel::<_, Complex<T>>(&chain, &space, first).unwrap();
                close(
                    Complex::new(psi.norm_squared().unwrap(), T::zero()),
                    c(1.0, 0.0),
                );
                let expected = if n == 0 {
                    T::zero()
                } else {
                    let mut e = j * real::<T>((n - 1) as f64 / 4.0);
                    if boundary == Boundary::Periodic && n > 2 {
                        e = e + j * real::<T>(if n % 2 == 0 { 0.25 } else { -0.25 });
                    }
                    e
                };
                close(
                    psi.expectation(&h).unwrap(),
                    Complex::new(expected, T::zero()),
                );
                close(
                    psi.expectation(&h)
                        .unwrap()
                        .divide(psi.inner_product(&psi).unwrap()),
                    Complex::new(expected, T::zero()),
                );
                let scaled = psi.scaled(c(2.0, 1.0));
                close(
                    scaled
                        .expectation(&h)
                        .unwrap()
                        .divide(scaled.inner_product(&scaled).unwrap()),
                    Complex::new(expected, T::zero()),
                );
                #[cfg(feature = "normalize")]
                {
                    close(
                        psi.normalized_expectation(&h).unwrap(),
                        Complex::new(expected, T::zero()),
                    );
                    close(
                        scaled.normalized_expectation(&h).unwrap(),
                        Complex::new(expected, T::zero()),
                    );
                }
                if n == 0 {
                    close(psi.inner_product(&State::unit()).unwrap(), c(1.0, 0.0));
                    close(h.trace().unwrap(), c(0.0, 0.0));
                }
            }
            assert_eq!(space.indices().copied().collect::<Vec<_>>(), before);
        }
    }
    for n in 0..=3 {
        for boundary in [Boundary::Open, Boundary::Periodic] {
            let chain = explicit_chain(n, boundary);
            let space = chain.to_space().unwrap();
            for j in [real::<T>(1.3), real::<T>(-0.7)] {
                let h = heisenberg(&chain, &space, j).unwrap();
                let adj = h.adjoint();
                for row in 0..(1 << n) {
                    for col in 0..(1 << n) {
                        let left = basis::<T>(&chain, row);
                        let right = basis::<T>(&chain, col);
                        let value = left.inner_product(&right.act(&h).unwrap()).unwrap();
                        close(value, matrix_element(n, boundary, row, col, [j; 3]));
                        close(
                            value,
                            left.inner_product(&right.act(&adj).unwrap()).unwrap(),
                        );
                        let reverse = right.inner_product(&left.act(&h).unwrap()).unwrap();
                        close(value, Complex::new(reverse.re, -reverse.im));
                    }
                }
            }
        }
    }
}
#[test]
fn models_match_analytic_energies_and_independent_matrices_both_precisions() {
    check::<f32>();
    check::<f64>();
}

#[cfg(feature = "normalize")]
fn check_norms<T: Real + Float>() {
    for n in 0..=6 {
        for first in [SpinHalfState::Up, SpinHalfState::Down] {
            let chain = Chain::<Id128>::new(n, 1, Boundary::Open).unwrap();
            let space = chain.to_space().unwrap();
            let psi = neel::<_, Complex<T>>(&chain, &space, first).unwrap();
            close(Complex::new(psi.norm().unwrap(), T::zero()), c(1.0, 0.0));
            let normalized = psi.scaled(c(2.0, -3.0)).normalized().unwrap();
            close(
                Complex::new(normalized.norm().unwrap(), T::zero()),
                c(1.0, 0.0),
            );
        }
    }
}
#[cfg(feature = "normalize")]
#[test]
fn default_facade_enables_norms_at_both_precisions() {
    check_norms::<f32>();
    check_norms::<f64>();
}

#[test]
fn utilities_reject_empty_and_multiple_factors_but_accept_nested_single_factors() {
    for k in [0, 2] {
        let chain = Chain::<Id128>::new(2, k, Boundary::Open).unwrap();
        let space = chain.to_space().unwrap();
        assert_eq!(
            neel::<_, Complex<f64>>(&chain, &space, SpinHalfState::Up).err(),
            Some(UtilityError::InvalidSite(0))
        );
        assert_eq!(
            heisenberg::<_, _, Complex<f64>>(&chain, &space, 1.0_f64).err(),
            Some(UtilityError::InvalidSite(0))
        );
        assert_eq!(
            xyz::<_, _, Complex<f64>>(&chain, &space, 1.0_f64, 2.0, -3.0).err(),
            Some(UtilityError::InvalidSite(0))
        );
    }
    let a = HilbertSpace::<_>::factor(Id128::new());
    let b = HilbertSpace::factor(Id128::new());
    let nested = HilbertSpace::tensor([
        HilbertSpace::unit(),
        HilbertSpace::tensor([a.clone()]).unwrap(),
    ])
    .unwrap();
    let full = HilbertSpace::tensor([
        HilbertSpace::tensor([b.clone(), HilbertSpace::unit()]).unwrap(),
        a,
    ])
    .unwrap();
    let chain = Chain::from_sites(
        [
            nested.indices().copied().collect::<Vec<_>>(),
            b.indices().copied().collect(),
        ],
        Boundary::Open,
    )
    .unwrap();
    let space = full.clone();
    let psi = neel::<_, Complex<f64>>(&chain, &space, SpinHalfState::Up).unwrap();
    assert!(psi.space().is_compatible_with(&full));
    close(
        psi.expectation(&heisenberg(&chain, &space, 1.0).unwrap())
            .unwrap(),
        c(0.25, 0.0),
    );
    // Explicit full order is Down then Up, although site order starts Up.
    let expected = State::tensor(
        full.indices()
            .copied()
            .zip([SpinHalfState::Down, SpinHalfState::Up])
            .map(|(index, k)| {
                LocalState::primitive(k)
                    .at(&HilbertSpace::factor(index))
                    .unwrap()
            }),
    )
    .unwrap();
    close(psi.inner_product(&expected).unwrap(), c(1.0, 0.0));
    let zero = Operator::<_, qmb_spin_half::SpinHalfOperator, Complex<f64>>::zero()
        .bind(&full)
        .unwrap();
    close(psi.expectation(&zero).unwrap(), c(0.0, 0.0));
    assert!(space.shares_node_with(&full));
    assert_eq!(full.children().next().unwrap().children().count(), 2);
}

#[test]
fn model_utilities_reject_missing_extra_and_wrong_space_indices() {
    let chain = Chain::from_sites([[0u8], [1]], Boundary::Open).unwrap();
    for indices in [vec![], vec![0], vec![0, 1, 2], vec![0, 2]] {
        let space = HilbertSpace::tensor(indices.into_iter().map(HilbertSpace::factor)).unwrap();
        assert_eq!(
            neel::<_, Complex<f64>>(&chain, &space, SpinHalfState::Up).err(),
            Some(UtilityError::IncompatibleSpace)
        );
        assert_eq!(
            heisenberg::<_, _, Complex<f64>>(&chain, &space, 1.0).err(),
            Some(UtilityError::IncompatibleSpace)
        );
        assert_eq!(
            xyz::<_, _, Complex<f64>>(&chain, &space, 1.0, 2.0, 3.0).err(),
            Some(UtilityError::IncompatibleSpace)
        );
    }
    let empty = Chain::<u8>::from_sites(Vec::<Vec<u8>>::new(), Boundary::Open).unwrap();
    let space = HilbertSpace::factor(0);
    assert_eq!(
        neel::<_, Complex<f64>>(&empty, &space, SpinHalfState::Down).err(),
        Some(UtilityError::IncompatibleSpace)
    );
}

#[test]
fn existing_space_models_need_no_index_clone_and_walk_deep_grouping_iteratively() {
    // No Clone, Default, Copy, Hash or Ord. Equal values denote the same identity.
    #[derive(PartialEq, Eq)]
    struct Index(u8);
    std::thread::Builder::new()
        .stack_size(128 * 1024)
        .spawn(|| {
            let chain = Chain::from_sites([[Index(0)], [Index(1)]], Boundary::Periodic).unwrap();
            let mut space = HilbertSpace::tensor([
                HilbertSpace::factor(Index(1)),
                HilbertSpace::factor(Index(0)),
            ])
            .unwrap();
            for _ in 0..20_000 {
                space = HilbertSpace::tensor([space]).unwrap();
            }
            let state = neel::<_, Complex<f64>>(&chain, &space, SpinHalfState::Up).unwrap();
            let h = heisenberg(&chain, &space, 1.0_f64).unwrap();
            let isotropic = xyz(&chain, &space, 1.0_f64, 1.0, 1.0).unwrap();
            assert!(h.space().unwrap().shares_node_with(&space));
            assert!(state.space().is_compatible_with(&space));
            close(state.expectation(&h).unwrap(), c(0.25, 0.0));
            close(state.expectation(&isotropic).unwrap(), c(0.25, 0.0));
            assert_eq!(chain.site(0).unwrap()[0].0, 0);
        })
        .unwrap()
        .join()
        .unwrap();
}

fn check_xyz<T: Real>() {
    for n in 0..=3 {
        for boundary in [Boundary::Open, Boundary::Periodic] {
            let chain = explicit_chain(n, boundary);
            let space = chain.to_space().unwrap();
            for values in [
                [0.7, -1.3, 2.1],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 0.0],
                [0.8, 0.8, 0.8],
            ] {
                let [jx, jy, jz] = values.map(real::<T>);
                let h = xyz(&chain, &space, jx, jy, jz).unwrap();
                assert!(h.space().unwrap().shares_node_with(&space));
                let isotropic = heisenberg(&chain, &space, jx).unwrap();
                for row in 0..(1 << n) {
                    for col in 0..(1 << n) {
                        let left = basis::<T>(&chain, row);
                        let right = basis::<T>(&chain, col);
                        let value = left.inner_product(&right.act(&h).unwrap()).unwrap();
                        close(value, matrix_element(n, boundary, row, col, [jx, jy, jz]));
                        let reverse = right.inner_product(&left.act(&h).unwrap()).unwrap();
                        close(value, Complex::new(reverse.re, -reverse.im));
                        if values[0] == values[1] && values[1] == values[2] {
                            close(
                                value,
                                left.inner_product(&right.act(&isotropic).unwrap()).unwrap(),
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn xyz_matches_anisotropic_matrices_and_heisenberg_limit_both_precisions() {
    check_xyz::<f32>();
    check_xyz::<f64>();
}

#[test]
fn utility_error_conversion_preserves_message_and_source() {
    let original = qmb::ExpressionError::NotAFactor;
    let error = UtilityError::from(original);
    assert_eq!(error, UtilityError::Expression(original));
    assert_eq!(error.to_string(), original.to_string());
    assert_eq!(
        std::error::Error::source(&error).unwrap().to_string(),
        original.to_string()
    );
    let invalid = UtilityError::InvalidSite(2);
    assert_eq!(
        invalid.to_string(),
        "site 2 must have exactly one spin factor"
    );
    assert!(std::error::Error::source(&invalid).is_none());
    assert_eq!(
        UtilityError::IncompatibleSpace.to_string(),
        "site indices do not partition the full space"
    );
    assert!(std::error::Error::source(&UtilityError::IncompatibleSpace).is_none());
}
