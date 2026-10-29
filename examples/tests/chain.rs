use num_complex::Complex;
#[cfg(feature = "normalize")]
use num_traits::Float;
use num_traits::{FromPrimitive, Num, ToPrimitive};
use qmb_basic::Id128;
use qmb_core::{Divide, HilbertSpace, LocalState, State};
use qmb_geometry::{Boundary, Chain};
#[cfg(feature = "normalize")]
use qmb_normalize::{Norm, Normalize, NormalizedExpectation};
use qmb_spin_half::{SpinHalf, SpinHalfState};
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
    factors: &[HilbertSpace<Id128>],
    bits: usize,
) -> qmb_core::StateOf<SpinHalf<Complex<T>>, Id128> {
    State::tensor_in(
        &HilbertSpace::tensor(factors.iter().cloned()).unwrap(),
        factors.iter().enumerate().map(|(i, factor)| {
            // Independent equal-index handles need not share the original nodes.
            LocalState::primitive(if (bits >> (factors.len() - 1 - i)) & 1 == 0 {
                SpinHalfState::Up
            } else {
                SpinHalfState::Down
            })
            .at(&HilbertSpace::factor(*factor.factor_index().unwrap()))
            .unwrap()
        }),
    )
    .unwrap()
}
fn factors(n: usize) -> Vec<HilbertSpace<Id128>> {
    (0..n).map(|_| HilbertSpace::factor(Id128::new())).collect()
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
            let chain = Chain::new(n, boundary);
            let factors = factors(n);
            let space = HilbertSpace::tensor(factors.iter().cloned()).unwrap();
            let j = real::<T>(1.3);
            let h = SpinHalf::<Complex<T>>::heisenberg(&chain, &factors, j).unwrap();
            assert!(h.space().unwrap().is_compatible_with(&space));
            let before: Vec<_> = space.indices().copied().collect();
            for first in [SpinHalfState::Up, SpinHalfState::Down] {
                let psi = SpinHalf::<Complex<T>>::neel(&chain, &factors, first).unwrap();
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
            let chain = Chain::new(n, boundary);
            let factors = factors(n);
            for j in [real::<T>(1.3), real::<T>(-0.7)] {
                let h = SpinHalf::<Complex<T>>::heisenberg(&chain, &factors, j).unwrap();
                let adj = h.adjoint();
                for row in 0..(1 << n) {
                    for col in 0..(1 << n) {
                        let left = basis::<T>(&factors, row);
                        let right = basis::<T>(&factors, col);
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

#[cfg(feature = "export")]
fn check_exports<T: Real + num_traits::NumAssign + Default + Send + Sync + 'static>() {
    use qmb_export::{
        Export, Exported, MatrixArgs,
        ndarray::{Array1, Array2},
        sprs::{CsMat, CsVec},
    };
    use qmb_spin_half::SzBasis;
    for n in 0..=4 {
        for boundary in [Boundary::Open, Boundary::Periodic] {
            let chain = Chain::new(n, boundary);
            let factors = factors(n);
            let space = HilbertSpace::tensor(factors.iter().cloned()).unwrap();
            let j = real::<T>(1.3);
            let h = SpinHalf::<Complex<T>>::heisenberg(&chain, &factors, j).unwrap();
            for first in [SpinHalfState::Up, SpinHalfState::Down] {
                let psi = SpinHalf::<Complex<T>>::neel(&chain, &factors, first).unwrap();
                let raw: Exported<Vec<Complex<T>>, _> = (&psi).export(&SzBasis).unwrap();
                let pointer = raw.data().as_ptr();
                let vector: Exported<Array1<_>, _> = raw.export(()).unwrap();
                assert_eq!(pointer, vector.data().as_ptr());
                let direct: Exported<Array1<_>, _> = (&psi).export(&SzBasis).unwrap();
                assert_eq!(vector.data(), direct.data());
                let dense: Exported<Array2<_>, _> = (&h)
                    .export(MatrixArgs::<SpinHalfState, _>::new(&SzBasis))
                    .unwrap();
                let raw: Exported<Vec<Complex<T>>, _> = (&h)
                    .export(MatrixArgs::<SpinHalfState, _>::new(&SzBasis))
                    .unwrap();
                assert_eq!(raw.data().as_slice(), dense.data().as_slice().unwrap());
                assert_eq!(raw.metadata().shape(), [1 << n; 2]);
                let sparse: Exported<CsMat<_>, _> = (&h)
                    .export(MatrixArgs::<SpinHalfState, _>::new(&SzBasis))
                    .unwrap();
                let sparse_vector: Exported<CsVec<_>, _> = (&psi).export(&SzBasis).unwrap();
                assert!(sparse.data().is_csc());
                assert!(dense.metadata().space().is_compatible_with(&space));
                assert!(sparse.metadata().space().is_compatible_with(&space));
                let bits = (0..n).fold(0, |bits, i| {
                    (bits << 1) | usize::from((i % 2 == 1) == (first == SpinHalfState::Up))
                });
                for i in 0..1 << n {
                    close(vector.data()[i], c(if i == bits { 1. } else { 0. }, 0.));
                    for k in 0..1 << n {
                        let reference = matrix_element(n, boundary, i, k, [j; 3]);
                        close(dense.data()[[i, k]], reference);
                        close(
                            sparse
                                .data()
                                .get(i, k)
                                .copied()
                                .unwrap_or_else(|| c(0., 0.)),
                            reference,
                        );
                        close(dense.data()[[i, k]], dense.data()[[k, i]].conj());
                    }
                }
                let dense_acted = dense.data().dot(vector.data());
                let sparse_acted = sparse.data() * sparse_vector.data();
                let algebra = psi.act(&h).unwrap();
                let algebra: Exported<Vec<Complex<T>>, _> = (&algebra).export(&SzBasis).unwrap();
                for i in 0..1 << n {
                    close(dense_acted[i], algebra.data()[i]);
                    close(
                        sparse_acted.get(i).copied().unwrap_or_else(|| c(0., 0.)),
                        algebra.data()[i],
                    );
                }
                let energy_dense = vector
                    .data()
                    .iter()
                    .zip(&dense_acted)
                    .fold(c(0., 0.), |s, (&x, &y)| s + x.conj() * y);
                let energy_sparse = sparse_vector.data().iter().fold(c(0., 0.), |s, (i, &x)| {
                    s + x.conj() * sparse_acted.get(i).copied().unwrap_or_else(|| c(0., 0.))
                });
                close(energy_dense, psi.expectation(&h).unwrap());
                close(energy_sparse, energy_dense);
                close(
                    h.trace().unwrap(),
                    (0..1 << n).fold(c(0., 0.), |s, i| s + dense.data()[[i, i]]),
                );
                if boundary == Boundary::Open {
                    close(
                        energy_dense,
                        Complex::new(
                            if n == 0 {
                                T::zero()
                            } else {
                                j * real::<T>((n - 1) as f64 / 4.)
                            },
                            T::zero(),
                        ),
                    );
                }
            }
        }
    }
}

#[cfg(feature = "export")]
#[test]
fn exports_complete_chain_scenario_at_both_precisions() {
    check_exports::<f32>();
    check_exports::<f64>();
}

#[cfg(feature = "normalize")]
fn check_norms<T: Real + Float>() {
    for n in 0..=6 {
        for first in [SpinHalfState::Up, SpinHalfState::Down] {
            let chain = Chain::new(n, Boundary::Open);
            let factors = factors(n);
            let psi = SpinHalf::<Complex<T>>::neel(&chain, &factors, first).unwrap();
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
fn normalization_extensions_work_at_both_precisions() {
    check_norms::<f32>();
    check_norms::<f64>();
}

fn check_xyz<T: Real>() {
    for n in 0..=3 {
        for boundary in [Boundary::Open, Boundary::Periodic] {
            let chain = Chain::new(n, boundary);
            let factors = factors(n);
            let space = HilbertSpace::tensor(factors.iter().cloned()).unwrap();
            for values in [
                [0.7, -1.3, 2.1],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 0.0],
                [0.8, 0.8, 0.8],
            ] {
                let [jx, jy, jz] = values.map(real::<T>);
                let h = SpinHalf::<Complex<T>>::xyz(&chain, &factors, jx, jy, jz).unwrap();
                assert!(h.space().unwrap().is_compatible_with(&space));
                let isotropic = SpinHalf::<Complex<T>>::heisenberg(&chain, &factors, jx).unwrap();
                for row in 0..(1 << n) {
                    for col in 0..(1 << n) {
                        let left = basis::<T>(&factors, row);
                        let right = basis::<T>(&factors, col);
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
