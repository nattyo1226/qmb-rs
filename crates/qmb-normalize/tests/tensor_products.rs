#[path = "../../qmb-core/tests/support/graded.rs"]
mod graded;
use graded::*;
use qmb_core::{HilbertSpace, TensorKind};
use qmb_normalize::{Norm, Normalize, NormalizedExpectation, NormalizedExpectationIn};
fn check<T: Real + num_traits::Float>() {
    let factors = [
        HilbertSpace::factor(0),
        HilbertSpace::factor(1),
        HilbertSpace::factor(2),
    ];
    let h = HilbertSpace::tensor_all_g(factors.clone()).unwrap();
    let op = placed::<T>(A, &factors[2]);
    let psi = (ket(7, &factors, TensorKind::Graded)
        .scaled(c(0.3, -0.6))
        .clone()
        + &ket(6, &factors, TensorKind::Graded).scaled(c(-0.7, 0.4)))
        .unwrap();
    let m = embedding(A, 2, 3, TensorKind::Graded);
    let v = vector(&psi);
    let av = act(&m, &v);
    let inner = |a: &[num_complex::Complex<T>], b: &[num_complex::Complex<T>]| {
        a.iter()
            .zip(b)
            .fold(c(0., 0.), |s, (a, b)| s + a.conj() * b)
    };
    close(psi.expectation(&op).unwrap(), inner(&v, &av));
    close(
        psi.normalized_expectation(&op).unwrap(),
        inner(&v, &av) / inner(&v, &v),
    );
    close(
        psi.normalized_expectation_in(&op, &h).unwrap(),
        inner(&v, &av) / inner(&v, &v),
    );
    let normalized = psi.normalized().unwrap();
    close(normalized.inner_product(&normalized).unwrap(), c(1., 0.));
    let norm = psi.norm().unwrap();
    values(
        &vector(&normalized),
        &v.iter().map(|x| *x / norm).collect::<Vec<_>>(),
    );
    let acted = psi.act(&op).unwrap();
    close(
        acted
            .normalized()
            .unwrap()
            .inner_product(&acted.normalized().unwrap())
            .unwrap(),
        c(1., 0.),
    );
    assert!(
        qmb_core::State::<u8, Ket, Op, num_complex::Complex<T>>::zero(&h)
            .normalized()
            .is_err()
    );
}
#[test]
fn normalize_and_expectation_f32() {
    check::<f32>();
}
#[test]
fn normalize_and_expectation_f64() {
    check::<f64>();
}
