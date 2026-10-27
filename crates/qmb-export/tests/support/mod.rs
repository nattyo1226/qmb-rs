#![allow(dead_code)]
use num_complex::Complex;
use num_traits::{FromPrimitive, Num, NumAssign, ToPrimitive};
use qmb_spin_half::SpinHalfOperator as Op;
use std::{fmt::Debug, ops::Neg};

pub trait Real:
    Num
    + NumAssign
    + Default
    + Send
    + Sync
    + Copy
    + Neg<Output = Self>
    + FromPrimitive
    + ToPrimitive
    + Debug
{
    const TOL: f64;
}
impl Real for f32 {
    const TOL: f64 = 2e-5;
}
impl Real for f64 {
    const TOL: f64 = 1e-12;
}
pub fn c<T: Real>(re: f64, im: f64) -> Complex<T> {
    Complex::new(T::from_f64(re).unwrap(), T::from_f64(im).unwrap())
}
pub fn close<T: Real>(a: Complex<T>, b: Complex<T>) {
    let magnitude = |z: Complex<T>| z.re.to_f64().unwrap().hypot(z.im.to_f64().unwrap());
    assert!(
        magnitude(a - b) <= T::TOL + T::TOL * magnitude(a).max(magnitude(b)),
        "{a:?} != {b:?}"
    );
}
pub fn values<T: Real>(a: &[Complex<T>], b: &[Complex<T>]) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) {
        close(a, b);
    }
}
pub fn local<T: Real>(op: Op) -> Vec<Complex<T>> {
    use Op::*;
    match op {
        I => vec![c(1., 0.), c(0., 0.), c(0., 0.), c(1., 0.)],
        Sx => vec![c(0., 0.), c(0.5, 0.), c(0.5, 0.), c(0., 0.)],
        Sy => vec![c(0., 0.), c(0., -0.5), c(0., 0.5), c(0., 0.)],
        Sz => vec![c(0.5, 0.), c(0., 0.), c(0., 0.), c(-0.5, 0.)],
        SPlus => vec![c(0., 0.), c(1., 0.), c(0., 0.), c(0., 0.)],
        SMinus => vec![c(0., 0.), c(0., 0.), c(1., 0.), c(0., 0.)],
    }
}
pub fn product<T: Real>(a: &[Complex<T>], b: &[Complex<T>], d: usize) -> Vec<Complex<T>> {
    (0..d * d)
        .map(|ij| (0..d).fold(c(0., 0.), |s, k| s + a[ij / d * d + k] * b[k * d + ij % d]))
        .collect()
}
pub fn adjoint<T: Real>(a: &[Complex<T>], d: usize) -> Vec<Complex<T>> {
    (0..d * d).map(|ij| a[ij % d * d + ij / d].conj()).collect()
}
pub fn act<T: Real>(a: &[Complex<T>], x: &[Complex<T>]) -> Vec<Complex<T>> {
    let d = x.len();
    (0..d)
        .map(|i| {
            x.iter()
                .enumerate()
                .fold(c(0., 0.), |s, (j, &value)| s + a[i * d + j] * value)
        })
        .collect()
}
pub fn dot<T: Real>(a: &[Complex<T>], b: &[Complex<T>]) -> Complex<T> {
    a.iter()
        .zip(b)
        .fold(c(0., 0.), |s, (&a, &b)| s + a.conj() * b)
}
pub fn kron<T: Real>(a: &[Complex<T>], b: &[Complex<T>]) -> Vec<Complex<T>> {
    (0..16)
        .map(|ij| a[(ij / 4 / 2) * 2 + (ij % 4 / 2)] * b[(ij / 4 % 2) * 2 + ij % 2])
        .collect()
}
