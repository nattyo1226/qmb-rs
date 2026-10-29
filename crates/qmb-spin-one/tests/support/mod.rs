//! Independent dense reference calculations; no product rules or basis expansion.
#![allow(dead_code)]

use num_complex::Complex;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use qmb_core::{LinearCombination, LocalOperator, LocalState};
use qmb_spin_one::{SpinOneOperator as O, SpinOneState as S};

pub trait Real: Float + FromPrimitive + ToPrimitive + core::fmt::Debug {
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

pub const STATES: [S; 3] = [S::Plus, S::Zero, S::Minus];
pub const OPS: [O; 15] = [
    O::I,
    O::Sx,
    O::Sy,
    O::Sz,
    O::SPlus,
    O::SMinus,
    O::Dxy,
    O::Q0,
    O::Qxy,
    O::Qyz,
    O::Qzx,
    O::Q1Plus,
    O::Q1Minus,
    O::Q2Plus,
    O::Q2Minus,
];
pub const CART: [O; 9] = [
    O::I,
    O::Sx,
    O::Sy,
    O::Sz,
    O::Dxy,
    O::Q0,
    O::Qxy,
    O::Qyz,
    O::Qzx,
];
pub const LADDER: [O; 9] = [
    O::I,
    O::SPlus,
    O::SMinus,
    O::Sz,
    O::Q0,
    O::Q1Plus,
    O::Q1Minus,
    O::Q2Plus,
    O::Q2Minus,
];
pub const COMMON: [O; 3] = [O::I, O::Sz, O::Q0];
pub type Matrix = Vec<Complex<f64>>;
pub type LocalOp<T> = LocalOperator<O, Complex<T>>;
pub type LocalKet<T> = LocalState<S, O, Complex<T>>;

pub fn c<T: Real>(re: f64, im: f64) -> Complex<T> {
    Complex::new(T::from_f64(re).unwrap(), T::from_f64(im).unwrap())
}
pub fn to64<T: Real>(z: Complex<T>) -> Complex<f64> {
    Complex::new(z.re.to_f64().unwrap(), z.im.to_f64().unwrap())
}
pub fn close<T: Real>(a: Complex<f64>, b: Complex<f64>) {
    assert!(
        (a - b).norm() <= T::ATOL + T::RTOL * a.norm().max(b.norm()),
        "{a:?} != {b:?}"
    );
}
pub fn values<T: Real>(a: &[Complex<f64>], b: &[Complex<f64>]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        close::<T>(*a, *b);
    }
}
pub fn add(a: &[Complex<f64>], b: &[Complex<f64>]) -> Matrix {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).map(|(a, b)| a + b).collect()
}
pub fn scale(a: &[Complex<f64>], z: Complex<f64>) -> Matrix {
    a.iter().map(|a| a * z).collect()
}
pub fn product(a: &[Complex<f64>], b: &[Complex<f64>]) -> Matrix {
    let n = (a.len() as f64).sqrt() as usize;
    assert_eq!(a.len(), n * n);
    assert_eq!(a.len(), b.len());
    (0..n * n)
        .map(|ij| {
            (0..n)
                .map(|k| a[(ij / n) * n + k] * b[k * n + ij % n])
                .sum()
        })
        .collect()
}
pub fn adjoint(a: &[Complex<f64>]) -> Matrix {
    let n = (a.len() as f64).sqrt() as usize;
    (0..n * n)
        .map(|ij| a[(ij % n) * n + ij / n].conj())
        .collect()
}
pub fn trace(a: &[Complex<f64>]) -> Complex<f64> {
    let n = (a.len() as f64).sqrt() as usize;
    (0..n).map(|i| a[i * n + i]).sum()
}
pub fn apply(a: &[Complex<f64>], x: &[Complex<f64>]) -> Matrix {
    let n = x.len();
    assert_eq!(a.len(), n * n);
    a.chunks(n)
        .map(|row| row.iter().zip(x).map(|(a, x)| a * x).sum())
        .collect()
}
pub fn dot(x: &[Complex<f64>], y: &[Complex<f64>]) -> Complex<f64> {
    assert_eq!(x.len(), y.len());
    x.iter().zip(y).map(|(x, y)| x.conj() * y).sum()
}
pub fn kron(a: &[Complex<f64>], b: &[Complex<f64>]) -> Matrix {
    let n = (a.len() as f64).sqrt() as usize;
    let m = (b.len() as f64).sqrt() as usize;
    (0..n * m)
        .flat_map(|i| (0..n * m).map(move |j| a[(i / m) * n + j / m] * b[(i % m) * m + j % m]))
        .collect()
}

// Spins are written directly, independently of the implementation. Quadrupoles
// are defined from these matrices, not from production match arms or helpers.
pub fn matrix(op: O) -> Matrix {
    let z = Complex::new(0., 0.);
    let one = Complex::new(1., 0.);
    let h = Complex::new(1. / 2_f64.sqrt(), 0.);
    let ih = Complex::new(0., 1. / 2_f64.sqrt());
    let x = vec![z, h, z, h, z, h, z, h, z];
    let y = vec![z, -ih, z, ih, z, -ih, z, ih, z];
    let sz = vec![one, z, z, z, z, z, z, z, -one];
    let plus = vec![
        z,
        2_f64.sqrt().into(),
        z,
        z,
        z,
        2_f64.sqrt().into(),
        z,
        z,
        z,
    ];
    let minus = adjoint(&plus);
    let anti = |a: &[Complex<f64>], b: &[Complex<f64>]| add(&product(a, b), &product(b, a));
    match op {
        O::I => vec![one, z, z, z, one, z, z, z, one],
        O::Sx => x,
        O::Sy => y,
        O::Sz => sz,
        O::SPlus => plus,
        O::SMinus => minus,
        O::Dxy => add(&product(&x, &x), &scale(&product(&y, &y), (-1.).into())),
        O::Q0 => vec![one, z, z, z, (-2.).into(), z, z, z, one],
        O::Qxy => anti(&x, &y),
        O::Qyz => anti(&y, &sz),
        O::Qzx => anti(&sz, &x),
        O::Q1Plus => anti(&sz, &plus),
        O::Q1Minus => anti(&sz, &minus),
        O::Q2Plus => product(&plus, &plus),
        O::Q2Minus => product(&minus, &minus),
    }
}
pub fn state_index(s: S) -> usize {
    match s {
        S::Plus => 0,
        S::Zero => 1,
        S::Minus => 2,
    }
}
pub fn ket(s: S) -> Matrix {
    (0..3)
        .map(|i| Complex::new(f64::from(i == state_index(s)), 0.))
        .collect()
}
pub fn combination<T: Real>(lc: LinearCombination<O, Complex<T>>) -> Matrix {
    lc.into_iter().fold(vec![0.0.into(); 9], |sum, t| {
        add(&sum, &scale(&matrix(t.operand), to64(t.coefficient)))
    })
}
pub fn state_combination<T: Real>(lc: LinearCombination<S, Complex<T>>) -> Matrix {
    let mut result = vec![0.0.into(); 3];
    for term in lc {
        result[state_index(term.operand)] += to64(term.coefficient);
    }
    result
}
pub fn expression<T: Real>(lc: LinearCombination<O, Complex<T>>) -> LocalOp<T> {
    LocalOperator::from_terms(lc.map(LocalOperator::primitive))
}
pub fn local_value<T: Real>(op: &LocalOp<T>) -> Matrix {
    let mut result = vec![0.0.into(); 9];
    for (col, s) in STATES.into_iter().enumerate() {
        let psi = qmb_spin_one::SpinOne::<Complex<T>>::state(s).act(op);
        for (row, t) in STATES.into_iter().enumerate() {
            result[row * 3 + col] =
                to64(qmb_spin_one::SpinOne::<Complex<T>>::state(t).inner_product(&psi));
        }
    }
    result
}
