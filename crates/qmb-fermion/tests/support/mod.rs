#![allow(dead_code)]

use num_complex::Complex;
use num_traits::{Float, FromPrimitive};
use qmb_core::{LinearCombination, LocalOperator};
use qmb_fermion::{Fermion, FermionOperator as O, FermionState as S, LadderBasis};

pub const OPS: [O; 8] = [
    O::I,
    O::Annihilate,
    O::Create,
    O::Number,
    O::Parity,
    O::MajoranaX,
    O::MajoranaY,
    O::MajoranaZ,
];
pub const STATES: [S; 2] = [S::Empty, S::Occupied];
pub type Matrix = Vec<Complex<f64>>;
pub trait Real: Float + FromPrimitive + core::fmt::Debug {
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
pub fn to64<T: Real>(x: Complex<T>) -> Complex<f64> {
    Complex::new(x.re.to_f64().unwrap(), x.im.to_f64().unwrap())
}
pub fn close<T: Real>(actual: Complex<f64>, expected: Complex<f64>) {
    assert!(
        (actual - expected).norm() <= T::TOL * (1. + expected.norm()),
        "{actual:?} != {expected:?}"
    );
}
pub fn values<T: Real>(actual: &[Complex<f64>], expected: &[Complex<f64>]) {
    assert_eq!(actual.len(), expected.len());
    for (a, b) in actual.iter().zip(expected) {
        close::<T>(*a, *b);
    }
}

// Literal reference matrices; independent of production rules, evaluation,
// basis conversion, and export. State order is Empty, Occupied.
pub fn matrix(op: O) -> Matrix {
    let i = Complex::new(0., 1.);
    match op {
        O::I => vec![1.0.into(), 0.0.into(), 0.0.into(), 1.0.into()],
        O::Annihilate => vec![0.0.into(), 1.0.into(), 0.0.into(), 0.0.into()],
        O::Create => vec![0.0.into(), 0.0.into(), 1.0.into(), 0.0.into()],
        O::Number => vec![0.0.into(), 0.0.into(), 0.0.into(), 1.0.into()],
        O::Parity => vec![1.0.into(), 0.0.into(), 0.0.into(), (-1.0).into()],
        O::MajoranaX => vec![0.0.into(), 0.5.into(), 0.5.into(), 0.0.into()],
        O::MajoranaY => vec![0.0.into(), -0.5 * i, 0.5 * i, 0.0.into()],
        O::MajoranaZ => vec![0.5.into(), 0.0.into(), 0.0.into(), (-0.5).into()],
    }
}
pub fn add(a: &[Complex<f64>], b: &[Complex<f64>]) -> Matrix {
    a.iter().zip(b).map(|(a, b)| a + b).collect()
}
pub fn scale(a: &[Complex<f64>], c: Complex<f64>) -> Matrix {
    a.iter().map(|a| a * c).collect()
}
pub fn product(a: &[Complex<f64>], b: &[Complex<f64>]) -> Matrix {
    let n = (a.len() as f64).sqrt() as usize;
    (0..n * n)
        .map(|k| (0..n).map(|j| a[(k / n) * n + j] * b[j * n + k % n]).sum())
        .collect()
}
pub fn adjoint(a: &[Complex<f64>]) -> Matrix {
    let n = (a.len() as f64).sqrt() as usize;
    (0..n * n).map(|k| a[(k % n) * n + k / n].conj()).collect()
}
pub fn trace(a: &[Complex<f64>]) -> Complex<f64> {
    let n = (a.len() as f64).sqrt() as usize;
    (0..n).map(|k| a[k * n + k]).sum()
}
pub fn combination<T: Real>(terms: LinearCombination<O, Complex<T>>) -> Matrix {
    terms.into_iter().fold(vec![0.0.into(); 4], |m, t| {
        add(&m, &scale(&matrix(t.operand), to64(t.coefficient)))
    })
}
pub fn expression<T: Real>(
    terms: LinearCombination<O, Complex<T>>,
) -> LocalOperator<O, Complex<T>> {
    LocalOperator::from_terms(terms.map(Fermion::<Complex<T>>::operator))
}
pub fn local_value<T: Real>(op: &LocalOperator<O, Complex<T>>) -> Matrix {
    let coords = op.coordinates(&LadderBasis).unwrap();
    let basis = [O::I, O::Annihilate, O::Create, O::Number];
    coords.flat_entries().fold(vec![0.0.into(); 4], |m, t| {
        let (j, c) = t.unwrap();
        add(&m, &scale(&matrix(basis[j]), to64(*c)))
    })
}

// Independent Fock-space action: bit order follows the input mode order,
// leftmost factor most significant. The sign counts earlier occupied modes.
pub fn fock_action(bits: usize, mode: usize, modes: usize, create: bool) -> Option<(usize, f64)> {
    let mask = 1 << (modes - 1 - mode);
    if (bits & mask != 0) == create {
        return None;
    }
    let sign = if (bits >> (modes - mode)).count_ones() % 2 == 0 {
        1.
    } else {
        -1.
    };
    Some((bits ^ mask, sign))
}
pub fn fock_matrix(mode: usize, modes: usize, create: bool) -> Matrix {
    let n = 1 << modes;
    let mut result = vec![0.0.into(); n * n];
    for col in 0..n {
        if let Some((row, sign)) = fock_action(col, mode, modes, create) {
            result[row * n + col] = sign.into();
        }
    }
    result
}

pub fn exported_matrix<T: Real, I: Eq>(op: &qmb_core::Operator<I, O, Complex<T>>) -> Matrix {
    use qmb_export::{Export, Exported, MatrixArgs};
    use qmb_fermion::OccupationBasis;
    let dense: Exported<Vec<Complex<T>>, _> = op
        .export(MatrixArgs::<S, _>::new(&OccupationBasis))
        .unwrap();
    let array: Exported<qmb_export::ndarray::Array2<Complex<T>>, _> = op
        .export(MatrixArgs::<S, _>::new(&OccupationBasis))
        .unwrap();
    let sparse: Exported<qmb_export::sprs::CsMat<Complex<T>>, _> = op
        .export(MatrixArgs::<S, _>::new(&OccupationBasis))
        .unwrap();
    let result: Matrix = dense.data().iter().map(|x| to64(*x)).collect();
    values::<T>(
        &array.data().iter().map(|x| to64(*x)).collect::<Vec<_>>(),
        &result,
    );
    values::<T>(
        &sparse
            .data()
            .to_dense()
            .iter()
            .map(|x| to64(*x))
            .collect::<Vec<_>>(),
        &result,
    );
    let space = op.space().unwrap();
    assert!(dense.metadata().space().is_compatible_with(space));
    assert!(array.metadata().space().is_compatible_with(space));
    assert!(sparse.metadata().space().is_compatible_with(space));
    result
}
