//! Test-only homogeneous two-level algebra and independent dense reference.
#![allow(dead_code)]
use num_complex::Complex;
use num_traits::{FromPrimitive, Num, ToPrimitive};
use qmb_core::*;
use std::{fmt::Debug, num::NonZeroUsize};

pub trait Real:
    Num + Copy + core::ops::Neg<Output = Self> + FromPrimitive + ToPrimitive + Debug + 'static
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
#[derive(Clone, Copy, Debug)]
pub struct Ket(pub u8);
#[derive(Clone, Copy, Debug)]
pub enum Op {
    I,
    P,
    U(u8, u8),
}
pub const A: Op = Op::U(0, 1);
pub const CREATE: Op = Op::U(1, 0);
pub const OPS: [Op; 6] = [Op::I, Op::P, Op::U(0, 0), A, CREATE, Op::U(1, 1)];
impl HomogeneousParity for Ket {
    fn parity(&self) -> Parity {
        if self.0 == 0 {
            Parity::Even
        } else {
            Parity::Odd
        }
    }
}
impl HomogeneousParity for Op {
    fn parity(&self) -> Parity {
        match self {
            Self::U(i, j) if i != j => Parity::Odd,
            _ => Parity::Even,
        }
    }
}
impl Identity for Op {
    fn identity() -> Self {
        Self::I
    }
}
impl ParityOperator for Op {
    fn parity_operator() -> Self {
        Self::P
    }
}
impl Dimension for Ket {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(2).unwrap())
    }
}
impl Dimension for Op {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(4).unwrap())
    }
}
impl<T: Real> Act<Complex<T>> for Ket {
    type Operator = Op;
    type Term = Self;
    fn act(&self, op: &Op) -> LinearCombination<Self, Complex<T>> {
        match *op {
            Op::I => LinearCombination::primitive(*self),
            Op::P => LinearCombination::term(c(if self.0 == 0 { 1. } else { -1. }, 0.), *self),
            Op::U(i, j) if self.0 == j => LinearCombination::primitive(Ket(i)),
            _ => LinearCombination::zero(),
        }
    }
}
impl<T: Real> Adjoint<Complex<T>> for Op {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, Complex<T>> {
        LinearCombination::primitive(match *self {
            Self::U(i, j) => Self::U(j, i),
            other => other,
        })
    }
}
impl<T: Real> Product<Complex<T>> for Op {
    type Term = Self;
    fn product(&self, rhs: &Self) -> LinearCombination<Self, Complex<T>> {
        use Op::*;
        match (*self, *rhs) {
            (I, x) | (x, I) => LinearCombination::primitive(x),
            (P, P) => LinearCombination::primitive(I),
            (P, U(i, j)) => LinearCombination::term(c(if i == 0 { 1. } else { -1. }, 0.), U(i, j)),
            (U(i, j), P) => LinearCombination::term(c(if j == 0 { 1. } else { -1. }, 0.), U(i, j)),
            (U(i, j), U(k, l)) if j == k => LinearCombination::primitive(U(i, l)),
            _ => LinearCombination::zero(),
        }
    }
}
impl<T: Real> Trace<Complex<T>> for Op {
    fn trace(&self) -> Complex<T> {
        c(
            match self {
                Self::I => 2.,
                Self::U(i, j) if i == j => 1.,
                _ => 0.,
            },
            0.,
        )
    }
}
impl<T: Real> InnerProduct<Complex<T>> for Ket {
    fn inner_product(&self, rhs: &Self) -> Complex<T> {
        c(if self.0 == rhs.0 { 1. } else { 0. }, 0.)
    }
}
pub struct States;
impl<T: Real> Basis<Ket, Complex<T>> for States {
    fn elements(&self) -> impl Iterator<Item = Ket> + '_ {
        [Ket(0), Ket(1)].into_iter()
    }
    fn index(&self, p: &Ket) -> Option<usize> {
        Some(p.0 as usize)
    }
    fn expand(&self, p: &Ket) -> LinearCombination<Ket, Complex<T>> {
        LinearCombination::primitive(*p)
    }
}
pub struct Units;
impl<T: Real> Basis<Op, Complex<T>> for Units {
    fn elements(&self) -> impl Iterator<Item = Op> + '_ {
        [Op::U(0, 0), A, CREATE, Op::U(1, 1)].into_iter()
    }
    fn index(&self, p: &Op) -> Option<usize> {
        if let Op::U(i, j) = *p {
            Some((2 * i + j) as usize)
        } else {
            None
        }
    }
    fn expand(&self, p: &Op) -> LinearCombination<Op, Complex<T>> {
        match p {
            Op::I => {
                LinearCombination::primitive(Op::U(0, 0))
                    + LinearCombination::primitive(Op::U(1, 1))
            }
            Op::P => {
                LinearCombination::primitive(Op::U(0, 0))
                    + LinearCombination::term(c(-1., 0.), Op::U(1, 1))
            }
            _ => LinearCombination::primitive(*p),
        }
    }
}
pub struct Parities;
impl<T: Real> Basis<Op, Complex<T>> for Parities {
    fn elements(&self) -> impl Iterator<Item = Op> + '_ {
        [Op::I, Op::P, A, CREATE].into_iter()
    }
    fn index(&self, p: &Op) -> Option<usize> {
        match *p {
            Op::I => Some(0),
            Op::P => Some(1),
            Op::U(0, 1) => Some(2),
            Op::U(1, 0) => Some(3),
            _ => None,
        }
    }
    fn expand(&self, p: &Op) -> LinearCombination<Op, Complex<T>> {
        match *p {
            Op::U(i, j) if i == j => {
                LinearCombination::term(c(0.5, 0.), Op::I)
                    + LinearCombination::term(c(if i == 0 { 0.5 } else { -0.5 }, 0.), Op::P)
            }
            _ => LinearCombination::primitive(*p),
        }
    }
}
pub type Psi<T> = State<u8, Ket, Op, Complex<T>>;
pub type OperatorOf<T> = Operator<u8, Op, Complex<T>>;
pub fn placed<T: Real>(op: Op, factor: &HilbertSpace<u8>) -> OperatorOf<T> {
    LocalOperator::primitive(op).at(factor).unwrap()
}
pub fn ket<T: Real>(bits: usize, factors: &[HilbertSpace<u8>], kind: TensorKind) -> Psi<T> {
    State::tensor_in(
        &HilbertSpace::tensor_all_with(kind, factors.iter().cloned()).unwrap(),
        factors.iter().enumerate().map(|(i, f)| {
            LocalState::primitive(Ket(((bits >> (factors.len() - 1 - i)) & 1) as u8))
                .at(f)
                .unwrap()
        }),
    )
    .unwrap()
}
pub fn vector<T: Real>(psi: &Psi<T>) -> Vec<Complex<T>> {
    let coords = psi.coordinates(&States).unwrap();
    let mut result = vec![c(0., 0.); coords.dim().unwrap()];
    for entry in coords.flat_entries() {
        let (i, v) = entry.unwrap();
        result[i] = *v;
    }
    result
}
pub fn close<T: Real>(a: Complex<T>, b: Complex<T>) {
    assert!(
        (a - b).norm_sqr().to_f64().unwrap() < T::TOL * T::TOL,
        "{a:?} != {b:?}"
    );
}
pub fn values<T: Real>(a: &[Complex<T>], b: &[Complex<T>]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        close(*x, *y);
    }
}
// All reference operations below are plain dense arithmetic, independent of
// primitive algebra rules, evaluator, basis expansion and export.
pub fn local<T: Real>(op: Op) -> Vec<Complex<T>> {
    let mut m = vec![c(0., 0.); 4];
    match op {
        Op::I => {
            m[0] = c(1., 0.);
            m[3] = c(1., 0.);
        }
        Op::P => {
            m[0] = c(1., 0.);
            m[3] = c(-1., 0.);
        }
        Op::U(i, j) => m[(2 * i + j) as usize] = c(1., 0.),
    };
    m
}
pub fn kron<T: Real>(a: &[Complex<T>], b: &[Complex<T>]) -> Vec<Complex<T>> {
    let n = (a.len() as f64).sqrt() as usize;
    let m = (b.len() as f64).sqrt() as usize;
    let mut out = vec![c(0., 0.); n * n * m * m];
    for i in 0..n {
        for j in 0..n {
            for k in 0..m {
                for l in 0..m {
                    out[(i * m + k) * (n * m) + j * m + l] = a[i * n + j] * b[k * m + l];
                }
            }
        }
    }
    out
}
pub fn embedding<T: Real>(op: Op, at: usize, n: usize, kind: TensorKind) -> Vec<Complex<T>> {
    let odd = matches!(op,Op::U(i,j) if i!=j);
    let mut out = vec![c(1., 0.)];
    for i in 0..n {
        out = kron(
            &out,
            &local(if i == at {
                op
            } else if i < at && odd && kind == TensorKind::Graded {
                Op::P
            } else {
                Op::I
            }),
        );
    }
    out
}
pub fn product<T: Real>(a: &[Complex<T>], b: &[Complex<T>]) -> Vec<Complex<T>> {
    let n = (a.len() as f64).sqrt() as usize;
    let mut out = vec![c(0., 0.); n * n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                out[i * n + j] = out[i * n + j] + a[i * n + k] * b[k * n + j];
            }
        }
    }
    out
}
pub fn adjoint<T: Real>(a: &[Complex<T>]) -> Vec<Complex<T>> {
    let n = (a.len() as f64).sqrt() as usize;
    (0..n * n).map(|k| a[(k % n) * n + k / n].conj()).collect()
}
pub fn trace<T: Real>(a: &[Complex<T>]) -> Complex<T> {
    let n = (a.len() as f64).sqrt() as usize;
    (0..n).fold(c(0., 0.), |sum, i| sum + a[i * n + i])
}
pub fn act<T: Real>(a: &[Complex<T>], v: &[Complex<T>]) -> Vec<Complex<T>> {
    let n = v.len();
    (0..n)
        .map(|i| (0..n).fold(c(0., 0.), |sum, j| sum + a[i * n + j] * v[j]))
        .collect()
}
