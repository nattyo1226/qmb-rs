//! Shared three-level algebra and independent matrix reference for tests.
#![allow(dead_code, unused_imports)]

use std::{fmt::Debug, num::NonZeroUsize, ops::Neg};

use num_complex::Complex;
use num_traits::{FromPrimitive, Num, ToPrimitive};
use qmb_core::{
    Act, Adjoint, Basis, Dimension, DimensionOverflow, Identity, InnerProduct, LinearCombination,
    Product, Term, Trace, operator_basis_len, state_basis_len,
};

pub(super) trait Real:
    Num + Copy + Neg<Output = Self> + Debug + FromPrimitive + ToPrimitive
{
    const TOLERANCE: f64;
}
impl Real for f32 {
    const TOLERANCE: f64 = 2e-5;
}
impl Real for f64 {
    const TOLERANCE: f64 = 1e-12;
}

pub(super) fn c<T: Real>(re: f64, im: f64) -> Complex<T> {
    Complex::new(T::from_f64(re).unwrap(), T::from_f64(im).unwrap())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Axis {
    A,
    B,
    C,
}
pub(super) const AXES: [Axis; 3] = [Axis::A, Axis::B, Axis::C];
impl Axis {
    pub(super) fn index(self) -> usize {
        match self {
            Self::A => 0,
            Self::B => 1,
            Self::C => 2,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Ket {
    Unit(Axis),
    Sum,   // |0> + |1>, dependent on the canonical primitives
    Phase, // |0> + i|1>
}

pub(super) const KETS: [Ket; 5] = [
    Ket::Unit(Axis::A),
    Ket::Unit(Axis::B),
    Ket::Unit(Axis::C),
    Ket::Sum,
    Ket::Phase,
];

impl Ket {
    pub(super) fn coordinates<T: Real>(self) -> [Complex<T>; 3] {
        let mut result = [c(0.0, 0.0); 3];
        match self {
            Self::Unit(axis) => result[axis.index()] = c(1.0, 0.0),
            Self::Sum => {
                result[0] = c(1.0, 0.0);
                result[1] = c(1.0, 0.0);
            }
            Self::Phase => {
                result[0] = c(1.0, 0.0);
                result[1] = c(0.0, 1.0);
            }
        }
        result
    }
}

impl<T: Real> InnerProduct<Complex<T>> for Ket {
    fn inner_product(&self, rhs: &Self) -> Complex<T> {
        dot(self.coordinates(), rhs.coordinates())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Operator {
    Unit(Axis, Axis),
    Identity, // A dependent primitive, outside the matrix-unit basis.
}

impl<T: Real> Product<Complex<T>> for Operator {
    type Term = Self;
    fn product(&self, rhs: &Self) -> LinearCombination<Self, Complex<T>> {
        match (*self, *rhs) {
            (Self::Identity, op) | (op, Self::Identity) => LinearCombination::primitive(op),
            (Self::Unit(i, j), Self::Unit(k, l)) if j == k => {
                LinearCombination::primitive(Self::Unit(i, l))
            }
            _ => LinearCombination::zero(),
        }
    }
}

impl<T: Real> Act<Complex<T>> for Ket {
    type Operator = Operator;
    type Term = Self;
    fn act(&self, operator: &Operator) -> LinearCombination<Self, Complex<T>> {
        match *operator {
            Operator::Identity => LinearCombination::primitive(*self),
            Operator::Unit(i, j) => {
                LinearCombination::term(self.coordinates::<T>()[j.index()], Self::Unit(i))
            }
        }
    }
}

impl<T: Real> Adjoint<Complex<T>> for Operator {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, Complex<T>> {
        LinearCombination::primitive(match *self {
            Self::Unit(i, j) => Self::Unit(j, i),
            Self::Identity => Self::Identity,
        })
    }
}

impl<T: Real> Trace<Complex<T>> for Operator {
    fn trace(&self) -> Complex<T> {
        match *self {
            Self::Identity => c(3.0, 0.0),
            Self::Unit(i, j) => c(if i == j { 1.0 } else { 0.0 }, 0.0),
        }
    }
}

impl Identity for Operator {
    fn identity() -> Self {
        Self::Identity
    }
}

pub(super) struct Canonical;
impl<T: Real> Basis<Ket, Complex<T>> for Canonical {
    fn expand(&self, primitive: &Ket) -> LinearCombination<Ket, Complex<T>> {
        LinearCombination::from_terms(primitive.coordinates::<T>().into_iter().zip(AXES).map(
            |(coefficient, axis)| Term {
                coefficient,
                operand: Ket::Unit(axis),
            },
        ))
    }
    fn index(&self, primitive: &Ket) -> Option<usize> {
        match primitive {
            Ket::Unit(i) => Some(i.index()),
            _ => None,
        }
    }
    fn elements(&self) -> impl Iterator<Item = Ket> + '_ {
        AXES.into_iter().map(Ket::Unit)
    }
}

pub(super) struct NonOrthogonal;
impl<T: Real> Basis<Ket, Complex<T>> for NonOrthogonal {
    fn expand(&self, primitive: &Ket) -> LinearCombination<Ket, Complex<T>> {
        let [a, b, d] = primitive.coordinates::<T>();
        LinearCombination::from_terms([
            Term {
                coefficient: a - b,
                operand: Ket::Unit(Axis::A),
            },
            Term {
                coefficient: b,
                operand: Ket::Sum,
            },
            Term {
                coefficient: d,
                operand: Ket::Unit(Axis::C),
            },
        ])
    }
    fn index(&self, primitive: &Ket) -> Option<usize> {
        match primitive {
            Ket::Unit(Axis::A) => Some(0),
            Ket::Sum => Some(1),
            Ket::Unit(Axis::C) => Some(2),
            _ => None,
        }
    }
    fn elements(&self) -> impl Iterator<Item = Ket> + '_ {
        [Ket::Unit(Axis::A), Ket::Sum, Ket::Unit(Axis::C)].into_iter()
    }
}

pub(super) struct MatrixUnits;
impl<T: Real> Basis<Operator, Complex<T>> for MatrixUnits {
    fn expand(&self, primitive: &Operator) -> LinearCombination<Operator, Complex<T>> {
        match primitive {
            Operator::Identity => LinearCombination::from_terms(AXES.map(|i| Term {
                coefficient: c(1.0, 0.0),
                operand: Operator::Unit(i, i),
            })),
            _ => LinearCombination::primitive(*primitive),
        }
    }
    fn index(&self, primitive: &Operator) -> Option<usize> {
        match primitive {
            Operator::Unit(i, j) => Some(i.index() * 3 + j.index()),
            _ => None,
        }
    }
    fn elements(&self) -> impl Iterator<Item = Operator> + '_ {
        AXES.into_iter()
            .flat_map(|i| AXES.into_iter().map(move |j| Operator::Unit(i, j)))
    }
}

pub(super) type Vector<T> = [Complex<T>; 3];
pub(super) type Matrix<T> = [[Complex<T>; 3]; 3];

pub(super) fn dot<T: Real>(x: Vector<T>, y: Vector<T>) -> Complex<T> {
    x.into_iter()
        .zip(y)
        .fold(c(0.0, 0.0), |sum, (a, b)| sum + a.conj() * b)
}

// Matrix reference computations are independent of the delta-based local rules.
pub(super) fn matrix<T: Real>(op: Operator) -> Matrix<T> {
    let mut result = [[c(0.0, 0.0); 3]; 3];
    match op {
        Operator::Unit(i, j) => result[i.index()][j.index()] = c(1.0, 0.0),
        Operator::Identity => {
            for (i, row) in result.iter_mut().enumerate() {
                row[i] = c(1.0, 0.0);
            }
        }
    }
    result
}

pub(super) fn matmul<T: Real>(a: Matrix<T>, b: Matrix<T>) -> Matrix<T> {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..3).fold(c(0.0, 0.0), |sum, k| sum + a[i][k] * b[k][j]))
    })
}

pub(super) fn matvec<T: Real>(a: Matrix<T>, x: Vector<T>) -> Vector<T> {
    a.map(|row| {
        row.into_iter()
            .zip(x)
            .fold(c(0.0, 0.0), |sum, (a, x)| sum + a * x)
    })
}

pub(super) fn operator_value<T: Real>(sum: &LinearCombination<Operator, Complex<T>>) -> Matrix<T> {
    let mut result = [[c(0.0, 0.0); 3]; 3];
    for term in sum {
        let m = matrix(term.operand);
        for i in 0..3 {
            for j in 0..3 {
                result[i][j] = result[i][j] + term.coefficient * m[i][j];
            }
        }
    }
    result
}

pub(super) fn state_value<T: Real>(sum: &LinearCombination<Ket, Complex<T>>) -> Vector<T> {
    let mut result = [c(0.0, 0.0); 3];
    for term in sum {
        for (i, value) in term.operand.coordinates::<T>().into_iter().enumerate() {
            result[i] = result[i] + term.coefficient * value;
        }
    }
    result
}

pub(super) fn product<T: Real>(
    a: &LinearCombination<Operator, Complex<T>>,
    b: &LinearCombination<Operator, Complex<T>>,
) -> LinearCombination<Operator, Complex<T>> {
    LinearCombination::from_terms(a.terms().flat_map(|left| {
        b.terms().flat_map(move |right| {
            left.operand.product(&right.operand).into_iter().map(
                move |term: Term<Operator, Complex<T>>| Term {
                    coefficient: left.coefficient * right.coefficient * term.coefficient,
                    operand: term.operand,
                },
            )
        })
    }))
}

pub(super) fn action<T: Real>(
    a: &LinearCombination<Operator, Complex<T>>,
    state: &Ket,
) -> Vector<T> {
    state_value(&LinearCombination::from_terms(a.terms().flat_map(|t| {
        state.act(&t.operand).scaled(t.coefficient).into_iter()
    })))
}

// Evaluation requests action explicitly, independently of metadata.
pub(super) fn act_primitive<S, O, T>(state: &S, operator: &O) -> LinearCombination<S, Complex<T>>
where
    S: Act<Complex<T>, Operator = O, Term = S>,
{
    state.act(operator)
}

pub(super) fn inner<T: Real>(
    a: &LinearCombination<Ket, Complex<T>>,
    b: &LinearCombination<Ket, Complex<T>>,
) -> Complex<T> {
    let mut result = c(0.0, 0.0);
    for x in a {
        for y in b {
            result =
                result + x.coefficient.conj() * y.coefficient * x.operand.inner_product(&y.operand);
        }
    }
    result
}

pub(super) fn assert_complex<T: Real>(actual: Complex<T>, expected: Complex<T>) {
    for (actual, expected) in [(actual.re, expected.re), (actual.im, expected.im)] {
        let actual = actual.to_f64().unwrap();
        let expected = expected.to_f64().unwrap();
        assert!(
            (actual - expected).abs() <= T::TOLERANCE * (1.0 + expected.abs()),
            "{actual} != {expected}"
        );
    }
}
pub(super) fn assert_vector<T: Real>(actual: Vector<T>, expected: Vector<T>) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert_complex(a, b);
    }
}
pub(super) fn assert_matrix<T: Real>(actual: Matrix<T>, expected: Matrix<T>) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert_vector(a, b);
    }
}

impl Dimension for Ket {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(3).unwrap())
    }
}
impl Dimension for Operator {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(9).unwrap())
    }
}

impl qmb_core::HomogeneousParity for Ket {
    fn parity(&self) -> qmb_core::Parity {
        qmb_core::Parity::Even
    }
}
impl qmb_core::HomogeneousParity for Operator {
    fn parity(&self) -> qmb_core::Parity {
        qmb_core::Parity::Even
    }
}
impl qmb_core::ParityOperator for Operator {
    fn parity_operator() -> Self {
        Self::identity()
    }
}
