//! An independent three-level algebra using only the public qmb-core API.

use std::{fmt::Debug, num::NonZeroUsize, ops::Neg};

use qmb_core::{
    Act, Adjoint, Basis, Complex, DimensionOverflow, InnerProduct, LinearCombination,
    LocalStructure, OperatorBasis, Product, StateBasis, Term, Trace,
    num_traits::{FromPrimitive, Num, ToPrimitive},
    operator_basis_len, state_basis_len,
};

trait Real: Num + Copy + Neg<Output = Self> + Debug + FromPrimitive + ToPrimitive {
    const TOLERANCE: f64;
}
impl Real for f32 {
    const TOLERANCE: f64 = 2e-5;
}
impl Real for f64 {
    const TOLERANCE: f64 = 1e-12;
}

fn c<T: Real>(re: f64, im: f64) -> Complex<T> {
    Complex::new(T::from_f64(re).unwrap(), T::from_f64(im).unwrap())
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Axis {
    A,
    B,
    C,
}
const AXES: [Axis; 3] = [Axis::A, Axis::B, Axis::C];
impl Axis {
    fn index(self) -> usize {
        match self {
            Self::A => 0,
            Self::B => 1,
            Self::C => 2,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Ket {
    Unit(Axis),
    Sum,   // |0> + |1>, dependent on the canonical primitives
    Phase, // |0> + i|1>
}

const KETS: [Ket; 5] = [
    Ket::Unit(Axis::A),
    Ket::Unit(Axis::B),
    Ket::Unit(Axis::C),
    Ket::Sum,
    Ket::Phase,
];

impl Ket {
    fn coordinates<T: Real>(self) -> [Complex<T>; 3] {
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

impl<T: Real> InnerProduct<T> for Ket {
    fn inner_product(&self, rhs: &Self) -> Complex<T> {
        dot(self.coordinates(), rhs.coordinates())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Operator {
    Unit(Axis, Axis),
    Identity, // A dependent primitive, outside the matrix-unit basis.
}

impl<T: Real> Product<T> for Operator {
    fn product(&self, rhs: &Self) -> LinearCombination<Self, T> {
        match (*self, *rhs) {
            (Self::Identity, op) | (op, Self::Identity) => LinearCombination::primitive(op),
            (Self::Unit(i, j), Self::Unit(k, l)) if j == k => {
                LinearCombination::primitive(Self::Unit(i, l))
            }
            _ => LinearCombination::zero(),
        }
    }
}

impl<T: Real> Act<Operator, T> for Ket {
    fn act(&self, operator: &Operator) -> LinearCombination<Self, T> {
        match *operator {
            Operator::Identity => LinearCombination::primitive(*self),
            Operator::Unit(i, j) => {
                LinearCombination::term(self.coordinates::<T>()[j.index()], Self::Unit(i))
            }
        }
    }
}

impl<T: Real> Adjoint<T> for Operator {
    fn adjoint(&self) -> LinearCombination<Self, T> {
        LinearCombination::primitive(match *self {
            Self::Unit(i, j) => Self::Unit(j, i),
            Self::Identity => Self::Identity,
        })
    }
}

impl<T: Real> Trace<T> for Operator {
    fn trace(&self) -> Complex<T> {
        match *self {
            Self::Identity => c(3.0, 0.0),
            Self::Unit(i, j) => c(if i == j { 1.0 } else { 0.0 }, 0.0),
        }
    }
}

struct ThreeLevel;
impl<T: Real> LocalStructure<T> for ThreeLevel {
    type State = Ket;
    type Operator = Operator;
    const DIMENSION: NonZeroUsize = NonZeroUsize::new(3).unwrap();
    fn identity() -> LinearCombination<Operator, T> {
        LinearCombination::from_terms(AXES.map(|i| Term {
            coefficient: c(1.0, 0.0),
            primitive: Operator::Unit(i, i),
        }))
    }
}

struct Canonical;
impl<T: Real> Basis<Ket, T> for Canonical {
    fn expand(&self, primitive: &Ket) -> LinearCombination<Ket, T> {
        LinearCombination::from_terms(primitive.coordinates::<T>().into_iter().zip(AXES).map(
            |(coefficient, axis)| Term {
                coefficient,
                primitive: Ket::Unit(axis),
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
impl<T: Real> StateBasis<ThreeLevel, T> for Canonical {}

struct NonOrthogonal;
impl<T: Real> Basis<Ket, T> for NonOrthogonal {
    fn expand(&self, primitive: &Ket) -> LinearCombination<Ket, T> {
        let [a, b, d] = primitive.coordinates::<T>();
        LinearCombination::from_terms([
            Term {
                coefficient: a - b,
                primitive: Ket::Unit(Axis::A),
            },
            Term {
                coefficient: b,
                primitive: Ket::Sum,
            },
            Term {
                coefficient: d,
                primitive: Ket::Unit(Axis::C),
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
impl<T: Real> StateBasis<ThreeLevel, T> for NonOrthogonal {}

struct MatrixUnits;
impl<T: Real> Basis<Operator, T> for MatrixUnits {
    fn expand(&self, primitive: &Operator) -> LinearCombination<Operator, T> {
        match primitive {
            Operator::Identity => ThreeLevel::identity(),
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
impl<T: Real> OperatorBasis<ThreeLevel, T> for MatrixUnits {}

type Vector<T> = [Complex<T>; 3];
type Matrix<T> = [[Complex<T>; 3]; 3];

fn dot<T: Real>(x: Vector<T>, y: Vector<T>) -> Complex<T> {
    x.into_iter()
        .zip(y)
        .fold(c(0.0, 0.0), |sum, (a, b)| sum + a.conj() * b)
}

// Matrix reference computations are independent of the delta-based local rules.
fn matrix<T: Real>(op: Operator) -> Matrix<T> {
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

fn matmul<T: Real>(a: Matrix<T>, b: Matrix<T>) -> Matrix<T> {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..3).fold(c(0.0, 0.0), |sum, k| sum + a[i][k] * b[k][j]))
    })
}

fn matvec<T: Real>(a: Matrix<T>, x: Vector<T>) -> Vector<T> {
    a.map(|row| {
        row.into_iter()
            .zip(x)
            .fold(c(0.0, 0.0), |sum, (a, x)| sum + a * x)
    })
}

fn operator_value<T: Real>(sum: &LinearCombination<Operator, T>) -> Matrix<T> {
    let mut result = [[c(0.0, 0.0); 3]; 3];
    for term in sum {
        let m = matrix(term.primitive);
        for i in 0..3 {
            for j in 0..3 {
                result[i][j] = result[i][j] + term.coefficient * m[i][j];
            }
        }
    }
    result
}

fn state_value<T: Real>(sum: &LinearCombination<Ket, T>) -> Vector<T> {
    let mut result = [c(0.0, 0.0); 3];
    for term in sum {
        for (i, value) in term.primitive.coordinates::<T>().into_iter().enumerate() {
            result[i] = result[i] + term.coefficient * value;
        }
    }
    result
}

fn product<T: Real>(
    a: &LinearCombination<Operator, T>,
    b: &LinearCombination<Operator, T>,
) -> LinearCombination<Operator, T> {
    LinearCombination::from_terms(a.terms().flat_map(|left| {
        b.terms().flat_map(move |right| {
            left.primitive.product(&right.primitive).into_iter().map(
                move |term: Term<Operator, T>| Term {
                    coefficient: left.coefficient * right.coefficient * term.coefficient,
                    primitive: term.primitive,
                },
            )
        })
    }))
}

fn action<T: Real>(a: &LinearCombination<Operator, T>, state: &Ket) -> Vector<T> {
    state_value(&LinearCombination::from_terms(a.terms().flat_map(|t| {
        state.act(&t.primitive).scaled(t.coefficient).into_iter()
    })))
}

// LocalStructure alone must expose the state's action capability.
fn act_primitive<L, T>(state: &L::State, operator: &L::Operator) -> LinearCombination<L::State, T>
where
    L: LocalStructure<T>,
{
    state.act(operator)
}

fn inner<T: Real>(a: &LinearCombination<Ket, T>, b: &LinearCombination<Ket, T>) -> Complex<T> {
    let mut result = c(0.0, 0.0);
    for x in a {
        for y in b {
            result = result
                + x.coefficient.conj() * y.coefficient * x.primitive.inner_product(&y.primitive);
        }
    }
    result
}

fn assert_complex<T: Real>(actual: Complex<T>, expected: Complex<T>) {
    for (actual, expected) in [(actual.re, expected.re), (actual.im, expected.im)] {
        let actual = actual.to_f64().unwrap();
        let expected = expected.to_f64().unwrap();
        assert!(
            (actual - expected).abs() <= T::TOLERANCE * (1.0 + expected.abs()),
            "{actual} != {expected}"
        );
    }
}
fn assert_vector<T: Real>(actual: Vector<T>, expected: Vector<T>) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert_complex(a, b);
    }
}
fn assert_matrix<T: Real>(actual: Matrix<T>, expected: Matrix<T>) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert_vector(a, b);
    }
}

fn check_rules<T: Real>() {
    let operators: Vec<_> = <MatrixUnits as Basis<Operator, T>>::elements(&MatrixUnits)
        .chain([Operator::Identity])
        .collect();
    let identity = <ThreeLevel as LocalStructure<T>>::identity();
    assert_eq!(identity.len(), 3);
    assert_matrix(operator_value(&identity), matrix(Operator::Identity));
    for &a in &operators {
        let a_sum = LinearCombination::<_, T>::primitive(a);
        let a_matrix = matrix::<T>(a);
        let adjoint = <Operator as Adjoint<T>>::adjoint(&a);
        let expected_adjoint =
            std::array::from_fn(|i| std::array::from_fn(|j| a_matrix[j][i].conj()));
        assert_matrix(operator_value(&adjoint), expected_adjoint);
        assert_matrix(operator_value(&product(&identity, &a_sum)), a_matrix);
        assert_matrix(operator_value(&product(&a_sum, &identity)), a_matrix);
        let expected_trace = (0..3).fold(c(0.0, 0.0), |sum, i| sum + a_matrix[i][i]);
        assert_complex(<Operator as Trace<T>>::trace(&a), expected_trace);
        for state in KETS {
            let applied = act_primitive::<ThreeLevel, T>(&state, &a);
            assert_vector(state_value(&applied), matvec(a_matrix, state.coordinates()));
            assert_vector(action(&identity, &state), state.coordinates());
            for x in KETS {
                assert_complex(
                    dot(x.coordinates(), state_value(&applied)),
                    dot(action(&adjoint, &x), state.coordinates()),
                );
            }
        }
        for &b in &operators {
            let ab = <Operator as Product<T>>::product(&a, &b);
            assert_matrix(operator_value(&ab), matmul(a_matrix, matrix(b)));
            for state in KETS {
                assert_vector(
                    action(&ab, &state),
                    matvec(
                        a_matrix,
                        state_value(&<Ket as Act<Operator, T>>::act(&state, &b)),
                    ),
                );
            }
            for &d in &operators {
                let bc = <Operator as Product<T>>::product(&b, &d);
                let d_sum = LinearCombination::primitive(d);
                assert_matrix(
                    operator_value(&product(&ab, &d_sum)),
                    operator_value(&product(&a_sum, &bc)),
                );
            }
            // Cyclicity of the ordinary trace.
            let trace = |sum: &LinearCombination<Operator, T>| {
                sum.terms().fold(c::<T>(0.0, 0.0), |v, t| {
                    v + t.coefficient * <Operator as Trace<T>>::trace(&t.primitive)
                })
            };
            assert_complex(
                trace(&ab),
                trace(&<Operator as Product<T>>::product(&b, &a)),
            );
        }
    }
    let x = LinearCombination::primitive(Ket::Sum).scaled(c::<T>(0.3, -0.7))
        + LinearCombination::primitive(Ket::Phase).scaled(c(1.2, 0.4));
    let y = LinearCombination::primitive(Ket::Unit(Axis::B)).scaled(c::<T>(-0.2, 0.6));
    assert_complex(inner(&x, &y), dot(state_value(&x), state_value(&y)));
    assert_complex(inner(&x, &y), inner(&y, &x).conj());
    assert_complex(
        <Ket as InnerProduct<T>>::inner_product(&Ket::Sum, &Ket::Phase),
        c(1.0, 1.0),
    );
    let dependent = LinearCombination::primitive(Ket::Sum)
        + LinearCombination::primitive(Ket::Unit(Axis::A)).scaled(c::<T>(-1.0, 0.0))
        + LinearCombination::primitive(Ket::Unit(Axis::B)).scaled(c(-1.0, 0.0));
    assert_eq!(dependent.len(), 3);
    assert_vector(state_value(&dependent), [c(0.0, 0.0); 3]);
    assert_complex(inner(&dependent, &dependent), c(0.0, 0.0));
}

fn check_state_basis<T: Real, B: StateBasis<ThreeLevel, T>>(basis: &B) {
    let elements: Vec<_> = basis.elements().collect();
    assert_eq!(elements.len(), state_basis_len::<ThreeLevel, T>());
    for (i, p) in elements.iter().enumerate() {
        assert_eq!(basis.index(p), Some(i));
    }
    for p in KETS {
        let expanded = basis.expand(&p);
        assert_vector(state_value(&expanded), p.coordinates());
        for term in &expanded {
            assert!(basis.index(&term.primitive).is_some());
        }
    }
}

fn check_bases<T: Real>() {
    check_state_basis::<T, _>(&Canonical);
    check_state_basis::<T, _>(&NonOrthogonal);
    assert_eq!(
        <Canonical as Basis<Ket, T>>::index(&Canonical, &Ket::Sum),
        None
    );
    assert_eq!(
        <NonOrthogonal as Basis<Ket, T>>::index(&NonOrthogonal, &Ket::Unit(Axis::B)),
        None
    );
    let expanded = <NonOrthogonal as Basis<Ket, T>>::expand(&NonOrthogonal, &Ket::Unit(Axis::B));
    assert_complex(inner(&expanded, &expanded), c(1.0, 0.0));
    // Squared coefficient magnitudes alone give 2, while the actual norm² is 1.
    let coefficient_squares = expanded
        .terms()
        .fold(T::zero(), |v, t| v + t.coefficient.norm_sqr());
    assert_complex(Complex::new(coefficient_squares, T::zero()), c(2.0, 0.0));
    fn require_operator_basis<T: Real, B: OperatorBasis<ThreeLevel, T>>(_: &B) {}
    require_operator_basis::<T, _>(&MatrixUnits);
    let elements: Vec<_> = <MatrixUnits as Basis<Operator, T>>::elements(&MatrixUnits).collect();
    assert_eq!(
        elements.len(),
        operator_basis_len::<ThreeLevel, T>().unwrap()
    );
    for (i, p) in elements.iter().enumerate() {
        assert_eq!(
            <MatrixUnits as Basis<Operator, T>>::index(&MatrixUnits, p),
            Some(i)
        );
        assert_matrix(
            operator_value(&<MatrixUnits as Basis<Operator, T>>::expand(
                &MatrixUnits,
                p,
            )),
            matrix(*p),
        );
    }
    assert_eq!(
        <MatrixUnits as Basis<Operator, T>>::index(&MatrixUnits, &Operator::Identity),
        None
    );
    assert_matrix(
        operator_value(&<MatrixUnits as Basis<Operator, T>>::expand(
            &MatrixUnits,
            &Operator::Identity,
        )),
        matrix(Operator::Identity),
    );
}

#[test]
fn three_level_rules_and_bases_work_at_both_precisions() {
    check_rules::<f32>();
    check_rules::<f64>();
    check_bases::<f32>();
    check_bases::<f64>();
    assert_eq!(
        <ThreeLevel as LocalStructure<f32>>::DIMENSION,
        <ThreeLevel as LocalStructure<f64>>::DIMENSION
    );
    fn same_primitives<
        L: LocalStructure<f32>
            + LocalStructure<
                f64,
                State = <L as LocalStructure<f32>>::State,
                Operator = <L as LocalStructure<f32>>::Operator,
            >,
    >() {
    }
    same_primitives::<ThreeLevel>();
}

struct Dimension<const D: usize>;
impl<const D: usize> LocalStructure<f64> for Dimension<D> {
    type State = Ket;
    type Operator = Operator;
    const DIMENSION: NonZeroUsize = NonZeroUsize::new(D).unwrap();
    fn identity() -> LinearCombination<Operator, f64> {
        // Only the metadata arithmetic is exercised by this test fixture.
        unimplemented!("dimension-only test fixture")
    }
}

#[test]
fn basis_dimensions_handle_unit_large_and_overflow_cases() {
    assert_eq!(state_basis_len::<Dimension<1>, f64>(), 1);
    assert_eq!(operator_basis_len::<Dimension<1>, f64>(), Ok(1));
    assert_eq!(state_basis_len::<ThreeLevel, f64>(), 3);
    assert_eq!(operator_basis_len::<ThreeLevel, f64>(), Ok(9));
    const FITS: usize = 1 << (usize::BITS / 2 - 1);
    const OVERFLOWS: usize = 1 << (usize::BITS / 2);
    assert_eq!(
        operator_basis_len::<Dimension<FITS>, f64>(),
        Ok(FITS * FITS)
    );
    assert_eq!(
        operator_basis_len::<Dimension<OVERFLOWS>, f64>(),
        Err(DimensionOverflow)
    );
    assert_eq!(
        state_basis_len::<Dimension<{ usize::MAX }>, f64>(),
        usize::MAX
    );
    assert_eq!(
        operator_basis_len::<Dimension<{ usize::MAX }>, f64>(),
        Err(DimensionOverflow)
    );
    let error: &dyn core::error::Error = &DimensionOverflow;
    assert_eq!(
        error.to_string(),
        "operator-space dimension overflows usize"
    );
    assert!(error.source().is_none());
}
