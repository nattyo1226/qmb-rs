//! Public contracts checked against independent two-dimensional linear algebra.

use std::{fmt::Debug, hash::Hash, ops::Neg};

use qmb_core::{
    ActOn, Adjoint, Basis, Complex, InnerProduct, LinearCombination as Lc, LocalStructure,
    OperatorBasis, Product, StateBasis, Trace,
    num_traits::{FromPrimitive, Num, ToPrimitive},
    operator_basis_len, state_basis_len,
};
use qmb_spin_half::{
    CartesianBasis, LadderBasis, SpinHalf, SpinHalfOperator as Op, SpinHalfState as State, SzBasis,
};

use Op::{I, SMinus, SPlus, Sx, Sy, Sz};
use State::{Down, Up};

const OPS: [Op; 6] = [I, Sx, Sy, Sz, SPlus, SMinus];
const STATES: [State; 2] = [Up, Down];

trait Real: Num + Copy + Neg<Output = Self> + Debug + FromPrimitive + ToPrimitive {
    const ATOL: f64;
    const RTOL: f64;
}
impl Real for f32 {
    const ATOL: f64 = 2e-6;
    const RTOL: f64 = 2e-6;
}
impl Real for f64 {
    const ATOL: f64 = 1e-12;
    const RTOL: f64 = 1e-12;
}

type Vector = [Complex<f64>; 2];
type Matrix = [Vector; 2];

fn c<T: Real>(re: f64, im: f64) -> Complex<T> {
    Complex::new(T::from_f64(re).unwrap(), T::from_f64(im).unwrap())
}
fn complex64<T: Real>(value: Complex<T>) -> Complex<f64> {
    Complex::new(value.re.to_f64().unwrap(), value.im.to_f64().unwrap())
}
fn vector(state: State) -> Vector {
    match state {
        Up => [c(1.0, 0.0), c(0.0, 0.0)],
        Down => [c(0.0, 0.0), c(1.0, 0.0)],
    }
}

// Literal reference matrices use no production rules or basis expansion.
fn matrix(op: Op) -> Matrix {
    match op {
        I => [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]],
        Sx => [[c(0.0, 0.0), c(0.5, 0.0)], [c(0.5, 0.0), c(0.0, 0.0)]],
        Sy => [[c(0.0, 0.0), c(0.0, -0.5)], [c(0.0, 0.5), c(0.0, 0.0)]],
        Sz => [[c(0.5, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(-0.5, 0.0)]],
        SPlus => [[c(0.0, 0.0), c(1.0, 0.0)], [c(0.0, 0.0), c(0.0, 0.0)]],
        SMinus => [[c(0.0, 0.0), c(0.0, 0.0)], [c(1.0, 0.0), c(0.0, 0.0)]],
    }
}
fn add(a: Matrix, b: Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] + b[i][j]))
}
fn scale(a: Matrix, value: Complex<f64>) -> Matrix {
    a.map(|row| row.map(|x| value * x))
}
fn multiply(a: Matrix, b: Matrix) -> Matrix {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..2).fold(c(0.0, 0.0), |v, k| v + a[i][k] * b[k][j]))
    })
}
fn apply(a: Matrix, x: Vector) -> Vector {
    a.map(|row| {
        row.into_iter()
            .zip(x)
            .fold(c(0.0, 0.0), |sum, (a, x)| sum + a * x)
    })
}
fn dagger(a: Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i].conj()))
}
fn dot(x: Vector, y: Vector) -> Complex<f64> {
    x.into_iter()
        .zip(y)
        .fold(c(0.0, 0.0), |sum, (x, y)| sum + x.conj() * y)
}
fn operator_value<T: Real>(sum: &Lc<Op, T>) -> Matrix {
    sum.terms().fold([[c(0.0, 0.0); 2]; 2], |v, t| {
        add(v, scale(matrix(t.primitive), complex64(t.coefficient)))
    })
}
fn state_value<T: Real>(sum: &Lc<State, T>) -> Vector {
    sum.terms().fold([c(0.0, 0.0); 2], |v, t| {
        let x = vector(t.primitive);
        std::array::from_fn(|i| v[i] + complex64(t.coefficient) * x[i])
    })
}
fn assert_complex<T: Real>(actual: Complex<f64>, expected: Complex<f64>) {
    for (a, b) in [(actual.re, expected.re), (actual.im, expected.im)] {
        assert!(
            (a - b).abs() <= T::ATOL + T::RTOL * b.abs(),
            "{actual:?} != {expected:?}"
        );
    }
}
fn assert_vector<T: Real>(actual: Vector, expected: Vector) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert_complex::<T>(a, b);
    }
}
fn assert_matrix<T: Real>(actual: Matrix, expected: Matrix) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert_vector::<T>(a, b);
    }
}

// These lifting helpers are test-only. Products remain ordered and adjoints
// conjugate coefficients; no expression evaluator is added to the crate.
fn product<T: Real>(a: &Lc<Op, T>, b: &Lc<Op, T>) -> Lc<Op, T> {
    Lc::from_terms(a.terms().flat_map(|x| {
        b.terms().flat_map(move |y| {
            <Op as Product<T>>::product(&x.primitive, &y.primitive)
                .scaled(x.coefficient * y.coefficient)
        })
    }))
}
fn action<T: Real>(op: &Lc<Op, T>, state: &Lc<State, T>) -> Lc<State, T> {
    Lc::from_terms(op.terms().flat_map(|x| {
        state.terms().flat_map(move |y| {
            <Op as ActOn<State, T>>::act_on(&x.primitive, &y.primitive)
                .scaled(x.coefficient * y.coefficient)
        })
    }))
}
fn adjoint<T: Real>(op: &Lc<Op, T>) -> Lc<Op, T> {
    Lc::from_terms(
        op.terms()
            .flat_map(|t| <Op as Adjoint<T>>::adjoint(&t.primitive).scaled(t.coefficient.conj())),
    )
}
fn trace<T: Real>(op: &Lc<Op, T>) -> Complex<T> {
    op.terms().fold(c::<T>(0.0, 0.0), |v, t| {
        v + t.coefficient * <Op as Trace<T>>::trace(&t.primitive)
    })
}
fn inner<T: Real>(x: &Lc<State, T>, y: &Lc<State, T>) -> Complex<T> {
    let mut result = c::<T>(0.0, 0.0);
    for x in x {
        for y in y {
            result = result
                + x.coefficient.conj()
                    * y.coefficient
                    * <State as InnerProduct<T>>::inner_product(&x.primitive, &y.primitive);
        }
    }
    result
}
fn expand<T: Real, B: Basis<Op, T>>(basis: &B, sum: &Lc<Op, T>) -> Lc<Op, T> {
    Lc::from_terms(
        sum.terms()
            .flat_map(|t| basis.expand(&t.primitive).scaled(t.coefficient)),
    )
}
fn coordinates<T: Real, B: Basis<Op, T>>(basis: &B, sum: &Lc<Op, T>) -> [Complex<T>; 4] {
    let mut result = [c(0.0, 0.0); 4];
    for t in expand(basis, sum) {
        let i = basis.index(&t.primitive).unwrap();
        result[i] = result[i] + t.coefficient;
    }
    result
}

fn check_all_local_rules<T: Real>() {
    // Expected primitive choices/order supplement the independent matrix oracle.
    // Once these independent output primitives are fixed, matrix equality also
    // determines their coefficients uniquely.
    let product_primitives: [[&[Op]; 6]; 6] = [
        [&[I], &[Sx], &[Sy], &[Sz], &[SPlus], &[SMinus]],
        [&[Sx], &[I], &[Sz], &[Sy], &[I, Sz], &[I, Sz]],
        [&[Sy], &[Sz], &[I], &[Sx], &[I, Sz], &[I, Sz]],
        [&[Sz], &[Sy], &[Sx], &[I], &[SPlus], &[SMinus]],
        [&[SPlus], &[I, Sz], &[I, Sz], &[SPlus], &[], &[I, Sz]],
        [&[SMinus], &[I, Sz], &[I, Sz], &[SMinus], &[I, Sz], &[]],
    ];
    let action_primitives: [[&[State]; 2]; 6] = [
        [&[Up], &[Down]],
        [&[Down], &[Up]],
        [&[Down], &[Up]],
        [&[Up], &[Down]],
        [&[], &[Up]],
        [&[Down], &[]],
    ];
    let adjoints = [I, Sx, Sy, Sz, SMinus, SPlus];
    let identity = <SpinHalf as LocalStructure<T>>::identity();
    assert_eq!(identity, Lc::primitive(I));
    for (i, a) in OPS.into_iter().enumerate() {
        for (j, b) in OPS.into_iter().enumerate() {
            let actual = <Op as Product<T>>::product(&a, &b);
            assert_eq!(
                actual.terms().map(|t| t.primitive).collect::<Vec<_>>(),
                product_primitives[i][j]
            );
            assert_matrix::<T>(operator_value(&actual), multiply(matrix(a), matrix(b)));
        }
        for (j, state) in STATES.into_iter().enumerate() {
            let actual = <Op as ActOn<State, T>>::act_on(&a, &state);
            assert_eq!(
                actual.terms().map(|t| t.primitive).collect::<Vec<_>>(),
                action_primitives[i][j]
            );
            assert_vector::<T>(state_value(&actual), apply(matrix(a), vector(state)));
            assert_vector::<T>(
                state_value(&action(&identity, &Lc::primitive(state))),
                vector(state),
            );
        }
        let actual_adjoint = <Op as Adjoint<T>>::adjoint(&a);
        assert_eq!(actual_adjoint, Lc::primitive(adjoints[i]));
        assert_matrix::<T>(operator_value(&actual_adjoint), dagger(matrix(a)));
        let m = matrix(a);
        assert_complex::<T>(complex64(<Op as Trace<T>>::trace(&a)), m[0][0] + m[1][1]);
    }
    for x in STATES {
        for y in STATES {
            assert_complex::<T>(
                complex64(<State as InnerProduct<T>>::inner_product(&x, &y)),
                dot(vector(x), vector(y)),
            );
        }
    }
}

fn check_algebra_laws<T: Real>() {
    let identity = <SpinHalf as LocalStructure<T>>::identity();
    for a in OPS {
        let pa = Lc::<_, T>::primitive(a);
        assert_matrix::<T>(operator_value(&product(&identity, &pa)), matrix(a));
        assert_matrix::<T>(operator_value(&product(&pa, &identity)), matrix(a));
        assert_matrix::<T>(operator_value(&adjoint(&adjoint(&pa))), matrix(a));
        for b in OPS {
            let pb = Lc::primitive(b);
            let ab = product(&pa, &pb);
            assert_matrix::<T>(
                operator_value(&adjoint(&ab)),
                operator_value(&product(&adjoint(&pb), &adjoint(&pa))),
            );
            assert_complex::<T>(complex64(trace(&ab)), complex64(trace(&product(&pb, &pa))));
            for d in OPS {
                let pd = Lc::primitive(d);
                assert_matrix::<T>(
                    operator_value(&product(&ab, &pd)),
                    operator_value(&product(&pa, &product(&pb, &pd))),
                );
            }
            for state in STATES {
                let state = Lc::primitive(state);
                assert_vector::<T>(
                    state_value(&action(&ab, &state)),
                    state_value(&action(&pa, &action(&pb, &state))),
                );
            }
        }
        for x in STATES {
            for y in STATES {
                let x = Lc::primitive(x);
                let y = Lc::primitive(y);
                assert_complex::<T>(
                    complex64(inner(&x, &action(&pa, &y))),
                    complex64(inner(&action(&adjoint(&pa), &x), &y)),
                );
            }
        }
    }
    // Non-real coefficients exercise conjugate-linearity as well as the rules
    // on the two basis states, whose own inner products are entirely real.
    let x = Lc::primitive(Up).scaled(c::<T>(0.3, -0.7)) + Lc::primitive(Down).scaled(c(1.2, 0.4));
    let y = Lc::primitive(Up).scaled(c::<T>(-0.2, 0.6)) + Lc::primitive(Down).scaled(c(0.1, -0.4));
    let op = Lc::primitive(SPlus).scaled(c::<T>(0.4, -0.8)) + Lc::primitive(Sy).scaled(c(0.3, 0.2));
    assert_complex::<T>(
        complex64(inner(&x, &y)),
        dot(state_value(&x), state_value(&y)),
    );
    assert_complex::<T>(
        complex64(inner(&x, &action(&op, &y))),
        complex64(inner(&action(&adjoint(&op), &x), &y)),
    );
    assert_matrix::<T>(operator_value(&adjoint(&op)), dagger(operator_value(&op)));
    assert_vector::<T>(
        state_value(&action(&op, &x)),
        apply(operator_value(&op), state_value(&x)),
    );
    let m = operator_value(&op);
    assert_complex::<T>(complex64(trace(&op)), m[0][0] + m[1][1]);
}

fn check_spin_relations<T: Real>() {
    for (a, b, expected, factor) in [
        (Sx, Sy, Sz, c(0.0, 1.0)),
        (Sy, Sz, Sx, c(0.0, 1.0)),
        (Sz, Sx, Sy, c(0.0, 1.0)),
        (SPlus, SMinus, Sz, c(2.0, 0.0)),
        (Sz, SPlus, SPlus, c(1.0, 0.0)),
        (Sz, SMinus, SMinus, c(-1.0, 0.0)),
    ] {
        let ab = <Op as Product<T>>::product(&a, &b);
        let ba = <Op as Product<T>>::product(&b, &a).scaled(c(-1.0, 0.0));
        assert_matrix::<T>(operator_value(&(ab + ba)), scale(matrix(expected), factor));
    }
    for op in [Sx, Sy, Sz] {
        assert_matrix::<T>(
            operator_value(&<Op as Product<T>>::product(&op, &op)),
            scale(matrix(I), c(0.25, 0.0)),
        );
    }
    for op in [SPlus, SMinus] {
        assert!(<Op as Product<T>>::product(&op, &op).is_empty());
    }
    for (primitive, sign) in [(SPlus, 1.0), (SMinus, -1.0)] {
        let direct = Lc::<_, T>::primitive(primitive);
        let decomposition = Lc::primitive(Sx) + Lc::primitive(Sy).scaled(c(0.0, sign));
        assert_matrix::<T>(operator_value(&decomposition), matrix(primitive));
        for op in OPS {
            let op = Lc::primitive(op);
            assert_matrix::<T>(
                operator_value(&product(&decomposition, &op)),
                operator_value(&product(&direct, &op)),
            );
            assert_matrix::<T>(
                operator_value(&product(&op, &decomposition)),
                operator_value(&product(&op, &direct)),
            );
        }
        assert_matrix::<T>(
            operator_value(&adjoint(&decomposition)),
            operator_value(&adjoint(&direct)),
        );
        assert_complex::<T>(complex64(trace(&decomposition)), complex64(trace(&direct)));
        for state in STATES {
            let state = Lc::primitive(state);
            assert_vector::<T>(
                state_value(&action(&decomposition, &state)),
                state_value(&action(&direct, &state)),
            );
        }
    }
}

fn check_operator_basis<T: Real, B: OperatorBasis<SpinHalf, T>>(basis: &B, expected: [Op; 4]) {
    let elements: Vec<_> = basis.elements().collect();
    assert_eq!(elements, expected);
    assert_eq!(elements.len(), operator_basis_len::<SpinHalf, T>().unwrap());
    for op in OPS {
        assert_eq!(basis.index(&op), expected.iter().position(|&p| p == op));
        let expanded = basis.expand(&op);
        assert_eq!(
            expanded.len(),
            if basis.index(&op).is_some() { 1 } else { 2 }
        );
        assert_matrix::<T>(operator_value(&expanded), matrix(op));
        let primitives: Vec<_> = expanded.terms().map(|t| t.primitive).collect();
        if basis.index(&op).is_some() {
            assert_eq!(expanded, Lc::primitive(op));
        } else {
            assert_eq!(primitives, expected[1..3]);
        }
        for term in &expanded {
            assert!(basis.index(&term.primitive).is_some());
        }
    }
}

fn check_bases<T: Real>() {
    fn state_basis<T: Real, B: StateBasis<SpinHalf, T>>(basis: &B) {
        let elements: Vec<_> = basis.elements().collect();
        assert_eq!(elements, STATES);
        assert_eq!(elements.len(), state_basis_len::<SpinHalf, T>());
        for (i, state) in STATES.into_iter().enumerate() {
            assert_eq!(basis.index(&state), Some(i));
            assert_eq!(basis.expand(&state), Lc::primitive(state));
        }
    }
    state_basis::<T, _>(&SzBasis);
    check_operator_basis::<T, _>(&CartesianBasis, [I, Sx, Sy, Sz]);
    check_operator_basis::<T, _>(&LadderBasis, [I, SPlus, SMinus, Sz]);
    for op in OPS {
        let original = Lc::<_, T>::primitive(op);
        let cartesian = expand(&CartesianBasis, &original);
        let ladder = expand(&LadderBasis, &original);
        let from_ladder = expand(&CartesianBasis, &ladder);
        let from_cartesian = expand(&LadderBasis, &cartesian);
        for (a, b) in coordinates(&CartesianBasis, &from_ladder)
            .into_iter()
            .zip(coordinates(&CartesianBasis, &original))
        {
            assert_complex::<T>(complex64(a), complex64(b));
        }
        for (a, b) in coordinates(&LadderBasis, &from_cartesian)
            .into_iter()
            .zip(coordinates(&LadderBasis, &original))
        {
            assert_complex::<T>(complex64(a), complex64(b));
        }
    }
    let sum = Lc::<_, T>::primitive(SPlus) + Lc::primitive(SMinus);
    for (a, b) in coordinates(&CartesianBasis, &sum).into_iter().zip([
        c(0.0, 0.0),
        c(2.0, 0.0),
        c(0.0, 0.0),
        c(0.0, 0.0),
    ]) {
        assert_complex::<T>(complex64(a), b);
    }
}

#[test]
fn all_local_rules_at_both_precisions() {
    check_all_local_rules::<f32>();
    check_all_local_rules::<f64>();
}
#[test]
fn algebra_laws_at_both_precisions() {
    check_algebra_laws::<f32>();
    check_algebra_laws::<f64>();
}
#[test]
fn spin_relations_at_both_precisions() {
    check_spin_relations::<f32>();
    check_spin_relations::<f64>();
}
#[test]
fn basis_contracts_and_round_trips_at_both_precisions() {
    check_bases::<f32>();
    check_bases::<f64>();
}
#[test]
fn precision_independent_primitives_and_public_type_contracts() {
    fn same_primitives<
        L: LocalStructure<f32>
            + LocalStructure<
                f64,
                State = <L as LocalStructure<f32>>::State,
                Operator = <L as LocalStructure<f32>>::Operator,
            >,
    >() {
    }
    same_primitives::<SpinHalf>();
    assert_eq!(
        <SpinHalf as LocalStructure<f32>>::DIMENSION,
        <SpinHalf as LocalStructure<f64>>::DIMENSION
    );
    assert_eq!(<SpinHalf as LocalStructure<f64>>::DIMENSION.get(), 2);
    fn discrete<T: Clone + Copy + Debug + Eq + Hash>() {}
    fn parameter_free<T: Default + Clone + Copy + Debug + Eq + Hash>() {}
    discrete::<State>();
    discrete::<Op>();
    parameter_free::<SpinHalf>();
    parameter_free::<SzBasis>();
    parameter_free::<CartesianBasis>();
    parameter_free::<LadderBasis>();
    let _: qmb_spin_half::local_structure::SpinHalf = SpinHalf;
    let _: qmb_spin_half::state::SpinHalfState = Up;
    let _: qmb_spin_half::operator::SpinHalfOperator = I;
    let _: qmb_spin_half::basis::CartesianBasis = CartesianBasis;
}
