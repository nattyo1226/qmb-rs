//! Public contracts checked against independent two-dimensional linear algebra.

use std::{fmt::Debug, hash::Hash, ops::Neg};

use num_complex::Complex;
use num_traits::{FromPrimitive, Num, ToPrimitive};
use qmb_core::{
    Act, Adjoint, Basis, Identity, InnerProduct, LinearCombination as Lc, Product, SpanDimension,
    Trace, operator_basis_len, state_basis_len,
};
use qmb_spin_half::{
    CartesianBasis, LadderBasis, SpinHalfOperator as Op, SpinHalfState as State, SzBasis,
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
fn operator_value<T: Real>(sum: &Lc<Op, Complex<T>>) -> Matrix {
    sum.terms().fold([[c(0.0, 0.0); 2]; 2], |v, t| {
        add(v, scale(matrix(t.operand), complex64(t.coefficient)))
    })
}
fn state_value<T: Real>(sum: &Lc<State, Complex<T>>) -> Vector {
    sum.terms().fold([c(0.0, 0.0); 2], |v, t| {
        let x = vector(t.operand);
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
// conjugate coefficients independently of the production expression evaluator.
fn product<T: Real>(a: &Lc<Op, Complex<T>>, b: &Lc<Op, Complex<T>>) -> Lc<Op, Complex<T>> {
    Lc::from_terms(a.terms().flat_map(|x| {
        b.terms().flat_map(move |y| {
            <Op as Product<Complex<T>>>::product(&x.operand, &y.operand)
                .scaled(x.coefficient * y.coefficient)
        })
    }))
}
fn action<T: Real>(
    op: &Lc<Op, Complex<T>>,
    state: &Lc<State, Complex<T>>,
) -> Lc<State, Complex<T>> {
    Lc::from_terms(op.terms().flat_map(|x| {
        state.terms().flat_map(move |y| {
            <State as Act<Complex<T>>>::act(&y.operand, &x.operand)
                .scaled(x.coefficient * y.coefficient)
        })
    }))
}
fn adjoint<T: Real>(op: &Lc<Op, Complex<T>>) -> Lc<Op, Complex<T>> {
    Lc::from_terms(op.terms().flat_map(|t| {
        <Op as Adjoint<Complex<T>>>::adjoint(&t.operand).scaled(t.coefficient.conj())
    }))
}
fn trace<T: Real>(op: &Lc<Op, Complex<T>>) -> Complex<T> {
    op.terms().fold(c::<T>(0.0, 0.0), |v, t| {
        v + t.coefficient * <Op as Trace<Complex<T>>>::trace(&t.operand)
    })
}
fn inner<T: Real>(x: &Lc<State, Complex<T>>, y: &Lc<State, Complex<T>>) -> Complex<T> {
    let mut result = c::<T>(0.0, 0.0);
    for x in x {
        for y in y {
            result = result
                + x.coefficient.conj()
                    * y.coefficient
                    * <State as InnerProduct<Complex<T>>>::inner_product(&x.operand, &y.operand);
        }
    }
    result
}
fn expand<T: Real, B: Basis<Op, Complex<T>>>(
    basis: &B,
    sum: &Lc<Op, Complex<T>>,
) -> Lc<Op, Complex<T>> {
    Lc::from_terms(
        sum.terms()
            .flat_map(|t| basis.expand(&t.operand).scaled(t.coefficient)),
    )
}
fn coordinates<T: Real, B: Basis<Op, Complex<T>>>(
    basis: &B,
    sum: &Lc<Op, Complex<T>>,
) -> [Complex<T>; 4] {
    let mut result = [c(0.0, 0.0); 4];
    for t in expand(basis, sum) {
        let i = basis.index(&t.operand).unwrap();
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
    let identity = Lc::<_, Complex<T>>::primitive(Op::identity());
    assert_eq!(identity, Lc::primitive(I));
    for (i, a) in OPS.into_iter().enumerate() {
        for (j, b) in OPS.into_iter().enumerate() {
            let actual = <Op as Product<Complex<T>>>::product(&a, &b);
            assert_eq!(
                actual.terms().map(|t| t.operand).collect::<Vec<_>>(),
                product_primitives[i][j]
            );
            assert_matrix::<T>(operator_value(&actual), multiply(matrix(a), matrix(b)));
        }
        for (j, state) in STATES.into_iter().enumerate() {
            let actual = <State as Act<Complex<T>>>::act(&state, &a);
            assert_eq!(
                actual.terms().map(|t| t.operand).collect::<Vec<_>>(),
                action_primitives[i][j]
            );
            assert_vector::<T>(state_value(&actual), apply(matrix(a), vector(state)));
            assert_vector::<T>(
                state_value(&action(&identity, &Lc::primitive(state))),
                vector(state),
            );
        }
        let actual_adjoint = <Op as Adjoint<Complex<T>>>::adjoint(&a);
        assert_eq!(actual_adjoint, Lc::primitive(adjoints[i]));
        assert_matrix::<T>(operator_value(&actual_adjoint), dagger(matrix(a)));
        let m = matrix(a);
        assert_complex::<T>(
            complex64(<Op as Trace<Complex<T>>>::trace(&a)),
            m[0][0] + m[1][1],
        );
    }
    for x in STATES {
        for y in STATES {
            assert_complex::<T>(
                complex64(<State as InnerProduct<Complex<T>>>::inner_product(&x, &y)),
                dot(vector(x), vector(y)),
            );
        }
    }
}

fn check_algebra_laws<T: Real>() {
    let identity = Lc::<_, Complex<T>>::primitive(Op::identity());
    for a in OPS {
        let pa = Lc::<_, Complex<T>>::primitive(a);
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
    let x = Lc::<_, Complex<T>>::primitive(Up).scaled(c::<T>(0.3, -0.7))
        + Lc::<_, Complex<T>>::primitive(Down).scaled(c::<T>(1.2, 0.4));
    let y = Lc::<_, Complex<T>>::primitive(Up).scaled(c::<T>(-0.2, 0.6))
        + Lc::<_, Complex<T>>::primitive(Down).scaled(c::<T>(0.1, -0.4));
    let op = Lc::<_, Complex<T>>::primitive(SPlus).scaled(c::<T>(0.4, -0.8))
        + Lc::<_, Complex<T>>::primitive(Sy).scaled(c::<T>(0.3, 0.2));
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
        let ab = <Op as Product<Complex<T>>>::product(&a, &b);
        let ba = <Op as Product<Complex<T>>>::product(&b, &a).scaled(c(-1.0, 0.0));
        assert_matrix::<T>(operator_value(&(ab + ba)), scale(matrix(expected), factor));
    }
    for op in [Sx, Sy, Sz] {
        assert_matrix::<T>(
            operator_value(&<Op as Product<Complex<T>>>::product(&op, &op)),
            scale(matrix(I), c(0.25, 0.0)),
        );
    }
    for op in [SPlus, SMinus] {
        assert!(<Op as Product<Complex<T>>>::product(&op, &op).is_empty());
    }
    for (primitive, sign) in [(SPlus, 1.0), (SMinus, -1.0)] {
        let direct = Lc::<_, Complex<T>>::primitive(primitive);
        let decomposition =
            Lc::primitive(Sx) + Lc::<_, Complex<T>>::primitive(Sy).scaled(c(0.0, sign));
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

fn check_operator_basis<T: Real, B: Basis<Op, Complex<T>>>(basis: &B, expected: [Op; 4]) {
    let elements: Vec<_> = basis.elements().collect();
    assert_eq!(elements, expected);
    assert_eq!(elements.len(), operator_basis_len::<Op>().unwrap());
    for op in OPS {
        assert_eq!(basis.index(&op), expected.iter().position(|&p| p == op));
        let expanded = basis.expand(&op);
        assert_eq!(
            expanded.len(),
            if basis.index(&op).is_some() { 1 } else { 2 }
        );
        assert_matrix::<T>(operator_value(&expanded), matrix(op));
        let primitives: Vec<_> = expanded.terms().map(|t| t.operand).collect();
        if basis.index(&op).is_some() {
            assert_eq!(expanded, Lc::primitive(op));
        } else {
            assert_eq!(primitives, expected[1..3]);
        }
        for term in &expanded {
            assert!(basis.index(&term.operand).is_some());
        }
    }
}

fn check_bases<T: Real>() {
    fn state_basis<T: Real, B: Basis<State, Complex<T>>>(basis: &B) {
        let elements: Vec<_> = basis.elements().collect();
        assert_eq!(elements, STATES);
        assert_eq!(elements.len(), state_basis_len::<State>().unwrap());
        for (i, state) in STATES.into_iter().enumerate() {
            assert_eq!(basis.index(&state), Some(i));
            assert_eq!(basis.expand(&state), Lc::primitive(state));
        }
    }
    state_basis::<T, _>(&SzBasis);
    check_operator_basis::<T, _>(&CartesianBasis, [I, Sx, Sy, Sz]);
    check_operator_basis::<T, _>(&LadderBasis, [I, SPlus, SMinus, Sz]);
    for op in OPS {
        let original = Lc::<_, Complex<T>>::primitive(op);
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
    let sum = Lc::<_, Complex<T>>::primitive(SPlus) + Lc::primitive(SMinus);
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
    assert_eq!(State::span_dimension().unwrap().get(), 2);
    assert_eq!(Op::span_dimension().unwrap().get(), 4);
    fn discrete<T: Clone + Copy + Debug + Eq + Hash>() {}
    fn parameter_free<T: Default + Clone + Copy + Debug + Eq + Hash>() {}
    discrete::<State>();
    discrete::<Op>();
    parameter_free::<SzBasis>();
    parameter_free::<CartesianBasis>();
    parameter_free::<LadderBasis>();
    let _: qmb_spin_half::state::SpinHalfState = Up;
    let _: qmb_spin_half::operator::SpinHalfOperator = I;
    let _: qmb_spin_half::basis::CartesianBasis = CartesianBasis;
}

// Expression evaluation is checked against the literal matrices above, not the
// primitive lifting helpers. Projections also exercise inner-product evaluation.
fn expression_vector<T: Real>(
    state: &qmb_core::LocalState<
        qmb_spin_half::SpinHalfState,
        qmb_spin_half::SpinHalfOperator,
        Complex<T>,
    >,
) -> Vector {
    use qmb_core::LocalState;
    STATES.map(|basis| {
        complex64(LocalState::<qmb_spin_half::SpinHalfState, qmb_spin_half::SpinHalfOperator, Complex<T>>::primitive(basis).inner_product(state))
    })
}
fn expression_matrix<T: Real>(
    op: &qmb_core::LocalOperator<qmb_spin_half::SpinHalfOperator, Complex<T>>,
) -> Matrix {
    use qmb_core::LocalState;
    let columns = STATES.map(|basis| expression_vector(&LocalState::primitive(basis).act(op)));
    std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i]))
}
fn matrix_trace(matrix: Matrix) -> Complex<f64> {
    matrix[0][0] + matrix[1][1]
}

fn check_local_expressions<T: Real>() {
    use qmb_core::{LocalOperator, LocalState};
    type O<T> = LocalOperator<qmb_spin_half::SpinHalfOperator, Complex<T>>;
    type S<T> =
        LocalState<qmb_spin_half::SpinHalfState, qmb_spin_half::SpinHalfOperator, Complex<T>>;
    let up = S::<T>::primitive(Up);
    let down = S::<T>::primitive(Down);
    let state = up.scaled(c(0.7, 0.2)) + down.scaled(c(-0.3, 0.5));
    let reference = [c::<f64>(0.7, 0.2), c(-0.3, 0.5)];
    assert_vector::<T>(expression_vector(&state), reference);
    assert_complex::<T>(
        complex64(state.inner_product(&state)),
        dot(reference, reference),
    );
    for a in OPS {
        let a_expr = O::<T>::primitive(a).scaled(c(0.4, -0.2)) + O::identity().scaled(c(-0.3, 0.1));
        let a_matrix = add(
            scale(matrix(a), c(0.4, -0.2)),
            scale(matrix(I), c(-0.3, 0.1)),
        );
        assert_matrix::<T>(expression_matrix(&a_expr), a_matrix);
        assert_complex::<T>(complex64(a_expr.trace()), matrix_trace(a_matrix));
        assert_complex::<T>(
            complex64(state.expectation(&a_expr)),
            dot(reference, apply(a_matrix, reference)),
        );
        for b in OPS {
            let b_expr = O::<T>::primitive(b) + O::primitive(Sy).scaled(c(0.2, 0.3));
            let b_matrix = add(matrix(b), scale(matrix(Sy), c(0.2, 0.3)));
            let product = a_expr.product(&b_expr);
            let product_matrix = multiply(a_matrix, b_matrix);
            assert_matrix::<T>(expression_matrix(&product), product_matrix);
            assert_matrix::<T>(
                expression_matrix(&product.adjoint()),
                dagger(product_matrix),
            );
            assert_complex::<T>(complex64(product.trace()), matrix_trace(product_matrix));
            assert_complex::<T>(
                complex64(product.adjoint().trace()),
                matrix_trace(dagger(product_matrix)),
            );
            assert_complex::<T>(
                complex64(a_expr.hilbert_schmidt_inner_product(&b_expr)),
                matrix_trace(multiply(dagger(a_matrix), b_matrix)),
            );
            assert_matrix::<T>(
                expression_matrix(&a_expr.commutator(&b_expr)),
                add(
                    product_matrix,
                    scale(multiply(b_matrix, a_matrix), c(-1.0, 0.0)),
                ),
            );
            assert_matrix::<T>(
                expression_matrix(&a_expr.anticommutator(&b_expr)),
                add(product_matrix, multiply(b_matrix, a_matrix)),
            );
            let chain = state.act(&b_expr).act(&a_expr);
            let result = chain.evaluate();
            assert_vector::<T>(expression_vector(&chain), apply(product_matrix, reference));
            assert_vector::<T>(expression_vector(&result), apply(product_matrix, reference));
            assert_vector::<T>(
                expression_vector(&result.evaluate()),
                apply(product_matrix, reference),
            );
            assert_complex::<T>(
                complex64(chain.inner_product(&state.act(&b_expr))),
                dot(apply(product_matrix, reference), apply(b_matrix, reference)),
            );
            assert_vector::<T>(expression_vector(&state), reference);
        }
    }
    for basis in [CartesianBasis] {
        for primitive in OPS {
            let expanded = <CartesianBasis as Basis<Op, Complex<T>>>::expand(&basis, &primitive);
            let expression = O::from_linear_combination(expanded.map(O::primitive));
            assert_matrix::<T>(expression_matrix(&expression), matrix(primitive));
        }
    }
    for primitive in OPS {
        let expanded = <LadderBasis as Basis<Op, Complex<T>>>::expand(&LadderBasis, &primitive);
        assert_matrix::<T>(
            expression_matrix(&O::from_linear_combination(expanded.map(O::primitive))),
            matrix(primitive),
        );
    }
    for primitive in STATES {
        let expanded = <SzBasis as Basis<State, Complex<T>>>::expand(&SzBasis, &primitive);
        assert_vector::<T>(
            expression_vector(&S::from_linear_combination(expanded.map(S::primitive))),
            vector(primitive),
        );
    }
    for op in [
        O::<T>::zero(),
        O::sum([]),
        O::primitive(SPlus).product(&O::primitive(SPlus)),
    ] {
        assert_matrix::<T>(expression_matrix(&op), [[c(0.0, 0.0); 2]; 2]);
        assert_complex::<T>(complex64(op.trace()), c(0.0, 0.0));
    }
    assert_complex::<T>(complex64(O::<T>::product_all([]).trace()), c(2.0, 0.0));
    assert_complex::<T>(
        complex64(O::<T>::identity().hilbert_schmidt_inner_product(&O::identity())),
        c(2.0, 0.0),
    );
    assert_vector::<T>(expression_vector(&S::<T>::zero()), [c(0.0, 0.0); 2]);
}

#[test]
fn local_expressions_match_independent_matrices_at_both_precisions() {
    check_local_expressions::<f32>();
    check_local_expressions::<f64>();
}

#[test]
#[allow(
    clippy::op_ref,
    reason = "Exercise the public borrowed scalar overloads"
)]
fn all_expression_operator_ownership_forms_agree() {
    use qmb_core::{LocalOperator, LocalState};
    type O = LocalOperator<qmb_spin_half::SpinHalfOperator, Complex<f64>>;
    type S =
        LocalState<qmb_spin_half::SpinHalfState, qmb_spin_half::SpinHalfOperator, Complex<f64>>;
    let a = O::primitive(SPlus);
    let b = O::primitive(SMinus);
    let x = S::primitive(Up);
    let y = S::primitive(Down);
    for expression in [
        a.clone() + b.clone(),
        a.clone() + &b,
        &a + b.clone(),
        &a + &b,
    ] {
        assert_matrix::<f64>(
            expression_matrix(&expression),
            add(matrix(SPlus), matrix(SMinus)),
        );
    }
    for expression in [
        a.clone() - b.clone(),
        a.clone() - &b,
        &a - b.clone(),
        &a - &b,
    ] {
        assert_matrix::<f64>(
            expression_matrix(&expression),
            add(matrix(SPlus), scale(matrix(SMinus), c(-1.0, 0.0))),
        );
    }
    for expression in [
        a.clone() * b.clone(),
        a.clone() * &b,
        &a * b.clone(),
        &a * &b,
    ] {
        assert_matrix::<f64>(
            expression_matrix(&expression),
            multiply(matrix(SPlus), matrix(SMinus)),
        );
    }
    for expression in [
        a.clone() * y.clone(),
        a.clone() * &y,
        &a * y.clone(),
        &a * &y,
    ] {
        assert_vector::<f64>(expression_vector(&expression), vector(Up));
    }
    for expression in [
        x.clone() + y.clone(),
        x.clone() + &y,
        &x + y.clone(),
        &x + &y,
    ] {
        assert_vector::<f64>(expression_vector(&expression), [c(1.0, 0.0); 2]);
    }
    for expression in [
        x.clone() - y.clone(),
        x.clone() - &y,
        &x - y.clone(),
        &x - &y,
    ] {
        assert_vector::<f64>(expression_vector(&expression), [c(1.0, 0.0), c(-1.0, 0.0)]);
    }
    let coefficient = c(0.3, 0.7);
    for expression in [
        a.clone() * coefficient,
        &a * coefficient,
        coefficient * a.clone(),
        coefficient * &a,
        a.clone() * &coefficient,
        &a * &coefficient,
        &coefficient * a.clone(),
        &coefficient * &a,
    ] {
        assert_matrix::<f64>(
            expression_matrix(&expression),
            scale(matrix(SPlus), coefficient),
        );
    }
    for expression in [
        x.clone() * coefficient,
        &x * coefficient,
        coefficient * x.clone(),
        coefficient * &x,
        x.clone() * &coefficient,
        &x * &coefficient,
        &coefficient * x.clone(),
        &coefficient * &x,
    ] {
        assert_vector::<f64>(expression_vector(&expression), [coefficient, c(0.0, 0.0)]);
    }
    for expression in [-a.clone(), -&a] {
        assert_matrix::<f64>(
            expression_matrix(&expression),
            scale(matrix(SPlus), c(-1.0, 0.0)),
        );
    }
    for expression in [-x.clone(), -&x] {
        assert_vector::<f64>(expression_vector(&expression), [c(-1.0, 0.0), c(0.0, 0.0)]);
    }
}

#[test]
fn long_operator_trace_and_destruction_fit_a_small_stack() {
    std::thread::Builder::new()
        .stack_size(128 * 1024)
        .spawn(|| {
            use qmb_core::LocalOperator;
            type O = LocalOperator<qmb_spin_half::SpinHalfOperator, Complex<f64>>;
            let mut expression = O::primitive(Sz);
            for _ in 0..20_000 {
                expression = expression
                    .product(&O::identity())
                    .adjoint()
                    .scaled(c(1.0, 0.0));
            }
            assert_complex::<f64>(expression.trace(), c(0.0, 0.0));
            drop(expression);

            let factor = GlobalSpace::factor(0);
            let mut placed = O::primitive(Sz).at(&factor).unwrap();
            for _ in 0..20_000 {
                placed = placed
                    .product(&GlobalOperator::identity())
                    .unwrap()
                    .adjoint()
                    .scaled(c(1.0, 0.0));
            }
            assert_complex::<f64>(placed.trace_in(&factor).unwrap(), c(0.0, 0.0));
            drop(placed);
        })
        .unwrap()
        .join()
        .unwrap();
}

type GlobalState<T> =
    qmb_core::State<u8, qmb_spin_half::SpinHalfState, qmb_spin_half::SpinHalfOperator, Complex<T>>;
type GlobalOperator<T> = qmb_core::Operator<u8, qmb_spin_half::SpinHalfOperator, Complex<T>>;
type GlobalSpace = qmb_core::HilbertSpace<u8>;
type Dense = Vec<Vec<Complex<f64>>>;

fn placed_matrix(operator: Op, site: usize) -> Dense {
    // Literal Kronecker reference, with the rightmost factor varying fastest.
    (0..8)
        .map(|i| {
            (0..8)
                .map(|j| {
                    (0..3).fold(c(1.0, 0.0), |value, k| {
                        value
                            * matrix(if site == k { operator } else { I })[(i >> (2 - k)) & 1]
                                [(j >> (2 - k)) & 1]
                    })
                })
                .collect()
        })
        .collect()
}
fn dense_add(a: &Dense, b: &Dense) -> Dense {
    a.iter()
        .zip(b)
        .map(|(a, b)| a.iter().zip(b).map(|(a, b)| a + b).collect())
        .collect()
}
fn dense_scale(a: &Dense, c: Complex<f64>) -> Dense {
    a.iter()
        .map(|row| row.iter().map(|a| c * a).collect())
        .collect()
}
fn dense_mul(a: &Dense, b: &Dense) -> Dense {
    (0..a.len())
        .map(|i| {
            (0..a.len())
                .map(|j| (0..a.len()).fold(c(0.0, 0.0), |s, k| s + a[i][k] * b[k][j]))
                .collect()
        })
        .collect()
}
fn dense_dagger(a: &Dense) -> Dense {
    (0..a.len())
        .map(|i| (0..a.len()).map(|j| a[j][i].conj()).collect())
        .collect()
}
fn dense_apply(a: &Dense, x: &[Complex<f64>]) -> Vec<Complex<f64>> {
    a.iter()
        .map(|row| row.iter().zip(x).fold(c(0.0, 0.0), |s, (a, x)| s + a * x))
        .collect()
}
fn dense_trace(a: &Dense) -> Complex<f64> {
    (0..a.len()).fold(c(0.0, 0.0), |s, i| s + a[i][i])
}
fn dense_dot(a: &[Complex<f64>], b: &[Complex<f64>]) -> Complex<f64> {
    a.iter()
        .zip(b)
        .fold(c(0.0, 0.0), |s, (a, b)| s + a.conj() * b)
}

fn global_basis<T: Real>(factors: &[GlobalSpace; 3], n: usize) -> GlobalState<T> {
    qmb_core::State::tensor(factors.iter().enumerate().map(|(k, factor)| {
        qmb_core::LocalState::primitive(if (n >> (2 - k)) & 1 == 0 { Up } else { Down })
            .at(factor)
            .unwrap()
    }))
    .unwrap()
}
fn global_vector<T: Real>(state: &GlobalState<T>, factors: &[GlobalSpace; 3]) -> Vec<Complex<f64>> {
    let state = state.evaluate().unwrap();
    (0..8)
        .map(|n| complex64(global_basis(factors, n).inner_product(&state).unwrap()))
        .collect()
}
fn assert_dense_vector<T: Real>(a: &[Complex<f64>], b: &[Complex<f64>]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_complex::<T>(*a, *b);
    }
}

fn check_composite_expressions<T: Real>() {
    use qmb_core::{LocalOperator as LO, LocalState as LS};
    let factors = [
        GlobalSpace::factor(9),
        GlobalSpace::factor(2),
        GlobalSpace::factor(7),
    ];
    let h = factors[0]
        .tensor_product(&factors[1])
        .unwrap()
        .tensor_product(&factors[2])
        .unwrap();
    let [a, b, d] = &factors;
    let local =
        LS::<qmb_spin_half::SpinHalfState, qmb_spin_half::SpinHalfOperator, Complex<T>>::primitive(
            Up,
        )
        .scaled(c(0.7, 0.2))
            + LS::primitive(Down).scaled(c(-0.3, 0.5));
    let s0 = local.at(a).unwrap();
    let s1 = LS::primitive(Down).at(b).unwrap();
    let s2 = LS::primitive(Up).at(d).unwrap();
    let left = s0.tensor_product(&s1).unwrap().tensor_product(&s2).unwrap();
    let right = s0.tensor_product(&s1.tensor_product(&s2).unwrap()).unwrap();
    assert!(left.space().is_compatible_with(right.space()));
    assert_dense_vector::<T>(
        &global_vector(&left, &factors),
        &global_vector(&right, &factors),
    );
    let state = left
        .add(&global_basis(&factors, 5).scaled(c(-0.4, 0.6)))
        .unwrap();
    let mut reference = vec![c(0.0, 0.0); 8];
    reference[2] = c(0.7, 0.2);
    reference[6] = c(-0.3, 0.5);
    reference[5] = c(-0.4, 0.6);
    assert_dense_vector::<T>(&global_vector(&state, &factors), &reference);
    assert_complex::<T>(
        complex64(state.inner_product(&state).unwrap()),
        dense_dot(&reference, &reference),
    );
    for oa in OPS {
        for ob in OPS {
            let a = LO::<qmb_spin_half::SpinHalfOperator, Complex<T>>::primitive(oa)
                .at(d)
                .unwrap()
                .scaled(c(0.4, -0.2))
                .add(&GlobalOperator::scalar(c(-0.3, 0.1)))
                .unwrap();
            let b = LO::primitive(ob)
                .at(&factors[0])
                .unwrap()
                .add(
                    &LO::primitive(Sy)
                        .at(&factors[1])
                        .unwrap()
                        .scaled(c(0.2, 0.3)),
                )
                .unwrap();
            let am = dense_add(
                &dense_scale(&placed_matrix(oa, 2), c(0.4, -0.2)),
                &dense_scale(&placed_matrix(I, 0), c(-0.3, 0.1)),
            );
            let bm = dense_add(
                &placed_matrix(ob, 0),
                &dense_scale(&placed_matrix(Sy, 1), c(0.2, 0.3)),
            );
            // Include repeated action on the same factor as well as separated ones.
            let same = LO::primitive(SPlus).at(d).unwrap();
            let expression = a.product(&b).unwrap().product(&same).unwrap();
            let expected = dense_mul(&dense_mul(&am, &bm), &placed_matrix(SPlus, 2));
            let result = state.act(&expression).unwrap().evaluate().unwrap();
            assert_dense_vector::<T>(
                &global_vector(&result, &factors),
                &dense_apply(&expected, &reference),
            );
            let sequential = state.act(&same).unwrap().act(&b).unwrap().act(&a).unwrap();
            assert_dense_vector::<T>(
                &global_vector(&sequential, &factors),
                &dense_apply(&expected, &reference),
            );
            assert_dense_vector::<T>(
                &global_vector(&state.act(&expression.adjoint()).unwrap(), &factors),
                &dense_apply(&dense_dagger(&expected), &reference),
            );
            assert_complex::<T>(
                complex64(expression.trace_in(&h).unwrap()),
                dense_trace(&expected),
            );
            assert_complex::<T>(
                complex64(expression.adjoint().trace_in(&h).unwrap()),
                dense_trace(&dense_dagger(&expected)),
            );
            assert_complex::<T>(
                complex64(state.expectation(&expression).unwrap()),
                dense_dot(&reference, &dense_apply(&expected, &reference)),
            );
            assert_complex::<T>(
                complex64(a.hilbert_schmidt_inner_product_in(&b, &h).unwrap()),
                dense_trace(&dense_mul(&dense_dagger(&am), &bm)),
            );
            let ab = dense_mul(&am, &bm);
            let ba = dense_mul(&bm, &am);
            for (op, matrix) in [
                (
                    a.commutator(&b).unwrap(),
                    dense_add(&ab, &dense_scale(&ba, c(-1.0, 0.0))),
                ),
                (a.anticommutator(&b).unwrap(), dense_add(&ab, &ba)),
            ] {
                assert_dense_vector::<T>(
                    &global_vector(&state.act(&op).unwrap(), &factors),
                    &dense_apply(&matrix, &reference),
                );
            }
        }
    }
    assert_dense_vector::<T>(&global_vector(&state, &factors), &reference);
}

#[test]
fn composite_expressions_match_independent_three_factor_matrices() {
    check_composite_expressions::<f32>();
    check_composite_expressions::<f64>();
}

fn check_subsystem_contractions<T: Real>() {
    use qmb_core::{LocalOperator, LocalState};
    let factors = [9, 2, 7].map(GlobalSpace::factor);
    let ket = |site: usize, state| {
        LocalState::<qmb_spin_half::SpinHalfState, qmb_spin_half::SpinHalfOperator, Complex<T>>::primitive(state)
            .at(&factors[site])
            .unwrap()
    };
    let ab = ket(0, Up)
        .tensor_product(&ket(1, Up))
        .unwrap()
        .add(
            &ket(0, Down)
                .tensor_product(&ket(1, Down))
                .unwrap()
                .scaled(c(0.7, 0.2)),
        )
        .unwrap();
    let last = ket(2, Up).add(&ket(2, Down).scaled(c(-0.3, 0.5))).unwrap();
    let lhs = ab.tensor_product(&last).unwrap();
    let first = ket(0, Up).add(&ket(0, Down).scaled(c(0.4, -0.6))).unwrap();
    let bc = ket(1, Up)
        .tensor_product(&ket(2, Up))
        .unwrap()
        .add(
            &ket(1, Down)
                .tensor_product(&ket(2, Down))
                .unwrap()
                .scaled(c(-0.2, 0.8)),
        )
        .unwrap();
    let rhs = qmb_core::State::tensor([
        qmb_core::State::unit(),
        first,
        qmb_core::State::tensor([bc]).unwrap(),
    ])
    .unwrap();

    // Literal coordinates independently describe (AB)C and A(BC). Both AB and
    // BC are entangled sums, so their subsystem boundaries overlap.
    let mut left = vec![c::<f64>(0.0, 0.0); 8];
    left[0] = c(1.0, 0.0);
    left[1] = c(-0.3, 0.5);
    left[6] = c(0.7, 0.2);
    left[7] = left[6] * left[1];
    let mut right = vec![c::<f64>(0.0, 0.0); 8];
    right[0] = c(1.0, 0.0);
    right[3] = c(-0.2, 0.8);
    right[4] = c(0.4, -0.6);
    right[7] = right[3] * right[4];
    let op = LocalOperator::primitive(Sx).at(&factors[1]).unwrap();
    let states = [lhs, rhs];
    let vectors = [left, right];
    for (i, state) in states.iter().enumerate() {
        for (j, other) in states.iter().enumerate() {
            assert_complex::<T>(
                complex64(state.inner_product(other).unwrap()),
                dense_dot(&vectors[i], &vectors[j]),
            );
        }
        assert_complex::<T>(
            complex64(state.expectation(&op).unwrap()),
            dense_dot(
                &vectors[i],
                &dense_apply(&placed_matrix(Sx, 1), &vectors[i]),
            ),
        );
    }
}

#[test]
fn subsystem_contractions_match_independent_vectors_at_both_precisions() {
    check_subsystem_contractions::<f32>();
    check_subsystem_contractions::<f64>();
}

#[test]
fn composite_empty_zero_spectator_and_overflow_contracts() {
    use qmb_core::{
        ExpressionError as E, LocalOperator as LO, LocalState as LS, Operator as O, State as S,
    };
    type Op = O<u32, qmb_spin_half::SpinHalfOperator, Complex<f64>>;
    type Ket = S<u32, qmb_spin_half::SpinHalfState, qmb_spin_half::SpinHalfOperator, Complex<f64>>;
    let unit = Ket::unit();
    let identity = Op::identity();
    assert_eq!(identity.trace_in(unit.space()).unwrap(), c(1.0, 0.0));
    assert_eq!(unit.inner_product(&unit).unwrap(), c(1.0, 0.0));
    assert_eq!(
        unit.act(&identity).unwrap().inner_product(&unit).unwrap(),
        c(1.0, 0.0)
    );
    assert_eq!(
        unit.act(&Op::zero()).unwrap().inner_product(&unit).unwrap(),
        c(0.0, 0.0)
    );
    assert!(matches!(identity.trace(), Err(E::UnboundOperator)));
    let factors: Vec<_> = (0..usize::BITS + 1)
        .map(qmb_core::HilbertSpace::<_>::factor)
        .collect();
    let space = qmb_core::HilbertSpace::tensor(factors.iter().cloned()).unwrap();
    let state = Ket::tensor(factors.iter().map(|f| LS::primitive(Up).at(f).unwrap())).unwrap();
    assert!(space.dimension::<State>().is_err());
    assert_eq!(
        identity.trace_in(&space).unwrap(),
        c(2f64.powi(factors.len() as i32), 0.0)
    );
    assert_eq!(state.inner_product(&state).unwrap(), c(1.0, 0.0));
    let local = LO::primitive(Sz).at(factors.last().unwrap()).unwrap();
    assert_eq!(state.expectation(&local).unwrap(), c(0.5, 0.0));
    let zero = state
        .act(&LO::primitive(SPlus).at(&factors[0]).unwrap())
        .unwrap()
        .evaluate()
        .unwrap();
    assert!(zero.space().is_compatible_with(&space));
    assert_eq!(zero.inner_product(&zero).unwrap(), c(0.0, 0.0));
    let scalar = Op::scalar(c(0.3, 0.2)).bind(&space).unwrap();
    assert_eq!(state.expectation(&scalar).unwrap(), c(0.3, 0.2));
    let nested = Op::identity().bind(&factors[0]).unwrap();
    assert_eq!(
        nested.trace_in(&space).unwrap(),
        identity.trace_in(&space).unwrap()
    );
    assert!(nested.space().unwrap().is_compatible_with(&factors[0]));

    let placed_identity = LO::<qmb_spin_half::SpinHalfOperator, Complex<f64>>::primitive(I)
        .at(&factors[0])
        .unwrap()
        .bind(&factors[0])
        .unwrap();
    assert_eq!(placed_identity.trace().unwrap(), c(2.0, 0.0));
    assert_eq!(
        placed_identity.trace_in(&space).unwrap(),
        identity.trace_in(&space).unwrap()
    );
    assert_eq!(placed_identity.trace().unwrap(), c(2.0, 0.0));
    let bound_zero = Op::sum([]).unwrap().bind(&space).unwrap();
    assert_eq!(bound_zero.trace().unwrap(), c(0.0, 0.0));
    assert!(bound_zero.domain().unwrap().shares_node_with(&space));
    assert!(bound_zero.codomain().unwrap().shares_node_with(&space));
    assert_eq!(
        Op::product_all([]).unwrap().trace_in(unit.space()).unwrap(),
        c(1.0, 0.0)
    );
}
