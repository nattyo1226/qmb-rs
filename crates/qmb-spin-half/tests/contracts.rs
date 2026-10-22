//! Public contracts checked against independent two-dimensional linear algebra.

use std::{fmt::Debug, hash::Hash, ops::Neg};

use qmb_core::{
    Act, Adjoint, Basis, Complex, Identity, InnerProduct, LinearCombination as Lc, LocalStructure,
    Product, Trace,
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
        add(v, scale(matrix(t.operand), complex64(t.coefficient)))
    })
}
fn state_value<T: Real>(sum: &Lc<State, T>) -> Vector {
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
fn product<T: Real>(a: &Lc<Op, T>, b: &Lc<Op, T>) -> Lc<Op, T> {
    Lc::from_terms(a.terms().flat_map(|x| {
        b.terms().flat_map(move |y| {
            <Op as Product<T>>::product(&x.operand, &y.operand)
                .scaled(x.coefficient * y.coefficient)
        })
    }))
}
fn action<T: Real>(op: &Lc<Op, T>, state: &Lc<State, T>) -> Lc<State, T> {
    Lc::from_terms(op.terms().flat_map(|x| {
        state.terms().flat_map(move |y| {
            <State as Act<T>>::act(&y.operand, &x.operand).scaled(x.coefficient * y.coefficient)
        })
    }))
}
fn adjoint<T: Real>(op: &Lc<Op, T>) -> Lc<Op, T> {
    Lc::from_terms(
        op.terms()
            .flat_map(|t| <Op as Adjoint<T>>::adjoint(&t.operand).scaled(t.coefficient.conj())),
    )
}
fn trace<T: Real>(op: &Lc<Op, T>) -> Complex<T> {
    op.terms().fold(c::<T>(0.0, 0.0), |v, t| {
        v + t.coefficient * <Op as Trace<T>>::trace(&t.operand)
    })
}
fn inner<T: Real>(x: &Lc<State, T>, y: &Lc<State, T>) -> Complex<T> {
    let mut result = c::<T>(0.0, 0.0);
    for x in x {
        for y in y {
            result = result
                + x.coefficient.conj()
                    * y.coefficient
                    * <State as InnerProduct<T>>::inner_product(&x.operand, &y.operand);
        }
    }
    result
}
fn expand<T: Real, B: Basis<Op, T>>(basis: &B, sum: &Lc<Op, T>) -> Lc<Op, T> {
    Lc::from_terms(
        sum.terms()
            .flat_map(|t| basis.expand(&t.operand).scaled(t.coefficient)),
    )
}
fn coordinates<T: Real, B: Basis<Op, T>>(basis: &B, sum: &Lc<Op, T>) -> [Complex<T>; 4] {
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
    let identity = <SpinHalf as Identity<T>>::identity();
    assert_eq!(identity, Lc::primitive(I));
    for (i, a) in OPS.into_iter().enumerate() {
        for (j, b) in OPS.into_iter().enumerate() {
            let actual = <Op as Product<T>>::product(&a, &b);
            assert_eq!(
                actual.terms().map(|t| t.operand).collect::<Vec<_>>(),
                product_primitives[i][j]
            );
            assert_matrix::<T>(operator_value(&actual), multiply(matrix(a), matrix(b)));
        }
        for (j, state) in STATES.into_iter().enumerate() {
            let actual = <State as Act<T>>::act(&state, &a);
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
    let identity = <SpinHalf as Identity<T>>::identity();
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

fn check_operator_basis<T: Real, B: Basis<Op, T>>(basis: &B, expected: [Op; 4]) {
    let elements: Vec<_> = basis.elements().collect();
    assert_eq!(elements, expected);
    assert_eq!(elements.len(), operator_basis_len::<SpinHalf>().unwrap());
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
    fn state_basis<T: Real, B: Basis<State, T>>(basis: &B) {
        let elements: Vec<_> = basis.elements().collect();
        assert_eq!(elements, STATES);
        assert_eq!(elements.len(), state_basis_len::<SpinHalf>());
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
    fn metadata<L: LocalStructure<State = State, Operator = Op>>() {}
    metadata::<SpinHalf>();
    assert_eq!(<SpinHalf as LocalStructure>::DIMENSION.get(), 2);
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

// Expression evaluation is checked against the literal matrices above, not the
// primitive lifting helpers. Projections also exercise inner-product evaluation.
fn expression_vector<T: Real>(state: &qmb_core::LocalState<SpinHalf, T>) -> Vector {
    use qmb_core::LocalState;
    STATES.map(|basis| complex64(LocalState::<SpinHalf, T>::primitive(basis).inner_product(state)))
}
fn expression_matrix<T: Real>(op: &qmb_core::LocalOperator<SpinHalf, T>) -> Matrix {
    use qmb_core::LocalState;
    let columns = STATES.map(|basis| expression_vector(&LocalState::primitive(basis).act(op)));
    std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i]))
}
fn matrix_trace(matrix: Matrix) -> Complex<f64> {
    matrix[0][0] + matrix[1][1]
}

fn check_local_expressions<T: Real>() {
    use qmb_core::{LocalOperator, LocalState};
    type O<T> = LocalOperator<SpinHalf, T>;
    type S<T> = LocalState<SpinHalf, T>;
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
            let expanded = <CartesianBasis as Basis<Op, T>>::expand(&basis, &primitive);
            let expression = O::from_linear_combination(expanded.map(O::primitive));
            assert_matrix::<T>(expression_matrix(&expression), matrix(primitive));
        }
    }
    for primitive in OPS {
        let expanded = <LadderBasis as Basis<Op, T>>::expand(&LadderBasis, &primitive);
        assert_matrix::<T>(
            expression_matrix(&O::from_linear_combination(expanded.map(O::primitive))),
            matrix(primitive),
        );
    }
    for primitive in STATES {
        let expanded = <SzBasis as Basis<State, T>>::expand(&SzBasis, &primitive);
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
    type O = LocalOperator<SpinHalf, f64>;
    type S = LocalState<SpinHalf, f64>;
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
            type O = LocalOperator<SpinHalf, f64>;
            let mut expression = O::primitive(Sz);
            for _ in 0..20_000 {
                expression = expression
                    .product(&O::identity())
                    .adjoint()
                    .scaled(c(1.0, 0.0));
            }
            assert_complex::<f64>(expression.trace(), c(0.0, 0.0));
            drop(expression);
        })
        .unwrap()
        .join()
        .unwrap();
}
