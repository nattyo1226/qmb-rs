//! An independent three-level algebra using only the public qmb-core API.

use std::num::NonZeroUsize;

use num_complex::Complex;
use qmb_core::{
    Act, Adjoint, Basis, Dimension, DimensionOverflow, Identity, InnerProduct, LinearCombination,
    Product, Trace, dim_op, dim_state,
};

#[path = "support/three_level.rs"]
mod fixture;
use fixture::*;

fn check_rules<T: Real>() {
    let operators: Vec<_> = <MatrixUnits as Basis<Operator, Complex<T>>>::elements(&MatrixUnits)
        .chain([Operator::Identity])
        .collect();
    let identity = LinearCombination::<_, Complex<T>>::primitive(Operator::identity());
    assert_eq!(identity.len(), 1);
    assert_matrix(operator_value(&identity), matrix(Operator::Identity));
    for &a in &operators {
        let a_sum = LinearCombination::<_, Complex<T>>::primitive(a);
        let a_matrix = matrix::<T>(a);
        let adjoint = <Operator as Adjoint<Complex<T>>>::adjoint(&a);
        let expected_adjoint =
            std::array::from_fn(|i| std::array::from_fn(|j| a_matrix[j][i].conj()));
        assert_matrix(operator_value(&adjoint), expected_adjoint);
        assert_matrix(operator_value(&product(&identity, &a_sum)), a_matrix);
        assert_matrix(operator_value(&product(&a_sum, &identity)), a_matrix);
        let expected_trace = (0..3).fold(c(0.0, 0.0), |sum, i| sum + a_matrix[i][i]);
        assert_complex(<Operator as Trace<Complex<T>>>::trace(&a), expected_trace);
        for state in KETS {
            let applied = act_primitive::<Ket, Operator, T>(&state, &a);
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
            let ab = <Operator as Product<Complex<T>>>::product(&a, &b);
            assert_matrix(operator_value(&ab), matmul(a_matrix, matrix(b)));
            for state in KETS {
                assert_vector(
                    action(&ab, &state),
                    matvec(
                        a_matrix,
                        state_value(&<Ket as Act<Complex<T>>>::act(&state, &b)),
                    ),
                );
            }
            for &d in &operators {
                let bc = <Operator as Product<Complex<T>>>::product(&b, &d);
                let d_sum = LinearCombination::primitive(d);
                assert_matrix(
                    operator_value(&product(&ab, &d_sum)),
                    operator_value(&product(&a_sum, &bc)),
                );
            }
            // Cyclicity of the ordinary trace.
            let trace = |sum: &LinearCombination<Operator, Complex<T>>| {
                sum.terms().fold(c::<T>(0.0, 0.0), |v, t| {
                    v + t.coefficient * <Operator as Trace<Complex<T>>>::trace(&t.operand)
                })
            };
            assert_complex(
                trace(&ab),
                trace(&<Operator as Product<Complex<T>>>::product(&b, &a)),
            );
        }
    }
    let x = LinearCombination::<_, Complex<T>>::primitive(Ket::Sum).scaled(c::<T>(0.3, -0.7))
        + LinearCombination::<_, Complex<T>>::primitive(Ket::Phase).scaled(c::<T>(1.2, 0.4));
    let y =
        LinearCombination::<_, Complex<T>>::primitive(Ket::Unit(Axis::B)).scaled(c::<T>(-0.2, 0.6));
    assert_complex(inner(&x, &y), dot(state_value(&x), state_value(&y)));
    assert_complex(inner(&x, &y), inner(&y, &x).conj());
    assert_complex(
        <Ket as InnerProduct<Complex<T>>>::inner_product(&Ket::Sum, &Ket::Phase),
        c(1.0, 1.0),
    );
    let dependent = LinearCombination::primitive(Ket::Sum)
        + LinearCombination::<_, Complex<T>>::primitive(Ket::Unit(Axis::A))
            .scaled(c::<T>(-1.0, 0.0))
        + LinearCombination::<_, Complex<T>>::primitive(Ket::Unit(Axis::B))
            .scaled(c::<T>(-1.0, 0.0));
    assert_eq!(dependent.len(), 3);
    assert_vector(state_value(&dependent), [c(0.0, 0.0); 3]);
    assert_complex(inner(&dependent, &dependent), c(0.0, 0.0));
}

fn check_state_basis<T: Real, B: Basis<Ket, Complex<T>>>(basis: &B) {
    let elements: Vec<_> = basis.elements().collect();
    assert_eq!(elements.len(), dim_state::<Ket>().unwrap());
    for (i, p) in elements.iter().enumerate() {
        assert_eq!(basis.index(p), Some(i));
    }
    for p in KETS {
        let expanded = basis.expand(&p);
        assert_vector(state_value(&expanded), p.coordinates());
        for term in &expanded {
            assert!(basis.index(&term.operand).is_some());
        }
    }
}

fn check_bases<T: Real>() {
    check_state_basis::<T, _>(&Canonical);
    check_state_basis::<T, _>(&NonOrthogonal);
    assert_eq!(
        <Canonical as Basis<Ket, Complex<T>>>::index(&Canonical, &Ket::Sum),
        None
    );
    assert_eq!(
        <NonOrthogonal as Basis<Ket, Complex<T>>>::index(&NonOrthogonal, &Ket::Unit(Axis::B)),
        None
    );
    let expanded =
        <NonOrthogonal as Basis<Ket, Complex<T>>>::expand(&NonOrthogonal, &Ket::Unit(Axis::B));
    assert_complex(inner(&expanded, &expanded), c(1.0, 0.0));
    // Squared coefficient magnitudes alone give 2, while the actual norm² is 1.
    let coefficient_squares = expanded
        .terms()
        .fold(T::zero(), |v, t| v + t.coefficient.norm_sqr());
    assert_complex(Complex::new(coefficient_squares, T::zero()), c(2.0, 0.0));
    fn require_operator_basis<T: Real, B: Basis<Operator, Complex<T>>>(_: &B) {}
    require_operator_basis::<T, _>(&MatrixUnits);
    let elements: Vec<_> =
        <MatrixUnits as Basis<Operator, Complex<T>>>::elements(&MatrixUnits).collect();
    assert_eq!(elements.len(), dim_op::<Operator>().unwrap());
    for (i, p) in elements.iter().enumerate() {
        assert_eq!(
            <MatrixUnits as Basis<Operator, Complex<T>>>::index(&MatrixUnits, p),
            Some(i)
        );
        assert_matrix(
            operator_value(&<MatrixUnits as Basis<Operator, Complex<T>>>::expand(
                &MatrixUnits,
                p,
            )),
            matrix(*p),
        );
    }
    assert_eq!(
        <MatrixUnits as Basis<Operator, Complex<T>>>::index(&MatrixUnits, &Operator::Identity),
        None
    );
    assert_matrix(
        operator_value(&<MatrixUnits as Basis<Operator, Complex<T>>>::expand(
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
    assert_eq!(Ket::dim().unwrap().get(), 3);
}

struct FixedDimension<const D: usize>;
struct SquaredDimension<const D: usize>;
impl<const D: usize> Dimension for FixedDimension<D> {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(D).unwrap())
    }
}
impl<const D: usize> Dimension for SquaredDimension<D> {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        D.checked_mul(D)
            .and_then(NonZeroUsize::new)
            .ok_or(DimensionOverflow)
    }
}

#[test]
fn basis_dimensions_handle_unit_large_and_overflow_cases() {
    assert_eq!(dim_state::<FixedDimension<1>>().unwrap(), 1);
    assert_eq!(dim_op::<SquaredDimension<1>>(), Ok(1));
    assert_eq!(dim_state::<Ket>().unwrap(), 3);
    assert_eq!(dim_op::<Operator>(), Ok(9));
    const FITS: usize = 1 << (usize::BITS / 2 - 1);
    const OVERFLOWS: usize = 1 << (usize::BITS / 2);
    assert_eq!(dim_op::<SquaredDimension<FITS>>(), Ok(FITS * FITS));
    assert_eq!(
        dim_op::<SquaredDimension<OVERFLOWS>>(),
        Err(DimensionOverflow)
    );
    assert_eq!(
        dim_state::<FixedDimension<{ usize::MAX }>>().unwrap(),
        usize::MAX
    );
    assert_eq!(
        dim_op::<SquaredDimension<{ usize::MAX }>>(),
        Err(DimensionOverflow)
    );
    let error: &dyn core::error::Error = &DimensionOverflow;
    assert_eq!(error.to_string(), "space dimension overflows usize");
    assert!(error.source().is_none());
}

fn expression_coordinates<T: Real>(
    state: &qmb_core::LocalState<Ket, Operator, Complex<T>>,
) -> Vector<T> {
    AXES.map(|axis| {
        qmb_core::LocalState::<Ket, Operator, Complex<T>>::primitive(Ket::Unit(axis))
            .inner_product(state)
    })
}
fn expression_operator_matrix<T: Real>(
    operator: &qmb_core::LocalOperator<Operator, Complex<T>>,
) -> Matrix<T> {
    let columns = AXES.map(|axis| {
        expression_coordinates(&qmb_core::LocalState::primitive(Ket::Unit(axis)).act(operator))
    });
    std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i]))
}
fn check_three_level_expressions<T: Real>() {
    use qmb_core::{LocalOperator, LocalState};
    type S<T> = LocalState<Ket, Operator, Complex<T>>;
    type O<T> = LocalOperator<Operator, Complex<T>>;
    let state = S::<T>::primitive(Ket::Sum).scaled(c(0.2, 0.3))
        + S::primitive(Ket::Phase).scaled(c(-0.4, 0.1));
    let reference = state_value(
        &(LinearCombination::term(c::<T>(0.2, 0.3), Ket::Sum)
            + LinearCombination::term(c(-0.4, 0.1), Ket::Phase)),
    );
    assert_vector(expression_coordinates(&state), reference);
    assert_complex(state.inner_product(&state), dot(reference, reference));
    let operators: Vec<_> = <MatrixUnits as Basis<Operator, Complex<T>>>::elements(&MatrixUnits)
        .chain([Operator::Identity])
        .collect();
    for a in &operators {
        let a_expr = O::<T>::primitive(*a).scaled(c(0.3, 0.2)) + O::identity();
        let mut a_matrix = matrix::<T>(*a);
        for (i, row) in a_matrix.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = *value * c(0.3, 0.2) + if i == j { c(1.0, 0.0) } else { c(0.0, 0.0) };
            }
        }
        for b in &operators {
            let product = a_expr.clone() * &O::primitive(*b);
            let expected = matmul(a_matrix, matrix(*b));
            let expected_adjoint =
                std::array::from_fn(|i| std::array::from_fn(|j| expected[j][i].conj()));
            assert_matrix(expression_operator_matrix(&product), expected);
            assert_matrix(
                expression_operator_matrix(&product.adjoint()),
                expected_adjoint,
            );
            assert_complex(
                product.trace(),
                (0..3).fold(c(0.0, 0.0), |sum, i| sum + expected[i][i]),
            );
            assert_vector(
                expression_coordinates(&state.act(&product).evaluate()),
                matvec(expected, reference),
            );
            assert_complex(
                state.expectation(&product),
                dot(reference, matvec(expected, reference)),
            );
            let a_adjoint = std::array::from_fn(|i| std::array::from_fn(|j| a_matrix[j][i].conj()));
            let hs = matmul(a_adjoint, matrix(*b));
            assert_complex(
                a_expr.hs_inner_product(&O::primitive(*b)),
                (0..3).fold(c::<T>(0.0, 0.0), |sum, i| sum + hs[i][i]) / c::<T>(3.0, 0.0),
            );
        }
    }
    // The nonorthogonal basis uses a dependent primitive; importing its terms
    // must preserve overlaps, rather than interpreting coefficients as coordinates.
    for primitive in KETS {
        let terms = <NonOrthogonal as Basis<Ket, Complex<T>>>::expand(&NonOrthogonal, &primitive);
        let imported = S::from_terms(terms.map(S::primitive));
        assert_vector(expression_coordinates(&imported), primitive.coordinates());
        assert_complex(
            imported.inner_product(&imported),
            dot(primitive.coordinates(), primitive.coordinates()),
        );
    }
    let identity = O::<T>::primitive(Operator::identity());
    assert_matrix(
        expression_operator_matrix(&identity),
        matrix(Operator::Identity),
    );
    assert_complex(O::<T>::identity().trace(), c(3.0, 0.0));
    let basis_operator = O::from_terms(
        <MatrixUnits as Basis<Operator, Complex<T>>>::expand(&MatrixUnits, &Operator::Identity)
            .map(O::primitive),
    );
    assert_matrix(
        expression_operator_matrix(&basis_operator),
        matrix(Operator::Identity),
    );
    assert_vector(expression_coordinates(&S::<T>::zero()), [c(0.0, 0.0); 3]);
    assert_complex(O::<T>::zero().trace(), c(0.0, 0.0));
}

#[test]
fn three_level_expressions_and_nonorthogonal_basis_match_matrices() {
    check_three_level_expressions::<f32>();
    check_three_level_expressions::<f64>();
}

type CompositeState<T> = qmb_core::State<u8, Ket, Operator, Complex<T>>;
type CompositeOperator<T> = qmb_core::Operator<u8, Operator, Complex<T>>;
type CompositeSpace = qmb_core::HilbertSpace<u8>;
type PairVector<T> = [Complex<T>; 9];
type PairMatrix<T> = [[Complex<T>; 9]; 9];

fn kronecker<T: Real>(a: Matrix<T>, b: Matrix<T>) -> PairMatrix<T> {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i / 3][j / 3] * b[i % 3][j % 3]))
}

fn tensor_vector<T: Real>(a: Vector<T>, b: Vector<T>) -> PairVector<T> {
    std::array::from_fn(|i| a[i / 3] * b[i % 3])
}

fn pair_multiply<T: Real>(a: PairMatrix<T>, b: PairMatrix<T>) -> PairMatrix<T> {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..9).fold(c(0.0, 0.0), |sum, k| sum + a[i][k] * b[k][j]))
    })
}

fn pair_adjoint<T: Real>(a: PairMatrix<T>) -> PairMatrix<T> {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i].conj()))
}

fn pair_apply<T: Real>(a: PairMatrix<T>, x: PairVector<T>) -> PairVector<T> {
    a.map(|row| {
        row.into_iter()
            .zip(x)
            .fold(c(0.0, 0.0), |sum, (a, x)| sum + a * x)
    })
}

fn pair_dot<T: Real>(a: PairVector<T>, b: PairVector<T>) -> Complex<T> {
    a.into_iter()
        .zip(b)
        .fold(c(0.0, 0.0), |sum, (a, b)| sum + a.conj() * b)
}

fn pair_trace<T: Real>(a: PairMatrix<T>) -> Complex<T> {
    (0..9).fold(c(0.0, 0.0), |sum, i| sum + a[i][i])
}

fn composite_coordinates<T: Real>(
    state: &CompositeState<T>,
    factors: &[CompositeSpace; 2],
) -> PairVector<T> {
    use qmb_core::LocalState;
    let evaluated = state.evaluate().unwrap();
    std::array::from_fn(|n| {
        let left = LocalState::primitive(Ket::Unit(AXES[n / 3]))
            .at(&factors[0])
            .unwrap();
        let right = LocalState::primitive(Ket::Unit(AXES[n % 3]))
            .at(&factors[1])
            .unwrap();
        qmb_core::State::tensor_in(
            &CompositeSpace::tensor_all(factors.iter().cloned()).unwrap(),
            [left, right],
        )
        .unwrap()
        .inner_product(&evaluated)
        .unwrap()
    })
}

fn assert_pair_vector<T: Real>(actual: PairVector<T>, expected: PairVector<T>) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert_complex(a, b);
    }
}

fn check_three_level_composites<T: Real>() {
    use qmb_core::{LocalOperator as O, LocalState as S};
    let factors = [CompositeSpace::factor(17), CompositeSpace::factor(3)];
    let [a, b] = &factors;
    let h = a.tensor(b).unwrap();
    let grouped = a
        .tensor(&CompositeSpace::tensor_all([CompositeSpace::unit(), b.clone()]).unwrap())
        .unwrap();
    assert!(h.is_compatible_with(&grouped));

    // Import a nonorthogonal expansion, then place the whole local expression.
    // The reference uses literal ket coordinates, not basis expansion or Act.
    let local = S::<Ket, Operator, Complex<T>>::from_terms(
        <NonOrthogonal as Basis<Ket, Complex<T>>>::expand(&NonOrthogonal, &Ket::Phase)
            .map(S::primitive),
    );
    let left = local.at(a).unwrap();
    let right = S::primitive(Ket::Sum).at(b).unwrap();
    let first = qmb_core::State::tensor_in(
        &CompositeSpace::tensor_all(factors.iter().cloned()).unwrap(),
        [left, right],
    )
    .unwrap();
    let second = qmb_core::State::tensor_in(
        &h,
        [
            S::primitive(Ket::Unit(Axis::C)).at(a).unwrap(),
            S::primitive(Ket::Phase).at(b).unwrap(),
        ],
    )
    .unwrap();
    let state = (first.clone() + &second.scaled(c(-0.2, 0.3))).unwrap();
    let x = tensor_vector(Ket::Phase.coordinates::<T>(), Ket::Sum.coordinates());
    let y = tensor_vector(Ket::Unit(Axis::C).coordinates(), Ket::Phase.coordinates());
    let reference = std::array::from_fn(|i| x[i] + c::<T>(-0.2, 0.3) * y[i]);
    assert_pair_vector(composite_coordinates(&state, &factors), reference);
    assert_complex(
        state.inner_product(&state).unwrap(),
        pair_dot(reference, reference),
    );

    let primitives: Vec<_> = AXES
        .into_iter()
        .flat_map(|i| AXES.map(|j| Operator::Unit(i, j)))
        .chain([Operator::Identity])
        .collect();
    let identity = matrix(Operator::Identity);
    for &oa in &primitives {
        let alpha = (O::<Operator, Complex<T>>::primitive(oa)
            .at(a)
            .unwrap()
            .scaled(c(0.3, -0.4))
            .clone()
            + &CompositeOperator::scalar(c(-0.1, 0.2)))
            .unwrap();
        let am: PairMatrix<T> = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                c::<T>(0.3, -0.4) * kronecker(matrix(oa), identity)[i][j]
                    + c::<T>(-0.1, 0.2) * kronecker(identity, identity)[i][j]
            })
        });
        for &ob in &primitives {
            let beta = (O::primitive(ob).at(b).unwrap().clone()
                + &O::primitive(Operator::Unit(Axis::B, Axis::A))
                    .at(a)
                    .unwrap()
                    .scaled(c(0.2, 0.1)))
                .unwrap();
            let bm: PairMatrix<T> = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    kronecker(identity, matrix(ob))[i][j]
                        + c::<T>(0.2, 0.1)
                            * kronecker(matrix(Operator::Unit(Axis::B, Axis::A)), identity)[i][j]
                })
            });
            let product = (alpha.clone() * &beta).unwrap();
            let pm = pair_multiply(am, bm);
            let applied = state.act(&product).unwrap();
            assert_pair_vector(
                composite_coordinates(&applied, &factors),
                pair_apply(pm, reference),
            );
            assert_pair_vector(
                composite_coordinates(&state.act(&beta).unwrap().act(&alpha).unwrap(), &factors),
                pair_apply(pm, reference),
            );
            let adjoint = product.adjoint();
            assert_pair_vector(
                composite_coordinates(&state.act(&adjoint).unwrap(), &factors),
                pair_apply(pair_adjoint(pm), reference),
            );
            assert_complex(product.trace_in(&grouped).unwrap(), pair_trace(pm));
            assert_complex(adjoint.trace_in(&h).unwrap(), pair_trace(pair_adjoint(pm)));
            assert_complex(
                state.expectation(&product).unwrap(),
                pair_dot(reference, pair_apply(pm, reference)),
            );
            assert_complex(
                alpha.hs_inner_product_in(&beta, &h).unwrap(),
                pair_trace(pair_multiply(pair_adjoint(am), bm)) / c(9.0, 0.0),
            );
        }
    }
    assert_pair_vector(composite_coordinates(&state, &factors), reference);
    // Identity is returned as three matrix-unit terms by this algebra.
    assert_complex(
        CompositeOperator::<T>::identity().trace_in(&h).unwrap(),
        c(9.0, 0.0),
    );
}

#[test]
fn three_level_composites_and_nonorthogonal_imports_match_kronecker_matrices() {
    check_three_level_composites::<f32>();
    check_three_level_composites::<f64>();
}

fn check_inner_products_and_squared_norms<T: Real>() {
    use qmb_core::{LocalOperator as O, LocalState as S};
    let x = S::<Ket, Operator, Complex<T>>::primitive(Ket::Sum);
    let y = S::<Ket, Operator, Complex<T>>::primitive(Ket::Phase);
    let z = S::<Ket, Operator, Complex<T>>::primitive(Ket::Unit(Axis::A));
    let alpha = c::<T>(2.0, 1.0);
    let beta = c::<T>(-1.0, 2.0);
    let conjugate = |v: Complex<T>| Complex::new(v.re, -v.im);
    assert_complex(Complex::new(x.norm_squared(), T::zero()), c(2.0, 0.0));

    // Overlapping primitives and a nonorthogonal expansion require actual
    // inner products, rather than squared coefficient magnitudes.
    let overlapping = &x + &y;
    assert_complex(
        Complex::new(overlapping.norm_squared(), T::zero()),
        c(6.0, 0.0),
    );
    let imported = S::<Ket, Operator, Complex<T>>::from_terms(
        <NonOrthogonal as Basis<Ket, Complex<T>>>::expand(&NonOrthogonal, &Ket::Unit(Axis::B))
            .map(S::primitive),
    );
    assert_complex(
        Complex::new(imported.norm_squared(), T::zero()),
        c(1.0, 0.0),
    );
    assert_complex(
        (&x.scaled(alpha) + &y.scaled(beta)).inner_product(&z),
        conjugate(alpha) * x.inner_product(&z) + conjugate(beta) * y.inner_product(&z),
    );
    assert_complex(
        z.inner_product(&(&x.scaled(alpha) + &y.scaled(beta))),
        alpha * z.inner_product(&x) + beta * z.inner_product(&y),
    );
    let op = O::primitive(Operator::Unit(Axis::B, Axis::A));
    assert_complex(
        x.inner_product(&y.act(&op)),
        dot(
            Ket::Sum.coordinates(),
            matvec(
                matrix(Operator::Unit(Axis::B, Axis::A)),
                Ket::Phase.coordinates(),
            ),
        ),
    );
    for zero in [S::zero(), &x - &x.clone()] {
        assert_complex(Complex::new(zero.norm_squared(), T::zero()), c(0.0, 0.0));
    }
    let a = CompositeSpace::factor(11);
    let b = CompositeSpace::factor(7);
    let h = a.tensor(&b).unwrap();
    let first = qmb_core::State::tensor_in(&h, [x.at(&a).unwrap(), y.at(&b).unwrap()]).unwrap();
    let second = qmb_core::State::tensor_in(&h, [y.at(&a).unwrap(), z.at(&b).unwrap()]).unwrap();
    assert_complex(
        Complex::new(first.norm_squared().unwrap(), T::zero()),
        c(4.0, 0.0),
    );

    let sum = (first.scaled(alpha).clone() + &second.scaled(beta)).unwrap();
    assert_complex(
        sum.inner_product(&first).unwrap(),
        conjugate(alpha) * first.inner_product(&first).unwrap()
            + conjugate(beta) * second.inner_product(&first).unwrap(),
    );
    assert_complex(
        first.inner_product(&sum).unwrap(),
        alpha * first.inner_product(&first).unwrap() + beta * first.inner_product(&second).unwrap(),
    );
    for zero in [
        qmb_core::State::zero(&h),
        (&first - &first.clone()).unwrap(),
    ] {
        assert_complex(
            Complex::new(zero.norm_squared().unwrap(), T::zero()),
            c(0.0, 0.0),
        );
    }
}

#[test]
fn inner_products_and_squared_norms_work_at_both_precisions_without_sqrt() {
    check_inner_products_and_squared_norms::<f32>();
    check_inner_products_and_squared_norms::<f64>();
}
