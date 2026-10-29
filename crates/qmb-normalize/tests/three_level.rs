use num_complex::Complex;
use num_traits::Float;
use qmb_core::{Basis, HilbertSpace};
use qmb_normalize::{
    Norm, NormalizationError, Normalize, NormalizedExpectation, NormalizedExpectationIn,
};
#[path = "../../qmb-core/tests/support/three_level.rs"]
mod fixture;
use fixture::*;
type CompositeSpace = HilbertSpace<u8>;
type CompositeOperator<T> = qmb_core::Operator<u8, Operator, Complex<T>>;

fn check_three_level_normalization<T: Real + Float>() {
    use qmb_core::LocalState as S;
    let x = S::<Ket, Operator, Complex<T>>::primitive(Ket::Sum).scaled(c(2.0, -1.0));
    let normalized = x.normalized().unwrap();
    assert_complex(Complex::new(normalized.norm(), T::zero()), c(1.0, 0.0));
    let a = CompositeSpace::factor(1);
    let state = x.at(&a).unwrap();
    assert_complex(
        Complex::new(state.normalized().unwrap().norm().unwrap(), T::zero()),
        c(1.0, 0.0),
    );
    for zero in [S::zero(), &x - &x.clone()] {
        assert_complex(Complex::new(zero.norm(), T::zero()), c(0.0, 0.0));
        assert_eq!(zero.normalized().err(), Some(NormalizationError::ZeroNorm));
        let zero = zero.at(&a).unwrap();
        assert_complex(Complex::new(zero.norm().unwrap(), T::zero()), c(0.0, 0.0));
        assert_eq!(zero.normalized().err(), Some(NormalizationError::ZeroNorm));
    }
}
#[test]
fn three_level_normalization_and_zero_norm_rejection_both_precisions() {
    check_three_level_normalization::<f32>();
    check_three_level_normalization::<f64>();
}

fn check_derived_state_operations<T: Real>() {
    use qmb_core::{ExpressionError, LocalOperator as O, LocalState as S};
    let x = S::<Ket, Operator, Complex<T>>::primitive(Ket::Sum);
    let y = S::<Ket, Operator, Complex<T>>::primitive(Ket::Phase);
    let z = S::<Ket, Operator, Complex<T>>::primitive(Ket::Unit(Axis::A));
    let alpha = c::<T>(2.0, 1.0);
    let beta = c::<T>(-1.0, 2.0);
    let conjugate = |v: Complex<T>| Complex::new(v.re, -v.im);
    let identity = O::identity();
    assert_complex(Complex::new(x.norm_squared(), T::zero()), c(2.0, 0.0));
    assert_complex(
        x.scaled(alpha).normalized_expectation(&identity).unwrap(),
        c(1.0, 0.0),
    );
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
        assert_eq!(
            zero.normalized_expectation(&identity),
            Err(NormalizationError::ZeroNorm)
        );
    }
    let a = CompositeSpace::factor(11);
    let b = CompositeSpace::factor(7);
    let h = a.tensor(&b).unwrap();
    let first = qmb_core::State::tensor_in(&h, [x.at(&a).unwrap(), y.at(&b).unwrap()]).unwrap();
    let second = qmb_core::State::tensor_in(&h, [y.at(&a).unwrap(), z.at(&b).unwrap()]).unwrap();
    let ci = CompositeOperator::identity().bind(&h).unwrap();
    assert_complex(
        Complex::new(first.norm_squared().unwrap(), T::zero()),
        c(4.0, 0.0),
    );
    assert_complex(
        first.scaled(alpha).normalized_expectation(&ci).unwrap(),
        c(1.0, 0.0),
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
        assert_eq!(
            zero.normalized_expectation(&ci),
            Err(NormalizationError::ZeroNorm)
        );
        assert_eq!(
            zero.normalized_expectation_in(&ci, &h),
            Err(NormalizationError::ZeroNorm)
        );
    }
    // An explicit ambient space overrides the incompatible stored binding.
    let stale = CompositeOperator::identity().bind(&a).unwrap();
    assert!(first.normalized_expectation(&stale).is_err());
    assert_complex(
        first.normalized_expectation_in(&stale, &h).unwrap(),
        c(1.0, 0.0),
    );
    assert_eq!(
        first.normalized_expectation_in(&ci, &a),
        Err(NormalizationError::Expression(
            ExpressionError::IncompatibleSpaces
        ))
    );
}

#[test]
fn derived_states_use_inner_products_for_both_precisions_without_sqrt() {
    check_derived_state_operations::<f32>();
    check_derived_state_operations::<f64>();
}
