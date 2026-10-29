mod support;

use num_complex::Complex;
use qmb_core::QuantumSystem;
use qmb_export::{Export, Exported, MatrixArgs};
use qmb_normalize::{Norm, NormalizationError, Normalize, NormalizedExpectation};
use qmb_spin_one::{
    CartesianBasis, LadderBasis, SpinOne, SpinOneOperator as O, SpinOneState as S, SzBasis,
};
use support::*;

fn vector<T: Real>(psi: &LocalKet<T>) -> Matrix {
    let out: Exported<Vec<Complex<T>>, _> = psi.export(&SzBasis).unwrap();
    assert_eq!(out.metadata().shape(), [3]);
    let array: Exported<qmb_export::ndarray::Array1<Complex<T>>, _> = psi.export(&SzBasis).unwrap();
    let sparse: Exported<qmb_export::sprs::CsVec<Complex<T>>, _> = psi.export(&SzBasis).unwrap();
    let result: Vec<_> = out.data().iter().map(|x| to64(*x)).collect();
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
    result
}
fn exported<T: Real>(op: &LocalOp<T>) -> Matrix {
    let out: Exported<Vec<Complex<T>>, _> = op.export(MatrixArgs::<S, _>::new(&SzBasis)).unwrap();
    assert_eq!(out.metadata().shape(), [3, 3]);
    let array: Exported<qmb_export::ndarray::Array2<Complex<T>>, _> =
        op.export(MatrixArgs::<S, _>::new(&SzBasis)).unwrap();
    let sparse: Exported<qmb_export::sprs::CsMat<Complex<T>>, _> =
        op.export(MatrixArgs::<S, _>::new(&SzBasis)).unwrap();
    let result: Vec<_> = out.data().iter().map(|x| to64(*x)).collect();
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
    result
}
fn check<T: Real>() {
    let state = SpinOne::<Complex<T>>::state;
    let op = SpinOne::<Complex<T>>::operator;
    let psi = state(S::Plus).scaled(c::<T>(2., 1.))
        + state(S::Zero).scaled(c::<T>(-1., 3.))
        + state(S::Minus).scaled(c::<T>(0.5, -0.5));
    let v = [
        Complex::new(2., 1.),
        Complex::new(-1., 3.),
        Complex::new(0.5, -0.5),
    ];
    let n2 = dot(&v, &v).re;
    values::<T>(&vector(&psi), &v);
    values::<T>(&vector(&psi.expand_in(&SzBasis).unwrap()), &v);
    let coords = psi.coordinates_in(&SzBasis).unwrap();
    assert_eq!(coords.dimension().unwrap(), 3);
    for (i, z) in coords.flat_entries().map(Result::unwrap) {
        close::<T>(to64(*z), v[i]);
    }
    close::<T>(psi.norm_squared().to_f64().unwrap().into(), n2.into());
    close::<T>(psi.norm().to_f64().unwrap().into(), n2.sqrt().into());
    close::<T>(
        psi.normalized().unwrap().norm().to_f64().unwrap().into(),
        1.0.into(),
    );
    values::<T>(
        &vector(&psi.normalized().unwrap()),
        &scale(&v, (1. / n2.sqrt()).into()),
    );
    for primitive in OPS {
        let a = op(primitive);
        let m = matrix(primitive);
        values::<T>(&exported(&a), &m);
        let expected = apply(&m, &v);
        values::<T>(&vector(&psi.act(&a)), &expected);
        values::<T>(&vector(&psi.act(&a).evaluate()), &expected);
        close::<T>(to64(psi.expectation(&a)), dot(&v, &expected));
        close::<T>(
            to64(psi.normalized_expectation(&a).unwrap()),
            dot(&v, &expected) / n2,
        );
        for expanded in [
            a.expand_in(&CartesianBasis).unwrap(),
            a.expand_in(&LadderBasis).unwrap(),
        ] {
            values::<T>(&exported(&expanded), &m);
            values::<T>(&vector(&psi.act(&expanded)), &expected);
            close::<T>(to64(expanded.trace()), trace(&m));
        }
    }
    // Adjoint reverses products and conjugates genuinely complex coefficients.
    let a =
        op(O::SPlus).scaled(c::<T>(0.2, 0.7)) + op(O::Qyz).scaled(c::<T>(-0.4, 0.1)) + op(O::Q0);
    let reference = add(
        &add(
            &scale(&matrix(O::SPlus), Complex::new(0.2, 0.7)),
            &scale(&matrix(O::Qyz), Complex::new(-0.4, 0.1)),
        ),
        &matrix(O::Q0),
    );
    values::<T>(&exported(&a), &reference);
    values::<T>(&exported(&a.adjoint()), &adjoint(&reference));
    values::<T>(
        &exported(&a.adjoint().product(&a)),
        &product(&adjoint(&reference), &reference),
    );
    close::<T>(to64(a.trace()), trace(&reference));
    let duplicate = psi.clone();
    for zero in [LocalKet::<T>::zero(), &psi - &duplicate] {
        assert_eq!(zero.normalized().err(), Some(NormalizationError::ZeroNorm));
        assert_eq!(
            zero.normalized_expectation(&a).err(),
            Some(NormalizationError::ZeroNorm)
        );
    }
}
#[test]
fn local_expressions_and_all_backends_f32() {
    check::<f32>();
}
#[test]
fn local_expressions_and_all_backends_f64() {
    check::<f64>();
}
