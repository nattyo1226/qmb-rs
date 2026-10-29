mod support;
use num_complex::Complex;
use qmb_core::{LocalOperator, LocalState, QuantumSystem};
use qmb_export::{Export, Exported, MatrixArgs};
use qmb_fermion::{
    Fermion, FermionOperator as O, FermionState as S, LadderBasis, MajoranaBasis, OccupationBasis,
};
use qmb_normalize::{Norm, Normalize, NormalizedExpectation};
use support::*;

fn vector<T: Real>(psi: &LocalState<S, O, Complex<T>>) -> Vec<Complex<f64>> {
    let out: Exported<Vec<Complex<T>>, _> = psi.export(&OccupationBasis).unwrap();
    assert_eq!(out.metadata().shape(), [2]);
    let array: Exported<qmb_export::ndarray::Array1<Complex<T>>, _> =
        psi.export(&OccupationBasis).unwrap();
    let sparse: Exported<qmb_export::sprs::CsVec<Complex<T>>, _> =
        psi.export(&OccupationBasis).unwrap();
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
fn exported<T: Real>(op: &LocalOperator<O, Complex<T>>) -> Matrix {
    let out: Exported<Vec<Complex<T>>, _> = op
        .export(MatrixArgs::<S, _>::new(&OccupationBasis))
        .unwrap();
    assert_eq!(out.metadata().shape(), [2, 2]);
    let array: Exported<qmb_export::ndarray::Array2<Complex<T>>, _> = op
        .export(MatrixArgs::<S, _>::new(&OccupationBasis))
        .unwrap();
    let sparse: Exported<qmb_export::sprs::CsMat<Complex<T>>, _> = op
        .export(MatrixArgs::<S, _>::new(&OccupationBasis))
        .unwrap();
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
    let state = Fermion::<Complex<T>>::state;
    let op = Fermion::<Complex<T>>::operator;
    let psi = state(S::Empty).scaled(c::<T>(2., 1.)) + state(S::Occupied).scaled(c::<T>(-1., 3.));
    let v = [Complex::new(2., 1.), Complex::new(-1., 3.)];
    values::<T>(&vector(&psi), &v);
    values::<T>(&vector(&psi.expand_in(&OccupationBasis).unwrap()), &v);
    close::<T>(psi.norm().to_f64().unwrap().into(), 15_f64.sqrt().into());
    close::<T>(psi.norm_squared().to_f64().unwrap().into(), 15_f64.into());
    close::<T>(
        psi.normalized().unwrap().norm().to_f64().unwrap().into(),
        1_f64.into(),
    );
    values::<T>(
        &vector(&psi.normalized().unwrap()),
        &v.map(|z| z / 15_f64.sqrt()),
    );
    for primitive in OPS {
        let a = op(primitive);
        let m = matrix(primitive);
        values::<T>(&exported(&a), &m);
        let applied = psi.act(&a);
        let expected = [m[0] * v[0] + m[1] * v[1], m[2] * v[0] + m[3] * v[1]];
        values::<T>(&vector(&applied.evaluate()), &expected);
        close::<T>(
            to64(psi.expectation(&a)),
            v[0].conj() * expected[0] + v[1].conj() * expected[1],
        );
        close::<T>(
            to64(psi.normalized_expectation(&a).unwrap()),
            (v[0].conj() * expected[0] + v[1].conj() * expected[1]) / 15.,
        );
        values::<T>(
            &exported(
                &a.expand_in(&LadderBasis)
                    .unwrap()
                    .expand_in(&MajoranaBasis)
                    .unwrap(),
            ),
            &m,
        );
    }
    // Complex coefficients ensure adjoints exercise scalar conjugation too.
    let a = op(O::Create).scaled(c::<T>(0.2, 0.7))
        + op(O::Parity).scaled(c::<T>(-0.4, 0.1))
        + op(O::Number);
    let reference = add(
        &add(
            &scale(&matrix(O::Create), Complex::new(0.2, 0.7)),
            &scale(&matrix(O::Parity), Complex::new(-0.4, 0.1)),
        ),
        &matrix(O::Number),
    );
    values::<T>(&exported(&a), &reference);
    values::<T>(&exported(&a.adjoint()), &adjoint(&reference));
    values::<T>(
        &exported(&a.adjoint().product(&a)),
        &product(&adjoint(&reference), &reference),
    );
    close::<T>(to64(a.trace()), trace(&reference));
    assert!(LocalState::<S, O, Complex<T>>::zero().normalized().is_err());
    assert!(
        LocalState::<S, O, Complex<T>>::zero()
            .normalized_expectation(&a)
            .is_err()
    );
}
#[test]
fn expressions_normalize_coordinates_and_all_exports_f32() {
    check::<f32>();
}
#[test]
fn expressions_normalize_coordinates_and_all_exports_f64() {
    check::<f64>();
}
