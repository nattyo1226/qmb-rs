mod support;

use num_complex::Complex;
use qmb_core::{Basis, HilbertSpace, QuantumSystem, TensorKind};
use qmb_export::{Export, Exported, MatrixArgs};
use qmb_normalize::{Norm, Normalize, NormalizedExpectation};
use qmb_spin_one::{
    CartesianBasis, LadderBasis, SpinOne, SpinOneOperator as O, SpinOneState as S, SzBasis,
};
use support::*;

type Op<T> = qmb_core::Operator<u8, O, Complex<T>>;
type Ket<T> = qmb_core::State<u8, S, O, Complex<T>>;

fn vector<T: Real>(psi: &Ket<T>) -> Matrix {
    let out: Exported<Vec<Complex<T>>, _> = psi.export(&SzBasis).unwrap();
    assert_eq!(out.metadata().shape(), [9]);
    assert_eq!(out.metadata().factor_count(), 2);
    assert_eq!(out.metadata().local_dimension(), Some(3));
    out.data().iter().map(|z| to64(*z)).collect()
}
fn exported<T: Real>(op: &Op<T>) -> Matrix {
    let args = || MatrixArgs::<S, _>::new(&SzBasis);
    let dense: Exported<Vec<Complex<T>>, _> = op.export(args()).unwrap();
    let array: Exported<qmb_export::ndarray::Array2<Complex<T>>, _> = op.export(args()).unwrap();
    let sparse: Exported<qmb_export::sprs::CsMat<Complex<T>>, _> = op.export(args()).unwrap();
    assert_eq!(dense.metadata().shape(), [9, 9]);
    assert_eq!(dense.metadata().factor_count(), 2);
    assert_eq!(dense.metadata().local_dimension(), Some(3));
    assert!(
        dense
            .metadata()
            .space()
            .is_compatible_with(op.space().unwrap())
    );
    let result: Vec<_> = dense.data().iter().map(|z| to64(*z)).collect();
    values::<T>(
        &array.data().iter().map(|z| to64(*z)).collect::<Vec<_>>(),
        &result,
    );
    values::<T>(
        &sparse
            .data()
            .to_dense()
            .iter()
            .map(|z| to64(*z))
            .collect::<Vec<_>>(),
        &result,
    );
    result
}

fn basis<T: Real, B: Basis<O, Complex<T>>>(
    space: &HilbertSpace<u8>,
    states: &[Ket<T>],
    basis: &B,
    elements: [O; 9],
) {
    let operators: Vec<_> = space
        .basis_operators::<O, Complex<T>, _>(basis)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(operators.len(), 81);
    for (i, op) in operators.into_iter().enumerate() {
        let expected = kron(&matrix(elements[i / 9]), &matrix(elements[i % 9]));
        values::<T>(&exported(&op), &expected);
        close::<T>(to64(op.trace().unwrap()), trace(&expected));
        let coordinates = op.coordinates_in(basis).unwrap();
        assert_eq!(coordinates.dimension().unwrap(), 81);
        assert_eq!(coordinates.len(), 1);
        let (index, coefficient) = coordinates.flat_entries().next().unwrap().unwrap();
        assert_eq!(index, i);
        close::<T>(to64(*coefficient), 1.0.into());
        for (col, psi) in states.iter().enumerate() {
            let actual = psi.act(&op).unwrap().evaluate().unwrap();
            let column: Vec<_> = (0..9).map(|row| expected[row * 9 + col]).collect();
            values::<T>(&vector(&actual), &column);
        }
    }
}

fn check<T: Real>() {
    let local = SpinOne::<Complex<T>>::operator;
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        // Numeric index order is reversed deliberately: factor order is explicit.
        let a = HilbertSpace::factor(11u8);
        let b = HilbertSpace::factor(7u8);
        let space = a.tensor_product_with_kind(kind, &b).unwrap();
        assert_eq!(space.dimension::<S>().unwrap(), 9);
        assert_eq!(space.dimension::<O>().unwrap(), 81);
        let states: Vec<_> = space
            .basis_states::<S, O, Complex<T>, _>(&SzBasis)
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(states.len(), 9);
        for (i, psi) in states.iter().enumerate() {
            let expected: Vec<_> = (0..9)
                .map(|j| Complex::new(f64::from(i == j), 0.))
                .collect();
            values::<T>(&vector(psi), &expected);
            let coordinates = psi.coordinates_in(&SzBasis).unwrap();
            assert_eq!(coordinates.dimension().unwrap(), 9);
            let (index, coefficient) = coordinates.flat_entries().next().unwrap().unwrap();
            assert_eq!(index, i);
            close::<T>(to64(*coefficient), 1.0.into());
        }
        if kind == TensorKind::Ordinary {
            basis::<T, _>(&space, &states, &CartesianBasis, CART);
            basis::<T, _>(&space, &states, &LadderBasis, LADDER);
        }
        let identity = Op::<T>::identity().bind(&space).unwrap();
        let im = kron(&matrix(O::I), &matrix(O::I));
        values::<T>(&exported(&identity), &im);
        close::<T>(to64(identity.trace().unwrap()), 9.0.into());
        for primitive in OPS {
            for (factor, left) in [(&a, true), (&b, false)] {
                let placed = local(primitive).at(factor).unwrap().bind(&space).unwrap();
                let expected = if left {
                    kron(&matrix(primitive), &matrix(O::I))
                } else {
                    kron(&matrix(O::I), &matrix(primitive))
                };
                values::<T>(&exported(&placed), &expected);
                close::<T>(to64(placed.trace().unwrap()), trace(&expected));
            }
        }
        let first = local(O::SPlus)
            .at(&a)
            .unwrap()
            .bind(&space)
            .unwrap()
            .add(&local(O::Qyz).at(&b).unwrap().bind(&space).unwrap())
            .unwrap();
        let second = SpinOne::<Complex<T>>::operator_product_at([(&a, O::Qzx), (&b, O::Sy)])
            .unwrap()
            .bind(&space)
            .unwrap();
        let op = first
            .adjoint()
            .product(&second)
            .unwrap()
            .add(&identity.scaled(c::<T>(0.2, -0.3)))
            .unwrap();
        let first_m = add(
            &kron(&matrix(O::SPlus), &matrix(O::I)),
            &kron(&matrix(O::I), &matrix(O::Qyz)),
        );
        let second_m = kron(&matrix(O::Qzx), &matrix(O::Sy));
        let expected = add(
            &product(&adjoint(&first_m), &second_m),
            &scale(&im, Complex::new(0.2, -0.3)),
        );
        values::<T>(&exported(&op), &expected);
        values::<T>(&exported(&op.adjoint()), &adjoint(&expected));
        close::<T>(to64(op.trace().unwrap()), trace(&expected));
        for expanded in [
            op.expand_in(&CartesianBasis).unwrap(),
            op.expand_in(&LadderBasis).unwrap(),
        ] {
            values::<T>(&exported(&expanded), &expected);
            assert!(expanded.space().unwrap().is_compatible_with(&space));
        }
        let psi = states[0]
            .scaled(c::<T>(1., 2.))
            .add(&states[4].scaled(c::<T>(0.3, -0.7)))
            .unwrap()
            .add(&states[6].scaled(c::<T>(-1., 0.4)))
            .unwrap();
        let mut v = vec![0.0.into(); 9];
        v[0] = Complex::new(1., 2.);
        v[4] = Complex::new(0.3, -0.7);
        v[6] = Complex::new(-1., 0.4);
        values::<T>(&vector(&psi), &v);
        let acted = psi.act(&op).unwrap();
        let y = apply(&expected, &v);
        values::<T>(&vector(&acted), &y);
        close::<T>(to64(psi.expectation(&op).unwrap()), dot(&v, &y));
        close::<T>(
            to64(psi.normalized_expectation(&op).unwrap()),
            dot(&v, &y) / dot(&v, &v),
        );
        close::<T>(
            psi.normalized()
                .unwrap()
                .norm()
                .unwrap()
                .to_f64()
                .unwrap()
                .into(),
            1.0.into(),
        );
        close::<T>(to64(psi.inner_product(&psi).unwrap()), dot(&v, &v));
        // Global vector adapters share exactly the same factor/state order.
        let array: Exported<qmb_export::ndarray::Array1<Complex<T>>, _> =
            acted.export(&SzBasis).unwrap();
        let sparse: Exported<qmb_export::sprs::CsVec<Complex<T>>, _> =
            acted.export(&SzBasis).unwrap();
        values::<T>(
            &array.data().iter().map(|z| to64(*z)).collect::<Vec<_>>(),
            &y,
        );
        values::<T>(
            &sparse
                .data()
                .to_dense()
                .iter()
                .map(|z| to64(*z))
                .collect::<Vec<_>>(),
            &y,
        );
    }
}

#[test]
fn tensor_bases_actions_and_all_backends_f32() {
    check::<f32>();
}
#[test]
fn tensor_bases_actions_and_all_backends_f64() {
    check::<f64>();
}
