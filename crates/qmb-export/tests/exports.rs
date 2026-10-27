use num_complex::Complex;
use qmb_core::{HilbertSpace, LocalOperator, LocalState, State};
use qmb_export::{Export, Exported, Layout, MatrixArgs, MatrixMetadata, VectorMetadata};
use qmb_spin_half::{
    CartesianBasis, LadderBasis, SpinHalfOperator as Op, SpinHalfState as Ket, SzBasis,
};
mod support;
use support::*;

type Ls<T> = LocalState<Ket, Op, Complex<T>>;
type Lo<T> = LocalOperator<Op, Complex<T>>;

fn check_local<T: Real>() {
    let basis = SzBasis;
    let psi = Ls::<T>::primitive(Ket::Up).scaled(c(0.3, -0.7))
        + Ls::primitive(Ket::Down).scaled(c(-0.2, 0.6));
    let reference = [c::<T>(0.3, -0.7), c(-0.2, 0.6)];
    let vector: Exported<Vec<_>, VectorMetadata<'_, _, (), Ket>> = (&psi).export(&basis).unwrap();
    values(vector.data(), &reference);
    assert_eq!(vector.metadata().shape(), [2]);
    assert_eq!(vector.metadata().factor_count(), 1);
    assert!(std::ptr::eq(vector.metadata().basis(), &basis));
    let coordinates = psi.coordinates_in(&basis).unwrap();
    let borrowed: Exported<Vec<_>, VectorMetadata<'_, _, (), ()>> =
        (&coordinates).export(()).unwrap();
    let owned: Exported<Vec<_>, VectorMetadata<'_, _, (), ()>> = coordinates.export(()).unwrap();
    values(borrowed.data(), owned.data());

    #[cfg(feature = "ndarray")]
    {
        let pointer = vector.data().as_ptr();
        let array: Exported<qmb_export::ndarray::Array1<_>, VectorMetadata<'_, _, (), Ket>> =
            vector.export(()).unwrap();
        let direct: Exported<qmb_export::ndarray::Array1<_>, VectorMetadata<'_, _, (), Ket>> =
            (&psi).export(&basis).unwrap();
        assert_eq!(pointer, array.data().as_ptr());
        values(array.data().as_slice().unwrap(), &reference);
        values(direct.data().as_slice().unwrap(), &reference);
        let coords = psi.coordinates_in(&basis).unwrap();
        let borrowed: Exported<qmb_export::ndarray::Array1<_>, VectorMetadata<'_, _, (), ()>> =
            (&coords).export(()).unwrap();
        let owned: Exported<qmb_export::ndarray::Array1<_>, VectorMetadata<'_, _, (), ()>> =
            coords.export(()).unwrap();
        values(borrowed.data().as_slice().unwrap(), &reference);
        values(owned.data().as_slice().unwrap(), &reference);
    }
    #[cfg(feature = "sprs")]
    {
        let sparse: Exported<qmb_export::sprs::CsVec<_>, VectorMetadata<'_, _, (), Ket>> =
            (&psi).export(&basis).unwrap();
        assert_eq!(sparse.data().indices(), &[0, 1]);
        values(sparse.data().data(), &reference);
    }

    let ops = [Op::I, Op::Sx, Op::Sy, Op::Sz, Op::SPlus, Op::SMinus];
    for a in ops {
        for b in ops {
            let op = Lo::<T>::primitive(a)
                .product(&Lo::primitive(b))
                .adjoint()
                .scaled(c(0.4, -0.2));
            let expected: Vec<_> = adjoint(&product(&local::<T>(a), &local(b), 2), 2)
                .into_iter()
                .map(|v| v * c(0.4, -0.2))
                .collect();
            let dense: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, ()>> =
                (&op).export(MatrixArgs::new(&basis)).unwrap();
            values(dense.data(), &expected);
            assert_eq!(dense.metadata().shape(), [2, 2]);
            assert_eq!(dense.metadata().layout(), Layout::RowMajor);
            let acted = psi.act(&op);
            let actual: Exported<Vec<_>, VectorMetadata<'_, _, (), Ket>> =
                (&acted).export(&basis).unwrap();
            values(actual.data(), &act(&expected, &reference));
            close(
                psi.expectation(&op),
                dot(&reference, &act(&expected, &reference)),
            );
            close(op.trace(), expected[0] + expected[3]);
            for expanded in [
                op.expand_in(&CartesianBasis).unwrap(),
                op.expand_in(&LadderBasis).unwrap(),
            ] {
                let data: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, ()>> =
                    (&expanded).export(MatrixArgs::new(&basis)).unwrap();
                values(data.data(), &expected);
            }
            #[cfg(feature = "ndarray")]
            {
                let ptr = dense.data().as_ptr();
                let array: Exported<
                    qmb_export::ndarray::Array2<_>,
                    MatrixMetadata<'_, Ket, _, ()>,
                > = dense.export(()).unwrap();
                let direct: Exported<
                    qmb_export::ndarray::Array2<_>,
                    MatrixMetadata<'_, Ket, _, ()>,
                > = (&op).export(MatrixArgs::new(&basis)).unwrap();
                assert_eq!(ptr, array.data().as_ptr());
                assert!(array.data().is_standard_layout());
                values(array.data().as_slice().unwrap(), &expected);
                values(direct.data().as_slice().unwrap(), &expected);
            }
            #[cfg(feature = "sprs")]
            {
                let sparse: Exported<qmb_export::sprs::CsMat<_>, MatrixMetadata<'_, Ket, _, ()>> =
                    (&op).export(MatrixArgs::new(&basis)).unwrap();
                assert!(sparse.data().is_csc());
                for row in 0..2 {
                    for col in 0..2 {
                        close(
                            sparse
                                .data()
                                .get(row, col)
                                .copied()
                                .unwrap_or_else(|| c(0., 0.)),
                            expected[row * 2 + col],
                        );
                    }
                }
            }
        }
    }
}

fn check_composite<T: Real>() {
    let f = HilbertSpace::factor(3u8);
    let g = HilbertSpace::factor(8u8);
    let space =
        HilbertSpace::tensor([HilbertSpace::tensor([f.clone()]).unwrap(), g.clone()]).unwrap();
    let basis = SzBasis;
    let inputs: Vec<_> = space
        .basis_states::<Ket, Op, Complex<T>, _>(&basis)
        .unwrap()
        .map(Result::unwrap)
        .collect();
    let x = [c(0.3, -0.7), c(-0.2, 0.6), c(0.4, 0.2), c(-0.9, -0.3)];
    let psi = State::sum(&space, inputs.iter().zip(x).map(|(s, c)| s.scaled(c))).unwrap();
    let a = Lo::<T>::primitive(Op::Sx)
        .at(&f)
        .unwrap()
        .scaled(c(0.2, 0.7))
        .add(&Lo::primitive(Op::Sy).at(&g).unwrap().scaled(c(-0.3, 0.4)))
        .unwrap();
    let b = Lo::<T>::primitive(Op::Sz)
        .at(&f)
        .unwrap()
        .add(
            &Lo::primitive(Op::SMinus)
                .at(&g)
                .unwrap()
                .scaled(c(0.1, -0.6)),
        )
        .unwrap();
    let am: Vec<_> = kron(&local::<T>(Op::Sx), &local(Op::I))
        .into_iter()
        .zip(kron(&local::<T>(Op::I), &local(Op::Sy)))
        .map(|(a, b)| a * c(0.2, 0.7) + b * c(-0.3, 0.4))
        .collect();
    let bm: Vec<_> = kron(&local::<T>(Op::Sz), &local(Op::I))
        .into_iter()
        .zip(kron(&local::<T>(Op::I), &local(Op::SMinus)))
        .map(|(a, b)| a + b * c(0.1, -0.6))
        .collect();
    let expected = adjoint(&product(&am, &bm, 4), 4);
    let op = a.product(&b).unwrap().adjoint(); // Unbound; explicit space is required.
    let dense: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> = (&op)
        .export(MatrixArgs::new(&basis).in_space(&space))
        .unwrap();
    values(dense.data(), &expected);
    assert!(dense.metadata().domain().shares_node_with(&space));
    assert!(dense.metadata().codomain().shares_node_with(&space));
    assert_eq!(
        dense
            .metadata()
            .space()
            .indices()
            .copied()
            .collect::<Vec<_>>(),
        vec![3, 8]
    );
    let vector: Exported<Vec<_>, VectorMetadata<'_, _, _, Ket>> = (&psi).export(&basis).unwrap();
    values(vector.data(), &x);
    assert!(vector.metadata().space().shares_node_with(&space));
    assert_eq!(vector.metadata().domain().indices().count(), 0);
    assert_eq!(vector.metadata().shape(), [4]);
    let acted = psi.act(&op).unwrap();
    let data: Exported<Vec<_>, VectorMetadata<'_, _, _, Ket>> = (&acted).export(&basis).unwrap();
    values(data.data(), &act(&expected, &x));
    close(psi.expectation(&op).unwrap(), dot(&x, &act(&expected, &x)));
    close(
        op.trace_in(&space).unwrap(),
        (0..4).fold(c(0., 0.), |s, i| s + expected[i * 4 + i]),
    );
    for expanded in [
        op.expand_in_space(&CartesianBasis, &space).unwrap(),
        op.expand_in_space(&LadderBasis, &space).unwrap(),
    ] {
        let d: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> =
            (&expanded).export(MatrixArgs::new(&basis)).unwrap();
        values(d.data(), &expected);
    }
    #[cfg(feature = "ndarray")]
    {
        let ptr = dense.data().as_ptr();
        let array: Exported<qmb_export::ndarray::Array2<_>, MatrixMetadata<'_, Ket, _, _>> =
            dense.export(()).unwrap();
        let direct: Exported<qmb_export::ndarray::Array2<_>, MatrixMetadata<'_, Ket, _, _>> = (&op)
            .export(MatrixArgs::new(&basis).in_space(&space))
            .unwrap();
        assert_eq!(ptr, array.data().as_ptr());
        values(
            array.data().as_slice().unwrap(),
            direct.data().as_slice().unwrap(),
        );
        values(array.data().as_slice().unwrap(), &expected);
        let ptr = vector.data().as_ptr();
        let array: Exported<qmb_export::ndarray::Array1<_>, VectorMetadata<'_, _, _, Ket>> =
            vector.export(()).unwrap();
        let direct: Exported<qmb_export::ndarray::Array1<_>, VectorMetadata<'_, _, _, Ket>> =
            (&psi).export(&basis).unwrap();
        assert_eq!(ptr, array.data().as_ptr());
        values(
            array.data().as_slice().unwrap(),
            direct.data().as_slice().unwrap(),
        );
    }
    #[cfg(feature = "sprs")]
    {
        let sparse: Exported<qmb_export::sprs::CsMat<_>, MatrixMetadata<'_, Ket, _, _>> = (&op)
            .export(MatrixArgs::new(&basis).in_space(&space))
            .unwrap();
        let v: Exported<qmb_export::sprs::CsVec<_>, VectorMetadata<'_, _, _, Ket>> =
            (&psi).export(&basis).unwrap();
        for row in 0..4 {
            for col in 0..4 {
                close(
                    sparse
                        .data()
                        .get(row, col)
                        .copied()
                        .unwrap_or_else(|| c(0., 0.)),
                    expected[row * 4 + col],
                );
            }
        }
        let result = sparse.data() * v.data();
        for (i, &expected) in act(&expected, &x).iter().enumerate() {
            close(
                result.get(i).copied().unwrap_or_else(|| c(0., 0.)),
                expected,
            );
        }
        let coordinates = psi.coordinates_in(&basis).unwrap();
        let borrowed: Exported<qmb_export::sprs::CsVec<_>, VectorMetadata<'_, _, _, ()>> =
            (&coordinates).export(()).unwrap();
        let moved: Exported<qmb_export::sprs::CsVec<_>, VectorMetadata<'_, _, _, ()>> =
            coordinates.export(()).unwrap();
        assert_eq!(borrowed.data(), moved.data());
    }
    // Grouping and a re-created compatible space preserve the flat order.
    let flat = f.tensor_product(&g).unwrap();
    let other: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> = (&op)
        .export(MatrixArgs::new(&basis).in_space(&flat))
        .unwrap();
    values(other.data(), &expected);
    // A source can be dropped while the result retains the basis and shared space.
    drop(psi);
    assert_eq!(
        data.metadata()
            .space()
            .indices()
            .copied()
            .collect::<Vec<_>>(),
        vec![3, 8]
    );
}

#[test]
fn local_exports_match_independent_spin_matrices_both_precisions() {
    check_local::<f32>();
    check_local::<f64>();
}
#[test]
fn composite_exports_match_independent_four_by_four_matrices_both_precisions() {
    check_composite::<f32>();
    check_composite::<f64>();
}

#[test]
fn vector_data_can_be_explicitly_converted_between_precisions() {
    let psi = Ls::<f32>::primitive(Ket::Down).scaled(c(0.3, -0.7));
    let output: Exported<Vec<Complex<f32>>, _> = (&psi).export(&SzBasis).unwrap();
    let (data, metadata) = output.into_parts();
    let converted: Vec<Complex<f64>> = data
        .into_iter()
        .map(|c| Complex::new(f64::from(c.re), f64::from(c.im)))
        .collect();
    assert_eq!(metadata.shape(), [2]);
    assert_eq!(
        converted[1],
        Complex::new(f64::from(0.3_f32), f64::from(-0.7_f32))
    );
}
