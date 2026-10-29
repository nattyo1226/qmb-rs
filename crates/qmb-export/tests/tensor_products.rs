#[path = "../../qmb-core/tests/support/graded.rs"]
mod graded;
use graded::*;
use qmb_core::{HilbertSpace, TensorKind};
use qmb_export::{Export, Exported, MatrixArgs, MatrixMetadata, VectorMetadata};
fn check<T: Real>() {
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        let factors = [
            HilbertSpace::factor(7),
            HilbertSpace::factor(1),
            HilbertSpace::factor(4),
        ];
        let h = HilbertSpace::tensor_with_kind(kind, factors.clone()).unwrap();
        let op = placed::<T>(A, &factors[2])
            .scaled(c(0.3, -0.7))
            .add(&placed(CREATE, &factors[0]))
            .unwrap();
        let expected: Vec<_> = embedding::<T>(A, 2, 3, kind)
            .into_iter()
            .zip(embedding(CREATE, 0, 3, kind))
            .map(|(a, b)| a * c::<T>(0.3, -0.7) + b)
            .collect();
        let dense: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> =
            (&op).export(MatrixArgs::new(&States).in_space(&h)).unwrap();
        values(dense.data(), &expected);
        assert_eq!(dense.metadata().space().effective_kind(), Some(kind));
        assert_eq!(
            dense
                .metadata()
                .space()
                .indices()
                .copied()
                .collect::<Vec<_>>(),
            [7, 1, 4]
        );
        let psi = ket(7, &factors, kind)
            .scaled(c(0.2, 0.8))
            .add(&ket(0, &factors, kind))
            .unwrap()
            .act(&op)
            .unwrap();
        let reference = vector(&psi);
        let dense_vector: Exported<Vec<_>, VectorMetadata<'_, _, _, Ket>> =
            (&psi).export(&States).unwrap();
        values(dense_vector.data(), &reference);
        assert_eq!(dense_vector.metadata().space().effective_kind(), Some(kind));
        #[cfg(feature = "ndarray")]
        {
            let matrix: Exported<qmb_export::ndarray::Array2<_>, MatrixMetadata<'_, Ket, _, _>> =
                (&op).export(MatrixArgs::new(&States).in_space(&h)).unwrap();
            values(matrix.data().as_slice().unwrap(), &expected);
            assert_eq!(matrix.metadata().space().effective_kind(), Some(kind));
            let v: Exported<qmb_export::ndarray::Array1<_>, VectorMetadata<'_, _, _, Ket>> =
                (&psi).export(&States).unwrap();
            values(v.data().as_slice().unwrap(), &reference);
        }
        #[cfg(feature = "sprs")]
        {
            let matrix: Exported<qmb_export::sprs::CsMat<_>, MatrixMetadata<'_, Ket, _, _>> =
                (&op).export(MatrixArgs::new(&States).in_space(&h)).unwrap();
            values(matrix.data().to_dense().as_slice().unwrap(), &expected);
            assert_eq!(matrix.metadata().space().effective_kind(), Some(kind));
            let v: Exported<qmb_export::sprs::CsVec<_>, VectorMetadata<'_, _, _, Ket>> =
                (&psi).export(&States).unwrap();
            values(v.data().to_dense().as_slice().unwrap(), &reference);
        }
        // Explicit space overrides the old ordinary binding, and invalid zero
        // references propagate before the matrix builder allocates storage.
        let old = op
            .bind(&HilbertSpace::tensor(factors.clone()).unwrap())
            .unwrap();
        let rebound: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> = (&old)
            .export(MatrixArgs::new(&States).in_space(&h))
            .unwrap();
        values(rebound.data(), &expected);
        let bad = placed::<T>(A, &HilbertSpace::factor(99)).scaled(c(0., 0.));
        let bad: Result<Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>>, _> =
            (&bad).export(MatrixArgs::new(&States).in_space(&h));
        assert!(matches!(
            bad,
            Err(qmb_export::ExportError::Expression(
                qmb_core::ExpressionError::MissingIndex
            ))
        ));
    }
}
#[test]
fn all_backends_f32() {
    check::<f32>();
}
#[test]
fn all_backends_f64() {
    check::<f64>();
}
