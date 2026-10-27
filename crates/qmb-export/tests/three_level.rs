use num_complex::Complex;
use qmb_core::{Basis, HilbertSpace, LocalOperator, LocalState};
use qmb_export::{Export, Exported, MatrixArgs, MatrixMetadata, VectorMetadata};
#[path = "../../qmb-core/tests/support/three_level.rs"]
mod fixture;
use fixture::*;

// Follow the required complex-magnitude atol/rtol comparison, independently of
// the fixture's component-wise assertions used by older core tests.
fn close<T: Real>(actual: Complex<T>, expected: Complex<T>) {
    let abs = |z: Complex<T>| z.re.to_f64().unwrap().hypot(z.im.to_f64().unwrap());
    assert!(
        abs(actual - expected) <= T::TOLERANCE + T::TOLERANCE * abs(actual).max(abs(expected)),
        "{actual:?} != {expected:?}"
    );
}

fn check<T: Real>() {
    type S<T> = LocalState<Ket, Operator, Complex<T>>;
    type O<T> = LocalOperator<Operator, Complex<T>>;
    // The change-of-basis matrix has columns |a>, |a>+|b>, |c>.
    let b = [
        [c::<T>(1., 0.), c(1., 0.), c(0., 0.)],
        [c(0., 0.), c(1., 0.), c(0., 0.)],
        [c(0., 0.), c(0., 0.), c(1., 0.)],
    ];
    let inverse = [
        [c::<T>(1., 0.), c(-1., 0.), c(0., 0.)],
        [c(0., 0.), c(1., 0.), c(0., 0.)],
        [c(0., 0.), c(0., 0.), c(1., 0.)],
    ];
    let gram = matmul(
        std::array::from_fn(|i| std::array::from_fn(|j| b[j][i].conj())),
        b,
    );
    let basis = NonOrthogonal;
    let psi = S::<T>::primitive(Ket::Phase).scaled(c(0.3, -0.7))
        + S::primitive(Ket::Sum).scaled(c(-0.2, 0.6));
    let x = std::array::from_fn(|i| {
        Ket::Phase.coordinates::<T>()[i] * c(0.3, -0.7)
            + Ket::Sum.coordinates::<T>()[i] * c(-0.2, 0.6)
    });
    let expected_coords = matvec(inverse, x);
    let coords: Exported<Vec<_>, VectorMetadata<'_, _, (), Ket>> = (&psi).export(&basis).unwrap();
    for (a, b) in coords.data().iter().zip(expected_coords) {
        close(*a, b);
    }
    close(
        dot(expected_coords, matvec(gram, expected_coords)),
        psi.inner_product(&psi),
    );
    for a in [
        Operator::Identity,
        Operator::Unit(Axis::A, Axis::B),
        Operator::Unit(Axis::C, Axis::A),
    ] {
        let op = O::<T>::primitive(a).scaled(c(0.7, -0.2));
        let physical =
            std::array::from_fn(|i| std::array::from_fn(|j| matrix::<T>(a)[i][j] * c(0.7, -0.2)));
        let expected = matmul(matmul(inverse, physical), b);
        let output: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, ()>> =
            (&op).export(MatrixArgs::new(&basis)).unwrap();
        for (i, row) in expected.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                close(output.data()[i * 3 + j], value);
            }
        }
        let y = matvec(expected, expected_coords);
        close(dot(expected_coords, matvec(gram, y)), psi.expectation(&op));
        close(expected[0][0] + expected[1][1] + expected[2][2], op.trace());
        let expanded = op.expand_in(&MatrixUnits).unwrap();
        let again: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, ()>> =
            (&expanded).export(MatrixArgs::new(&basis)).unwrap();
        for (&a, &b) in again.data().iter().zip(output.data()) {
            close(a, b);
        }
        #[cfg(feature = "ndarray")]
        {
            let direct: Exported<qmb_export::ndarray::Array2<_>, MatrixMetadata<'_, Ket, _, ()>> =
                (&op).export(MatrixArgs::new(&basis)).unwrap();
            for (i, row) in expected.iter().enumerate() {
                for (j, &value) in row.iter().enumerate() {
                    close(direct.data()[[i, j]], value);
                }
            }
            let converted: Exported<
                qmb_export::ndarray::Array2<_>,
                MatrixMetadata<'_, Ket, _, ()>,
            > = output.export(()).unwrap();
            assert_eq!(direct.data(), converted.data());
        }
        #[cfg(feature = "sprs")]
        {
            let sparse: Exported<qmb_export::sprs::CsMat<_>, MatrixMetadata<'_, Ket, _, ()>> =
                (&op).export(MatrixArgs::new(&basis)).unwrap();
            for (i, row) in expected.iter().enumerate() {
                for (j, &value) in row.iter().enumerate() {
                    close(
                        sparse
                            .data()
                            .get(i, j)
                            .copied()
                            .unwrap_or_else(|| c(0., 0.)),
                        value,
                    );
                }
            }
        }
    }
    // Non-orthogonal tensor-product bases, with a spectator identity.
    let f = HilbertSpace::factor(0u8);
    let g = HilbertSpace::factor(1u8);
    let space = f.tensor_product(&g).unwrap();
    let op = O::<T>::primitive(Operator::Unit(Axis::A, Axis::B))
        .at(&f)
        .unwrap();
    let local = matmul(matmul(inverse, matrix(Operator::Unit(Axis::A, Axis::B))), b);
    let output: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> = (&op)
        .export(MatrixArgs::new(&basis).in_space(&space))
        .unwrap();
    for row in 0..9 {
        for col in 0..9 {
            close(
                output.data()[row * 9 + col],
                if row % 3 == col % 3 {
                    local[row / 3][col / 3]
                } else {
                    c(0., 0.)
                },
            );
        }
    }
    let inputs = space
        .basis_states::<Ket, Operator, Complex<T>, _>(&basis)
        .unwrap();
    for (col, input) in inputs.enumerate() {
        let acted = input.unwrap().act(&op).unwrap();
        let coords: Exported<Vec<_>, VectorMetadata<'_, _, _, Ket>> =
            (&acted).export(&basis).unwrap();
        for row in 0..9 {
            close(coords.data()[row], output.data()[row * 9 + col]);
        }
    }
    // Identity is not itself a MatrixUnits basis element; export remains valid.
    let identity = qmb_core::Operator::<u8, Operator, Complex<T>>::identity()
        .bind(&space)
        .unwrap();
    let output: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> =
        (&identity).export(MatrixArgs::new(&basis)).unwrap();
    for row in 0..9 {
        for col in 0..9 {
            close(
                output.data()[row * 9 + col],
                c(if row == col { 1. } else { 0. }, 0.),
            );
        }
    }
    let canonical: Vec<_> = <Canonical as Basis<Ket, Complex<T>>>::elements(&Canonical).collect();
    assert_eq!(canonical.len(), 3);
}

#[test]
fn three_level_coordinates_actions_and_gram_expectations_both_precisions() {
    check::<f32>();
    check::<f64>();
}
