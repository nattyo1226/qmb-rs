//! Spin-1 dipoles, quadrupoles, chain models and common expression APIs.

use num_complex::Complex;
use qmb_core::{HilbertSpace, QuantumSystem};
use qmb_geometry::{Boundary, Chain};
use qmb_spin_one::{
    CartesianBasis, LadderBasis, SpinOne, SpinOneOperator as Op, SpinOneState as Ket,
};

fn assert_close(actual: Complex<f64>, expected: Complex<f64>) {
    let delta = actual - expected;
    assert!(delta.re.hypot(delta.im) < 1e-12);
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    type Spin = SpinOne<Complex<f64>>;
    let psi = Spin::state(Ket::Plus)
        + Spin::state(Ket::Zero).scaled(Complex::new(0.0, 1.0))
        + Spin::state(Ket::Minus);
    let z = Spin::operator(Op::Sz);
    // Axial and rhombic single-ion anisotropy, plus a transverse field.
    let h = z.product(&z).scaled(0.4)
        + Spin::operator(Op::Dxy).scaled(0.2)
        + Spin::operator(Op::Sx).scaled(0.3);
    let energy = psi.expectation(&h);
    assert_close(energy, Complex::new(1.2, 0.0));
    println!("Unnormalized energy: {energy}");

    let plus = Spin::operator(Op::SPlus);
    let squared = plus.product(&plus);
    let coordinates = squared.coordinates_in(&LadderBasis)?;
    assert_eq!(coordinates.dimension()?, 9);
    assert_eq!(coordinates.len(), 1);
    for entry in coordinates.flat_entries() {
        let (index, coefficient) = entry?;
        assert_eq!(index, 7); // Q2Plus, with coefficient one.
        assert_eq!(*coefficient, Complex::new(1.0, 0.0));
        println!("SPlus²: ladder coordinate {index} = {coefficient}");
    }
    let round = squared
        .expand_in(&CartesianBasis)?
        .expand_in(&LadderBasis)?;
    assert_eq!(round.coordinates_in(&LadderBasis)?.len(), 1);

    // The same primitives work on independently identified tensor factors.
    let a = HilbertSpace::factor("a");
    let b = HilbertSpace::factor("b");
    let space = a.tensor_product(&b)?;
    let state = Spin::product_state(&space, [Ket::Plus, Ket::Minus])?;
    let zz = Spin::operator_product_at([(&a, Op::Sz), (&b, Op::Sz)])?;
    assert_eq!(state.space().dimension::<Ket>()?, 9);
    assert_eq!(state.expectation(&zz)?, Complex::new(-1.0, 0.0));

    // The model helpers use -J for bilinear coupling, as in spin-half.
    // A two-site periodic chain has one bond, counted once.
    let sites = [a, b];
    let chain = Chain::new(2, Boundary::Periodic);
    let xyz = Spin::xyz(&chain, &sites, 1.0, 1.0, 1.0)?;
    let heisenberg = Spin::heisenberg(&chain, &sites, 1.0)?;
    let blbq = Spin::blbq(&chain, &sites, 1.0, 0.25)?;
    assert_close(state.expectation(&xyz)?, Complex::new(1.0, 0.0));
    assert_close(state.expectation(&heisenberg)?, Complex::new(1.0, 0.0));
    let model_energy = state.expectation(&blbq)?;
    assert_close(model_energy, Complex::new(0.5, 0.0));
    println!("BLBQ energy in |+,->: {model_energy}");

    #[cfg(feature = "normalize")]
    {
        use qmb_normalize::{Norm, Normalize, NormalizedExpectation};
        let normalized = psi.normalized()?;
        assert!((normalized.norm() - 1.0).abs() < 1e-12);
        let energy = psi.normalized_expectation(&h)?;
        assert_close(energy, Complex::new(0.4, 0.0));
        println!("Normalized energy: {energy}");
    }
    #[cfg(feature = "export")]
    {
        use qmb_export::{Export, Exported, MatrixArgs};
        use qmb_spin_one::SzBasis;
        let dense: Exported<Vec<Complex<f64>>, _> =
            (&h).export(MatrixArgs::<Ket, _>::new(&SzBasis))?;
        let sparse: Exported<qmb_export::sprs::CsMat<Complex<f64>>, _> =
            (&h).export(MatrixArgs::<Ket, _>::new(&SzBasis))?;
        assert_eq!(dense.metadata().shape(), [3, 3]);
        for (left, right) in dense.data().iter().zip(sparse.data().to_dense().iter()) {
            assert_close(*left, *right);
        }
        println!(
            "Hamiltonian: 3×3 matrix, {} sparse entries",
            sparse.data().nnz()
        );

        let dense: Exported<Vec<Complex<f64>>, _> =
            blbq.export(MatrixArgs::<Ket, _>::new(&SzBasis))?;
        let sparse: Exported<qmb_export::sprs::CsMat<Complex<f64>>, _> =
            blbq.export(MatrixArgs::<Ket, _>::new(&SzBasis))?;
        assert_eq!(dense.metadata().shape(), [9, 9]);
        for (left, right) in dense.data().iter().zip(sparse.data().to_dense().iter()) {
            assert_close(*left, *right);
        }
        println!("BLBQ: 9×9 matrix, {} sparse entries", sparse.data().nnz());
    }
    Ok(())
}
