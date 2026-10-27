use num_complex::Complex;
#[cfg(not(feature = "normalize"))]
use qmb::Divide;
use qmb::spin_half::utility::{heisenberg, neel};
use qmb::{Boundary, Chain, Id128, SpinHalfState, ToSpace};
#[cfg(feature = "normalize")]
use qmb::{Norm, NormalizedExpectation};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let chain = Chain::<Id128>::new(6, 1, Boundary::Open)?;
    let space = chain.to_space()?;
    // H = -J sum(Sx Sx + Sy Sy + Sz Sz); positive J is ferromagnetic.
    let h = heisenberg(&chain, &space, 1.0_f64)?;
    let psi = neel::<_, Complex<f64>>(&chain, &space, SpinHalfState::Up)?;
    #[cfg(feature = "normalize")]
    let energy = psi.normalized_expectation(&h)?;
    // This Neel state has nonzero norm; core's scalar quotient needs no
    // normalization extension or square-root implementation.
    #[cfg(not(feature = "normalize"))]
    let energy = psi.expectation(&h)?.divide(psi.inner_product(&psi)?);
    println!("Neel energy: {}", energy);
    #[cfg(feature = "normalize")]
    println!("Neel norm: {}", psi.norm()?);
    #[cfg(feature = "export")]
    {
        use qmb::export::{
            Export, Exported, MatrixArgs,
            ndarray::{Array1, Array2},
            sprs::{CsMat, CsVec},
        };
        use qmb::spin_half::SzBasis;
        // The Vec output is also usable with other numerical libraries.
        let raw: Exported<Vec<_>, _> = (&psi).export(&SzBasis)?;
        println!("Vector shape: {:?}", raw.metadata().shape());
        let vector: Exported<Array1<_>, _> = raw.export(())?;
        // Direct ndarray export internally follows the public Vec output path.
        let matrix: Exported<Array2<_>, _> =
            (&h).export(MatrixArgs::<SpinHalfState, _>::new(&SzBasis))?;
        let sparse_vector: Exported<CsVec<_>, _> = (&psi).export(&SzBasis)?;
        let sparse_matrix: Exported<CsMat<_>, _> =
            (&h).export(MatrixArgs::<SpinHalfState, _>::new(&SzBasis))?;
        let dense_action = matrix.data().dot(vector.data());
        let sparse_action = sparse_matrix.data() * sparse_vector.data();
        let dense_energy = vector.data().mapv(|z| z.conj()).dot(&dense_action);
        let sparse_energy = sparse_vector
            .data()
            .iter()
            .fold(Complex::new(0., 0.), |s, (i, &x)| {
                s + x.conj() * sparse_action.get(i).copied().unwrap_or_default()
            });
        assert!((dense_energy - energy).norm() < 1e-12);
        assert!((sparse_energy - energy).norm() < 1e-12);
        for (i, &value) in dense_action.iter().enumerate() {
            assert!((value - sparse_action.get(i).copied().unwrap_or_default()).norm() < 1e-12);
        }
        println!("Dense energy: {dense_energy}; sparse energy: {sparse_energy}");
        println!(
            "CSC shape: {:?}; nonzero entries: {}",
            sparse_matrix.data().shape(),
            sparse_matrix.data().nnz()
        );
    }
    Ok(())
}
