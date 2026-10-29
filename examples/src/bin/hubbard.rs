//! Standard and symmetric Hubbard models on a two-site graded chain.
use num_complex::Complex;
use qmb_basic::Id128;
use qmb_core::HilbertSpace;
use qmb_fermion::{Fermion, FermionOperator as O, FermionState as S, MajoranaBasis};
use qmb_geometry::{Boundary, Chain};
#[cfg(feature = "normalize")]
use qmb_normalize::{Normalize, NormalizedExpectation};

type F = Fermion<Complex<f64>>;
type Models = qmb_fermion::models::FermionModels<F>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let chain = Chain::new(2, Boundary::Open);
    // Site-major ordering: up_0, down_0, up_1, down_1.
    let modes: Vec<_> = (0..4).map(|_| HilbertSpace::factor(Id128::new())).collect();
    let sites: Vec<_> = modes
        .chunks_exact(2)
        .map(|pair| HilbertSpace::tensor_all_g(pair.iter().cloned()))
        .collect::<Result<_, _>>()?;
    let h = Models::hubbard(&chain, &sites, 1., 4.)?;
    let symmetric = Models::symmetric_hubbard(&chain, &sites, 1., 4.)?;
    let psi = F::product_state(
        h.space().expect("model is bound"),
        [S::Occupied, S::Empty, S::Empty, S::Occupied],
    )?;
    assert!((psi.norm_squared()? - 1.).abs() < 1e-12);
    let energy = psi.expectation(&h)?;
    let shifted_energy = psi.expectation(&symmetric)?;
    assert!(energy.norm_sqr() < 1e-24);
    assert!((shifted_energy - Complex::new(-2., 0.)).norm_sqr() < 1e-24);
    println!("|up, down>: Hubbard energy {energy}, symmetric energy {shifted_energy}");
    #[cfg(feature = "normalize")]
    {
        let normalized = psi.scaled(Complex::new(0., 2.)).normalized()?;
        assert!(
            (normalized.normalized_expectation(&symmetric)? - shifted_energy).norm_sqr() < 1e-24
        );
    }
    let p = F::operator(O::Parity).coordinates(&MajoranaBasis)?;
    for entry in p.flat_entries() {
        let (index, coefficient) = entry?;
        println!("Parity in MajoranaBasis: coordinate {index} = {coefficient}");
    }
    #[cfg(feature = "export")]
    {
        use qmb_export::{Export, Exported, MatrixArgs};
        use qmb_fermion::OccupationBasis;
        let dense: Exported<Vec<Complex<f64>>, _> =
            (&h).export(MatrixArgs::<S, _>::new(&OccupationBasis))?;
        let sparse: Exported<qmb_export::sprs::CsMat<Complex<f64>>, _> =
            (&h).export(MatrixArgs::<S, _>::new(&OccupationBasis))?;
        assert_eq!(dense.metadata().shape(), [16, 16]);
        for (left, right) in dense.data().iter().zip(sparse.data().to_dense().iter()) {
            assert!((*left - *right).norm_sqr() < 1e-24);
        }
        println!(
            "Hubbard matrix: {:?}, {} nonzero entries",
            dense.metadata().shape(),
            sparse.data().nnz()
        );
    }
    Ok(())
}
