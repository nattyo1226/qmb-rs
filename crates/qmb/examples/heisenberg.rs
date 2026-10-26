use num_complex::Complex;
#[cfg(not(feature = "normalize"))]
use qmb::Divide;
use qmb::spin_half::utility::{heisenberg, neel};
use qmb::{Boundary, Chain, Id128, SpinHalfState};
#[cfg(feature = "normalize")]
use qmb::{Norm, NormalizedExpectation};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let chain = Chain::<Id128>::new(6, 1, Boundary::Open)?;
    // H = -J sum(Sx Sx + Sy Sy + Sz Sz); positive J is ferromagnetic.
    let h = heisenberg(&chain, 1.0_f64)?;
    let psi = neel::<_, Complex<f64>>(&chain, SpinHalfState::Up)?;
    #[cfg(feature = "normalize")]
    let energy = psi.normalized_expectation(&h)?;
    // This Neel state has nonzero norm; core's scalar quotient needs no
    // normalization extension or square-root implementation.
    #[cfg(not(feature = "normalize"))]
    let energy = psi.expectation(&h)?.divide(psi.inner_product(&psi)?);
    println!("Neel energy: {}", energy);
    #[cfg(feature = "normalize")]
    println!("Neel norm: {}", psi.norm()?);
    Ok(())
}
