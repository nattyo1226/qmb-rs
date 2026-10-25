use qmb::spin_half::utility::{heisenberg, neel};
use qmb::{Boundary, Chain, Id128, SpinHalf, SpinHalfState};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let chain = Chain::<Id128, SpinHalf>::new(6, 1, Boundary::Open)?;
    // H = -J sum(Sx Sx + Sy Sy + Sz Sz); positive J is ferromagnetic.
    let h = heisenberg(&chain, 1.0_f64)?;
    let psi = neel(&chain, SpinHalfState::Up)?;
    println!("Neel energy: {}", psi.normalized_expectation(&h)?);
    #[cfg(feature = "libm")]
    println!("Neel norm: {}", psi.norm()?);
    Ok(())
}
