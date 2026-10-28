//! Build singlet/triplet states and an exchange Hamiltonian using QuantumSystem.

use num_complex::Complex;
use qmb_core::{HilbertSpace, OperatorOf, QuantumSystem};
use qmb_spin_half::{SpinHalf, SpinHalfOperator as Op, SpinHalfState as Ket};

// Choose primitives and coefficient precision once for states and operators.
type Spin = SpinHalf<Complex<f64>>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = HilbertSpace::factor("a");
    let b = HilbertSpace::factor("b");
    let space = a.tensor_product(&b)?;

    let up = Spin::state(Ket::Up);
    let down = Spin::state(Ket::Down);
    let up_down = up.at(&a)?.tensor_product(&down.at(&b)?)?;
    let down_up = down.at(&a)?.tensor_product(&up.at(&b)?)?;
    let amplitude = Complex::new(1.0 / 2.0_f64.sqrt(), 0.0);
    // |singlet> = (|up,down> - |down,up>) / sqrt(2)
    // |triplet> = (|up,down> + |down,up>) / sqrt(2)
    let singlet = up_down
        .add(&down_up.scaled(Complex::new(-1.0, 0.0)))?
        .scaled(amplitude);
    let triplet = up_down.add(&down_up)?.scaled(amplitude);

    // H = S_a . S_b, with positive exchange coupling J = 1.
    // Build the interaction explicitly from local operators.
    let mut terms = Vec::new();
    for axis in [Op::Sx, Op::Sy, Op::Sz] {
        let local = Spin::operator(axis);
        terms.push(local.at(&a)?.product(&local.at(&b)?)?);
    }
    let h = OperatorOf::<Spin, _>::sum(terms)?.bind(&space)?;

    for (name, state, expected_energy) in
        [("singlet", &singlet, -0.75), ("triplet", &triplet, 0.25)]
    {
        let norm_squared = state.norm_squared()?;
        let energy = state.expectation(&h)?;
        // act constructs H|state>; the residual checks the eigenstate relation.
        let residual = state
            .act(&h)?
            .add(&state.scaled(Complex::new(-expected_energy, 0.0)))?;
        let residual_norm_squared = residual.norm_squared()?;
        println!(
            "{name}: norm squared = {norm_squared:.6}, <H> = {:.6}, \
             ||H psi - E psi||^2 = {residual_norm_squared:.6}",
            energy.re,
        );
        assert!((norm_squared - 1.0).abs() < 1e-12);
        assert!((energy - Complex::new(expected_energy, 0.0)).norm_sqr() < 1e-24);
        assert!(residual_norm_squared.abs() < 1e-12);
    }
    let overlap = singlet.inner_product(&triplet)?;
    println!("<singlet|triplet> = {overlap}");
    assert!(overlap.norm_sqr() < 1e-24);

    Ok(())
}
