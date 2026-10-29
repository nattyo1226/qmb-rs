//! Build singlet/triplet states and an exchange Hamiltonian using QuantumSystem.

use num_complex::Complex;
use qmb_core::{HilbertSpace, QuantumSystem};
use qmb_spin_half::{SpinHalf, SpinHalfOperator as Op, SpinHalfState as Ket};

// Choose primitives and coefficient precision once for states and operators.
type Spin = SpinHalf<Complex<f64>>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = HilbertSpace::factor("a");
    let b = HilbertSpace::factor("b");
    let space = a.tensor_product(&b)?;

    let up_down = Spin::product_state(&space, [Ket::Up, Ket::Down])?;
    let down_up = Spin::product_state(&space, [Ket::Down, Ket::Up])?;
    let amplitude = 1.0 / 2.0_f64.sqrt();
    // |singlet> = (|up,down> - |down,up>) / sqrt(2)
    // |triplet> = (|up,down> + |down,up>) / sqrt(2)
    let singlet = (&up_down - &down_up)?.scaled(amplitude);
    let triplet = (&up_down + &down_up)?.scaled(amplitude);

    // H = S_a . S_b, with positive exchange coupling J = 1.
    // Each axis builds a placed product; try_operator_sum collects the terms.
    let h = Spin::try_operator_sum(
        [Op::Sx, Op::Sy, Op::Sz]
            .into_iter()
            .map(|axis| Spin::operator_product_at([(&a, axis), (&b, axis)])),
    )?
    .bind(&space)?;

    for (name, state, expected_energy) in
        [("singlet", &singlet, -0.75), ("triplet", &triplet, 0.25)]
    {
        let norm_squared = state.norm_squared()?;
        let energy = state.expectation(&h)?;
        // act constructs H|state>; the residual checks the eigenstate relation.
        let residual = (state.act(&h)? - state.scaled(expected_energy))?;
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
