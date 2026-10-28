//! Expand ordinary expressions, then obtain sparse coordinates for numeric export.

use num_complex::Complex;
use qmb_core::{HilbertSpace, QuantumSystem};
use qmb_spin_half::{
    CartesianBasis, LadderBasis, SpinHalf, SpinHalfOperator as Op, SpinHalfState as Ket, SzBasis,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    type Spin = SpinHalf<Complex<f64>>;
    let a = HilbertSpace::factor("a");
    let b = HilbertSpace::factor("b");
    let up_down = Spin::state(Ket::Up)
        .at(&a)?
        .tensor_product(&Spin::state(Ket::Down).at(&b)?)?;
    let down_up = Spin::state(Ket::Down)
        .at(&a)?
        .tensor_product(&Spin::state(Ket::Up).at(&b)?)?;
    let state = up_down
        .scaled(Complex::new(2.0, 0.0))
        .add(&up_down.scaled(Complex::new(3.0, 0.0)))?
        .add(&down_up.scaled(Complex::new(-1.0, 0.0)))?;

    // The result remains an ordinary State with the original space.
    let expanded = state.expand_in(&SzBasis)?;
    let coordinates = expanded.coordinates_in(&SzBasis)?;
    println!(
        "State coordinates (dimension {}):",
        coordinates.dimension()?
    );
    for entry in coordinates.flat_entries() {
        let (index, coefficient) = entry?;
        println!("  {index}: {coefficient}"); // 1: 5+0i, 2: -1+0i
    }

    // Re-expansion needs no record of the previous basis in the expression.
    let raising = Spin::operator(Op::SPlus).expand_in(&CartesianBasis)?;
    let raising = raising.expand_in(&LadderBasis)?;
    let coordinates = raising.coordinates_in(&LadderBasis)?;
    println!(
        "Raising operator has {} ladder-basis coordinate.",
        coordinates.len()
    );
    Ok(())
}
