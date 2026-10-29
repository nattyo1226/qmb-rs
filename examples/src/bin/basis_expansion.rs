//! Expand ordinary expressions, then obtain sparse coordinates for numeric export.

use num_complex::Complex;
use qmb_core::HilbertSpace;
use qmb_spin_half::{
    CartesianBasis, LadderBasis, SpinHalf, SpinHalfOperator as Op, SpinHalfState as Ket, SzBasis,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    type Spin = SpinHalf<Complex<f64>>;
    let a = HilbertSpace::factor("a");
    let b = HilbertSpace::factor("b");
    let space = a.tensor(&b)?;
    let up_down = Spin::product_state(&space, [Ket::Up, Ket::Down])?;
    let down_up = Spin::product_state(&space, [Ket::Down, Ket::Up])?;
    let state = Spin::state_sum(
        &space,
        [
            up_down.scaled(2.0),
            up_down.scaled(3.0),
            down_up.scaled(-1.0),
        ],
    )?;

    // The result remains an ordinary State with the original space.
    let expanded = state.expand(&SzBasis)?;
    let coordinates = expanded.coordinates(&SzBasis)?;
    println!("State coordinates (dimension {}):", coordinates.dim()?);
    for entry in coordinates.flat_entries() {
        let (index, coefficient) = entry?;
        println!("  {index}: {coefficient}"); // 1: 5+0i, 2: -1+0i
    }

    // Re-expansion needs no record of the previous basis in the expression.
    let raising = Spin::operator(Op::SPlus).expand(&CartesianBasis)?;
    let raising = raising.expand(&LadderBasis)?;
    let coordinates = raising.coordinates(&LadderBasis)?;
    println!(
        "Raising operator has {} ladder-basis coordinate.",
        coordinates.len()
    );
    Ok(())
}
