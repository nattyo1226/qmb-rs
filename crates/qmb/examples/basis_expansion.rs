//! Expand ordinary expressions, then obtain sparse coordinates for numeric export.

use num_complex::Complex;
use qmb::spin_half::{CartesianBasis, LadderBasis, SzBasis};
use qmb::{HilbertSpace, LocalOperator, LocalState, SpinHalfOperator as Op, SpinHalfState as Ket};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    type S = LocalState<Ket, Op, Complex<f64>>;
    type O = LocalOperator<Op, Complex<f64>>;
    let a = HilbertSpace::factor("a");
    let b = HilbertSpace::factor("b");
    let up_down = S::primitive(Ket::Up)
        .at(&a)?
        .tensor_product(&S::primitive(Ket::Down).at(&b)?)?;
    let down_up = S::primitive(Ket::Down)
        .at(&a)?
        .tensor_product(&S::primitive(Ket::Up).at(&b)?)?;
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
    let raising = O::primitive(Op::SPlus).expand_in(&CartesianBasis)?;
    let raising = raising.expand_in(&LadderBasis)?;
    let coordinates = raising.coordinates_in(&LadderBasis)?;
    println!(
        "Raising operator has {} ladder-basis coordinate.",
        coordinates.len()
    );
    Ok(())
}
