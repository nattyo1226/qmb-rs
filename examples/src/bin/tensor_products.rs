//! Select a tensor rule at construction and use the common evaluation API.
use num_complex::Complex;
use qmb_core::{HilbertSpace, HomogeneousParity, Parity, ParityOperator, TensorKind};
use qmb_spin_half::{SpinHalf, SpinHalfOperator as Op, SpinHalfState as Ket};
type Spin = SpinHalf<Complex<f64>>;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = HilbertSpace::factor("a");
    let b = HilbertSpace::factor("b");
    // Spin-half has trivial grading: every primitive is even and P = I.
    assert_eq!(Op::Sx.parity(), Parity::Even);
    assert_eq!(Ket::Down.parity(), Parity::Even);
    assert_eq!(Op::parity_operator(), Op::I);
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        let space = HilbertSpace::tensor_all_with(kind, [a.clone(), b.clone()])?;
        let psi = Spin::product_state(&space, [Ket::Up, Ket::Down])?;
        let op = Spin::operator(Op::SPlus).at(&b)?.bind(&space)?;
        assert!((psi.act(&op)?.norm_squared()? - 1.0).abs() < 1e-12);
        assert!(((op.adjoint().clone() * &op)?.trace()?.re - 2.0).abs() < 1e-12);
        println!(
            "{kind:?}: action norm squared = {}",
            psi.act(&op)?.norm_squared()?
        );
    }
    Ok(())
}
