//! Select a tensor rule at construction and use the common evaluation API.
use num_complex::Complex;
use qmb_core::{
    HilbertSpace, HomogeneousParity, Parity, ParityOperator, QuantumSystem, TensorKind,
};
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
        let space = HilbertSpace::tensor_with_kind(kind, [a.clone(), b.clone()])?;
        let psi = Spin::state_tensor_at_with_kind(kind, [(&a, Ket::Up), (&b, Ket::Down)])?;
        let op = Spin::operator(Op::SPlus).at(&b)?.bind(&space)?;
        assert!((psi.act(&op)?.norm_squared()? - 1.0).abs() < 1e-12);
        assert!((op.adjoint().product(&op)?.trace()?.re - 2.0).abs() < 1e-12);
        println!(
            "{kind:?}: action norm squared = {}",
            psi.act(&op)?.norm_squared()?
        );
    }
    Ok(())
}
