use core::{
    cell::Cell,
    ops::{Add, Mul},
};
use num_complex::Complex;
use num_traits::Zero;
use qmb_core::{
    Act, Adjoint, Conjugate, HilbertSpace, InnerProduct, LinearCombination, LocalOperator,
    LocalState, RealPart,
};
use qmb_normalize::{Norm, Normalize, NormalizedExpectation};

struct Primitive<'a>(&'a Cell<usize>);

// A coefficient with a floating real projection, but without One, From<f64>,
// division, Copy, PartialEq or Float. Norm alone must not need normalization's
// real-embedding capability.
#[derive(Clone)]
struct NormOnly(f64);
impl Add for NormOnly {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}
impl Mul for NormOnly {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0)
    }
}
impl Zero for NormOnly {
    fn zero() -> Self {
        Self(0.0)
    }
    fn is_zero(&self) -> bool {
        self.0 == 0.0
    }
}
impl Conjugate for NormOnly {
    type Output = Self;
    fn conjugate(&self) -> Self {
        self.clone()
    }
}
impl RealPart for NormOnly {
    type Real = f64;
    fn real_part(self) -> f64 {
        self.0
    }
}
impl InnerProduct<NormOnly> for () {
    fn inner_product(&self, _: &Self) -> NormOnly {
        NormOnly(1.0)
    }
}
impl Act<NormOnly> for () {
    type Operator = Self;
    type Term = Self;
    fn act(&self, _: &Self) -> LinearCombination<Self, NormOnly> {
        LinearCombination::term(NormOnly(1.0), ())
    }
}
impl Adjoint<NormOnly> for () {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, NormOnly> {
        LinearCombination::term(NormOnly(1.0), ())
    }
}

#[test]
fn norm_does_not_require_inverse_real_embedding_or_coefficient_float() {
    let state = LocalState::<(), (), NormOnly>::primitive(()).scaled(NormOnly(3.0));
    let norm: f64 = state.norm();
    assert_eq!(norm, 3.0);
    assert_eq!(LocalState::<(), (), NormOnly>::zero().norm(), 0.0);
}

struct Op;

impl qmb_core::Identity for Op {
    fn identity() -> Self {
        Self
    }
}

impl InnerProduct<Complex<f64>> for Primitive<'_> {
    fn inner_product(&self, _: &Self) -> Complex<f64> {
        self.0.set(self.0.get() + 1);
        Complex::new(1.0, 0.0)
    }
}
impl<'a> Act<Complex<f64>> for Primitive<'a> {
    type Operator = Op;
    type Term = Self;
    fn act(&self, _: &Op) -> LinearCombination<Self, Complex<f64>> {
        LinearCombination::primitive(Self(self.0))
    }
}
impl Adjoint<Complex<f64>> for Op {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, Complex<f64>> {
        LinearCombination::primitive(Op)
    }
}

#[test]
fn normalization_preserves_phase_and_contracts_once() {
    let calls = Cell::new(0);
    let basis = LocalState::<Primitive<'_>, Op, Complex<f64>>::primitive(Primitive(&calls));
    let state = basis.scaled(Complex::new(3.0, 4.0));
    let normalized = state.normalized().unwrap();
    assert_eq!(calls.get(), 1);
    let amplitude = basis.inner_product(&normalized);
    assert!((amplitude.re - 0.6).abs() < 1e-12);
    assert!((amplitude.im - 0.8).abs() < 1e-12);
    let space = HilbertSpace::factor(0u8);
    let state = state.at(&space).unwrap();
    calls.set(0);
    let normalized = state.normalized().unwrap();
    assert_eq!(calls.get(), 1);
    assert!((normalized.norm().unwrap() - 1.0).abs() < 1e-12);
    // Normalization has not evaluated, replaced or mutated the original input.
    assert_eq!(state.norm().unwrap(), 5.0);
    let identity = LocalOperator::identity().at(&space).unwrap();
    assert_eq!(
        normalized.normalized_expectation(&identity).unwrap(),
        Complex::new(1.0, 0.0)
    );
}
