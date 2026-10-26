//! Coefficient capabilities are independent of storage and of each other.
use core::ops::{Add, Mul};
use num_complex::Complex;
use num_traits::{One, Zero};
use qmb_core::{
    Act, Adjoint, Conjugate, Divide, HilbertSpace, InnerProduct, LinearCombination, LocalOperator,
    LocalState,
};

struct Primitive;

impl<C: One> Act<C> for Primitive {
    type Operator = Self;
    type Term = Self;
    fn act(&self, _: &Self) -> LinearCombination<Self, C> {
        LinearCombination::primitive(Self)
    }
}
impl<C: One> Adjoint<C> for Primitive {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, C> {
        LinearCombination::primitive(Self)
    }
}
impl<C: One> InnerProduct<C> for Primitive {
    fn inner_product(&self, _: &Self) -> C {
        C::one()
    }
}

#[test]
fn opaque_whole_coefficients_move_through_map_and_shared_local_and_composite_syntax() {
    struct Opaque(u8);
    let terms = LinearCombination::term(Opaque(7), Primitive).map(|_| 9);
    let term = terms.into_iter().next().unwrap();
    assert_eq!((term.coefficient.0, term.operand), (7, 9));
    let state = LocalState::<Primitive, Primitive, Opaque>::primitive(Primitive).scaled(Opaque(2));
    let operator = LocalOperator::primitive(Primitive).scaled(Opaque(3));
    let _ = state.clone().act(&operator.adjoint());
    let space = HilbertSpace::factor(0u8);
    let state = state.at(&space).unwrap();
    let operator = operator.at(&space).unwrap();
    let _ = state.clone().scaled(Opaque(4)).act(&operator).unwrap();
}

// No Num, Float, Copy, PartialEq, Neg, Div or RealPart. Real-part extraction is
// deliberately absent: scalar quotients do not need a real projection.
#[derive(Clone)]
struct Weight(f64);
impl Add for Weight {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}
impl Mul for Weight {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0)
    }
}
impl Zero for Weight {
    fn zero() -> Self {
        Self(0.0)
    }
    fn is_zero(&self) -> bool {
        self.0 == 0.0
    }
}
impl One for Weight {
    fn one() -> Self {
        Self(1.0)
    }
}
impl Conjugate for Weight {
    type Output = Self;
    fn conjugate(&self) -> Self {
        self.clone()
    }
}
impl Divide for Weight {
    type Output = Self;
    fn divide(self, rhs: Self) -> Self {
        Self(self.0 / rhs.0)
    }
}

#[test]
fn expectation_quotients_require_neither_real_part_nor_float_or_general_comparison() {
    let state =
        LocalState::<Primitive, Primitive, Weight>::primitive(Primitive).scaled(Weight(2.0));
    let op = LocalOperator::primitive(Primitive).scaled(Weight(3.0));
    assert_eq!(
        state.expectation(&op).divide(state.inner_product(&state)).0,
        3.0
    );
    let space = HilbertSpace::factor(0u8);
    let state = state.at(&space).unwrap();
    let op = op.at(&space).unwrap();
    assert_eq!(
        state
            .expectation(&op)
            .unwrap()
            .divide(state.inner_product(&state).unwrap())
            .0,
        3.0
    );
}

macro_rules! complex_regressions {
    ($name:ident, $real:ty, $large:expr, $small:expr, $tol:expr) => {
        #[test]
        fn $name() {
            // These self-inner-products are representable, but squaring them
            // during complex division would overflow or underflow.
            for scale in [$large, $small] {
                let state =
                    LocalState::<Primitive, Primitive, Complex<$real>>::primitive(Primitive)
                        .scaled(Complex::new(scale, 0.0));
                let op = LocalOperator::primitive(Primitive);
                let value = state.expectation(&op).divide(state.inner_product(&state));
                assert!((value.re - 1.0).abs() < $tol && value.im == 0.0);
                let space = HilbertSpace::factor(0u8);
                let state = state.at(&space).unwrap();
                let value = state
                    .expectation(&op.at(&space).unwrap())
                    .unwrap()
                    .divide(state.inner_product(&state).unwrap());
                assert!((value.re - 1.0).abs() < $tol && value.im == 0.0);
            }
            let quotient = Complex::<$real>::new(3.0, 4.0).divide(Complex::new(1.0, 2.0));
            assert!((quotient.re - 2.2).abs() < $tol);
            assert!((quotient.im + 0.4).abs() < $tol);
        }
    };
}
complex_regressions!(
    complex32_division_preserves_extreme_real_denominators,
    f32,
    1e15_f32,
    1e-15_f32,
    1e-5
);
complex_regressions!(
    complex64_division_preserves_extreme_real_denominators,
    f64,
    1e100_f64,
    1e-100_f64,
    1e-12
);
