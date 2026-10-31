mod support;

use core::ops::{Add, Div, Mul, Neg};
use num_complex::Complex;
use num_traits::One;
use qmb_core::{
    Adjoint, Basis, HilbertSpace, Identity, ImaginaryUnit, InnerProduct, LinearCombination as Lc,
    Trace,
};
use qmb_spin_one::{
    CartesianBasis, LadderBasis, SpinOne, SpinOneOperator as O, SpinOneState as S, SzBasis,
};
use support::*;

// Exact arithmetic is test-owned. The outer coefficient deliberately implements
// neither Clone, Copy, Zero, Float, Num, NumCast, PartialEq, Debug nor Conjugate.
#[derive(Clone, Copy)]
struct Rational(i64, i64);
impl Rational {
    fn new(n: i64, d: i64) -> Self {
        assert_ne!(d, 0);
        let (mut a, mut b) = (n.abs(), d.abs());
        while b != 0 {
            (a, b) = (b, a % b);
        }
        let sign = if d < 0 { -1 } else { 1 };
        Self(n / a * sign, d / a * sign)
    }
    fn value(self) -> f64 {
        self.0 as f64 / self.1 as f64
    }
}
impl Add for Rational {
    type Output = Self;
    fn add(self, r: Self) -> Self {
        Self::new(self.0 * r.1 + r.0 * self.1, self.1 * r.1)
    }
}
impl Neg for Rational {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0, self.1)
    }
}
impl Mul for Rational {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        Self::new(self.0 * r.0, self.1 * r.1)
    }
}
impl Div for Rational {
    type Output = Self;
    fn div(self, r: Self) -> Self {
        Self::new(self.0 * r.1, self.1 * r.0)
    }
}
struct Exact {
    re: Rational,
    im: Rational,
}
impl Exact {
    fn integers(re: i64, im: i64) -> Self {
        Self {
            re: Rational(re, 1),
            im: Rational(im, 1),
        }
    }
    fn value(self) -> Complex<f64> {
        Complex::new(self.re.value(), self.im.value())
    }
}
impl Add for Exact {
    type Output = Self;
    fn add(self, r: Self) -> Self {
        Self {
            re: self.re + r.re,
            im: self.im + r.im,
        }
    }
}
impl Neg for Exact {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}
impl Mul for Exact {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        Self {
            re: self.re * r.re + -(self.im * r.im),
            im: self.re * r.im + self.im * r.re,
        }
    }
}
impl Div for Exact {
    type Output = Self;
    fn div(self, r: Self) -> Self {
        let norm = r.re * r.re + r.im * r.im;
        Self {
            re: (self.re * r.re + self.im * r.im) / norm,
            im: (self.im * r.re + -(self.re * r.im)) / norm,
        }
    }
}
impl One for Exact {
    fn one() -> Self {
        Self::integers(1, 0)
    }
}
impl ImaginaryUnit for Exact {
    fn imaginary_unit() -> Self {
        Self::integers(0, 1)
    }
}
fn exact_value(lc: Lc<O, Exact>) -> Matrix {
    lc.into_iter().fold(vec![0.0.into(); 9], |sum, t| {
        add(&sum, &scale(&matrix(t.operand), t.coefficient.value()))
    })
}

#[test]
fn basis_conversions_accept_nonclone_exact_coefficients() {
    for a in OPS {
        values::<f64>(
            &exact_value(<CartesianBasis as Basis<O, Exact>>::expand(
                &CartesianBasis,
                &a,
            )),
            &matrix(a),
        );
        values::<f64>(
            &exact_value(<LadderBasis as Basis<O, Exact>>::expand(&LadderBasis, &a)),
            &matrix(a),
        );
        assert_eq!(
            <CartesianBasis as Basis<O, Exact>>::elements(&CartesianBasis).count(),
            9
        );
        assert_eq!(
            <LadderBasis as Basis<O, Exact>>::elements(&LadderBasis).count(),
            9
        );
    }
}

struct OneOnly(u8);
impl Mul for OneOnly {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        Self(self.0 * r.0)
    }
}
impl One for OneOnly {
    fn one() -> Self {
        Self(1)
    }
}

// No addition or division: Cartesian expansion emits independent terms rather
// than requiring arithmetic that belongs to a later coefficient aggregation.
struct SignedComplex(i32, i32);
impl Mul for SignedComplex {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        Self(self.0 * r.0 - self.1 * r.1, self.0 * r.1 + self.1 * r.0)
    }
}
impl One for SignedComplex {
    fn one() -> Self {
        Self(1, 0)
    }
}
impl Neg for SignedComplex {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0, -self.1)
    }
}
impl ImaginaryUnit for SignedComplex {
    fn imaginary_unit() -> Self {
        Self(0, 1)
    }
}

#[test]
fn individual_capabilities_keep_their_small_bounds() {
    for a in OPS {
        let term = <O as Adjoint<OneOnly>>::adjoint(&a)
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(term.coefficient.0, 1);
        assert_eq!(<O as Trace<u32>>::trace(&a), if a == O::I { 3 } else { 0 });
        let expanded = <CartesianBasis as Basis<O, SignedComplex>>::expand(&CartesianBasis, &a);
        let actual = expanded.into_iter().fold(vec![0.0.into(); 9], |sum, t| {
            add(
                &sum,
                &scale(
                    &matrix(t.operand),
                    Complex::new(f64::from(t.coefficient.0), f64::from(t.coefficient.1)),
                ),
            )
        });
        values::<f64>(&actual, &matrix(a));
    }
    for s in STATES {
        let term = <SzBasis as Basis<S, OneOnly>>::expand(&SzBasis, &s)
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(term.coefficient.0, 1);
        assert_eq!(term.operand, s);
        for t in STATES {
            assert_eq!(
                <S as InnerProduct<u32>>::inner_product(&s, &t),
                u32::from(s == t)
            );
        }
    }
}

#[test]
fn storage_identity_placement_and_delayed_operations_need_no_numbers() {
    struct Opaque;
    type Spin = SpinOne<Opaque>;
    let psi = Spin::state(S::Zero);
    let op = Spin::operator(O::identity());
    let _ = psi.act(&(op.clone() * &op.adjoint()));
    let _ = qmb_core::LocalOperator::<O, Opaque>::identity();
    let a = HilbertSpace::factor("a");
    let b = HilbertSpace::factor("b");
    let psi = Spin::product_state(
        &qmb_core::HilbertSpace::tensor_all([(a).clone(), (b).clone()]).unwrap(),
        [S::Plus, S::Minus],
    )
    .unwrap();
    let _ = psi.act(&op.at(&a).unwrap()).unwrap();
    let _ = Spin::product_state(
        &qmb_core::HilbertSpace::tensor_all_g([(a).clone(), (b).clone()]).unwrap(),
        [S::Zero, S::Minus],
    )
    .unwrap();
}
