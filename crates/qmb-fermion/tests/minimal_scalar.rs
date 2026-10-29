mod support;
use core::ops::{Add, Div, Mul, Neg};
use num_complex::Complex;
use num_traits::{One, Zero};
use qmb_core::{
    Act, Adjoint, Basis, HilbertSpace, Identity, ImaginaryUnit, InnerProduct,
    LinearCombination as Lc, Product, Trace,
};
use qmb_fermion::{
    Fermion, FermionOperator as O, FermionState as S, LadderBasis, MajoranaBasis, OccupationBasis,
};
use qmb_geometry::{Boundary, Chain};
use support::*;

// A small exact rational implementation used only by these capability tests.
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

// Deliberately no Clone, Copy, Float, Num, Debug, PartialEq, or conjugation.
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
impl Zero for Exact {
    fn zero() -> Self {
        Self::integers(0, 0)
    }
    fn is_zero(&self) -> bool {
        self.re.0 == 0 && self.im.0 == 0
    }
}
impl ImaginaryUnit for Exact {
    fn imaginary_unit() -> Self {
        Self::integers(0, 1)
    }
}

fn exact_matrix(lc: Lc<O, Exact>) -> Matrix {
    lc.into_iter().fold(vec![Complex::new(0., 0.); 4], |m, t| {
        add(&m, &scale(&matrix(t.operand), t.coefficient.value()))
    })
}
#[test]
fn all_primitive_rules_and_bases_work_without_clone_copy_float_or_conjugation() {
    for a in OPS {
        for b in OPS {
            let actual = exact_matrix(<O as Product<Exact>>::product(&a, &b));
            assert_eq!(actual, product(&matrix(a), &matrix(b)));
        }
        for (j, s) in STATES.into_iter().enumerate() {
            let mut v = [Complex::new(0., 0.); 2];
            for t in <S as Act<Exact>>::act(&s, &a) {
                v[usize::from(t.operand == S::Occupied)] += t.coefficient.value();
            }
            let m = matrix(a);
            assert_eq!(v, [m[j], m[2 + j]]);
        }
        assert_eq!(
            exact_matrix(<LadderBasis as Basis<O, Exact>>::expand(&LadderBasis, &a)),
            matrix(a)
        );
        assert_eq!(
            exact_matrix(<MajoranaBasis as Basis<O, Exact>>::expand(
                &MajoranaBasis,
                &a
            )),
            matrix(a)
        );
    }
    assert!(<O as Product<Exact>>::product(&O::Annihilate, &O::Annihilate).is_empty());
    assert!(<S as Act<Exact>>::act(&S::Empty, &O::Annihilate).is_empty());
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
#[test]
fn independent_capabilities_do_not_inherit_majorana_bounds() {
    for a in OPS {
        let term = <O as Adjoint<OneOnly>>::adjoint(&a)
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(term.coefficient.0, 1);
        let tr = <O as Trace<u32>>::trace(&a);
        assert_eq!(
            tr,
            if a == O::I {
                2
            } else if a == O::Number {
                1
            } else {
                0
            }
        );
    }
    for s in STATES {
        let expanded = <OccupationBasis as Basis<S, OneOnly>>::expand(&OccupationBasis, &s);
        assert_eq!(expanded.into_iter().next().unwrap().coefficient.0, 1);
        for t in STATES {
            assert_eq!(
                <S as InnerProduct<u32>>::inner_product(&s, &t),
                u32::from(s == t)
            );
        }
    }
    struct Opaque;
    let ket = Fermion::<Opaque>::state(S::Empty);
    let op = Fermion::<Opaque>::operator(O::identity());
    let _ = ket.act(&(op.clone() * &op.adjoint()));
    let a = HilbertSpace::factor(0);
    let b = HilbertSpace::factor(1);
    let psi = Fermion::<Opaque>::product_state(&a.tensor_g(&b).unwrap(), [S::Empty, S::Occupied])
        .unwrap();
    let _ = psi.act(&op.at(&a).unwrap()).unwrap();
}

// Models need coefficient cloning and negation, but no arithmetic or evaluator.
#[derive(Clone)]
struct ModelCoefficient(i32);
impl From<i32> for ModelCoefficient {
    fn from(n: i32) -> Self {
        Self(n)
    }
}
impl Neg for ModelCoefficient {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}
#[derive(PartialEq, Eq)]
struct Index(u8);
#[test]
fn hubbard_construction_needs_neither_primitive_evaluation_nor_index_cloning() {
    let site = HilbertSpace::tensor_all_g([
        HilbertSpace::factor(Index(9)),
        HilbertSpace::factor(Index(1)),
    ])
    .unwrap();
    let other = HilbertSpace::tensor_all_g([
        HilbertSpace::factor(Index(4)),
        HilbertSpace::factor(Index(2)),
    ])
    .unwrap();
    let chain = Chain::new(2, Boundary::Periodic);
    let sites = [site, other];
    for model in [
        qmb_fermion::models::FermionModels::<Fermion<ModelCoefficient>>::hubbard(
            &chain, &sites, 1, 4,
        ),
        qmb_fermion::models::FermionModels::<Fermion<ModelCoefficient>>::symmetric_hubbard(
            &chain, &sites, 1, 4,
        ),
    ] {
        assert!(
            model
                .unwrap()
                .space()
                .unwrap()
                .children()
                .next()
                .unwrap()
                .shares_node_with(&sites[0])
        );
    }
}
