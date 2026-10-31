//! Spin-1 dipoles, quadrupoles, and direct multiplication rules.

use crate::coefficient::{four, spin_scale, sqrt_two, three, two};
use num_complex::Complex;
use num_traits::{Float, One, Zero};
use qmb_core::{
    Adjoint, Dimension, DimensionOverflow, Identity, LinearCombination, Product, Trace,
};

/// Fifteen primitives spanning the nine-dimensional spin-1 operator space.
///
/// Cartesian components are `S = c J`, with physical spin `J` and `c = sqrt(3/2)`.
/// The nine Hermitian elements have `Tr(A† A)/3 = 1`.
/// `Sp = Sx + i Sy`, `Sm = Sx - i Sy`; no spherical-tensor phases are applied.
/// Primitive equality compares variants, not operator values after expansion.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpinOneOperator {
    /// The identity.
    I,
    /// The normalized x component `sqrt(3/2) Jx`.
    Sx,
    /// The normalized y component `sqrt(3/2) Jy`.
    Sy,
    /// The normalized z component, with eigenvalues `sqrt(3/2), 0, -sqrt(3/2)`.
    Sz,
    /// The raising operator `Sx + i Sy`.
    Sp,
    /// The lowering operator `Sx - i Sy`.
    Sm,
    /// The Cartesian quadrupole `(Sx² - Sy²)/c`.
    Dxy,
    /// The traceless diagonal quadrupole `sqrt(2) (Sz² - I)`.
    Q0,
    /// The normalized quadrupole `{Sx,Sy}/c`.
    Qxy,
    /// The normalized quadrupole `{Sy,Sz}/c`.
    Qyz,
    /// The normalized quadrupole `{Sz,Sx}/c`.
    Qzx,
    /// The quadrupole `{Sz, Sp}/c = Qzx + i Qyz`.
    Q1p,
    /// The quadrupole `{Sz, Sm}/c = Qzx - i Qyz`.
    Q1m,
    /// The double raising operator `Sp²/c = Dxy + i Qxy`.
    Q2p,
    /// The double lowering operator `Sm²/c = Dxy - i Qxy`.
    Q2m,
}

impl<R: Float> Product<Complex<R>> for SpinOneOperator {
    type Term = Self;
    fn product(&self, rhs: &Self) -> LinearCombination<Self, Complex<R>> {
        use LinearCombination as Lc;
        use SpinOneOperator::*;
        let real = |value| Complex::new(value, R::zero());
        let imaginary = |value| Complex::new(R::zero(), value);
        // Direct rules: analytically zero terms are omitted; each result keeps basis order.
        match (*self, *rhs) {
            (I, op) | (op, I) => Lc::primitive(op),
            (Sx, Sx) | (Qyz, Qyz) => {
                Lc::primitive(I)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
            }
            (Sx, Sy) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sx, Sz) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sx, Dxy) | (Sy, Qxy) | (Qzx, Sz) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qyz)
            }
            (Sx, Q0) => {
                Lc::term(real(-sqrt_two::<R>() / four::<R>()), Sx)
                    + Lc::term(
                        imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()),
                        Qyz,
                    )
            }
            (Sx, Qxy) | (Qyz, Sz) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sx, Qyz) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
            }
            (Sx, Qzx) | (Qyz, Sy) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sx, Sp) => {
                Lc::primitive(I)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sx, Sm) => {
                Lc::primitive(I)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sx, Q1p) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sx, Q1m) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sx, Q2p) | (Qzx, Q2p) | (Sm, Dxy) | (Q1m, Dxy) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sx, Q2m) | (Sp, Dxy) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sy, Sx) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sy, Sy) | (Qzx, Qzx) => {
                Lc::primitive(I)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
            }
            (Sy, Sz) | (Qyz, Dxy) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qyz)
            }
            (Sy, Dxy) => {
                Lc::term(real(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sy, Q0) => {
                Lc::term(real(-sqrt_two::<R>() / four::<R>()), Sy)
                    + Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Qzx)
            }
            (Sy, Qyz) | (Qzx, Sx) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sy, Qzx) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
            }
            (Sy, Sp) => {
                Lc::term(imaginary(R::one()), I)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sy, Sm) => {
                Lc::term(imaginary(-R::one()), I)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sy, Q1p) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sy, Q1m) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sy, Q2p) | (Qyz, Q2p) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sy, Q2m) | (Q1p, Qxy) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sz, Sx) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sz, Sy) | (Dxy, Qyz) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qyz)
            }
            (Sz, Sz) | (Dxy, Dxy) | (Qxy, Qxy) => {
                Lc::primitive(I) + Lc::term(real(sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Sz, Dxy) => Lc::term(imaginary(spin_scale::<R>()), Qxy),
            (Sz, Q0) | (Q0, Sz) => Lc::term(real(sqrt_two::<R>() / two::<R>()), Sz),
            (Sz, Qxy) => Lc::term(imaginary(-spin_scale::<R>()), Dxy),
            (Sz, Qyz) | (Qxy, Sx) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sz, Qzx) | (Dxy, Sx) | (Qxy, Sy) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qyz)
            }
            (Sz, Sp) | (Sz, Q1p) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sp)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Q1p)
            }
            (Sz, Sm) => {
                Lc::term(real(-spin_scale::<R>() / two::<R>()), Sm)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Q1m)
            }
            (Sz, Q1m) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sm)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Q1m)
            }
            (Sz, Q2p) | (Sp, Sp) | (Q1p, Sp) => Lc::term(real(spin_scale::<R>()), Q2p),
            (Q0, Q2p) | (Q2p, Q0) => Lc::term(real(sqrt_two::<R>() / two::<R>()), Q2p),
            (Sz, Q2m) | (Q1m, Sm) | (Q1m, Q1m) => Lc::term(real(-spin_scale::<R>()), Q2m),
            (Dxy, Sy) => {
                Lc::term(real(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Dxy, Sz) => Lc::term(imaginary(-spin_scale::<R>()), Qxy),
            (Dxy, Q0) | (Q0, Dxy) => Lc::term(real(sqrt_two::<R>() / two::<R>()), Dxy),
            (Dxy, Qxy) => Lc::term(imaginary(spin_scale::<R>()), Sz),
            (Dxy, Qzx) | (Qxy, Qyz) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Dxy, Sp) | (Dxy, Q1p) | (Q2m, Sx) | (Q2m, Qzx) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Dxy, Sm) | (Q2p, Sx) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Dxy, Q1m) | (Q2p, Qzx) => {
                Lc::term(real(-spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Dxy, Q2p) | (Q2m, Dxy) => {
                Lc::primitive(I)
                    + Lc::term(real(-spin_scale::<R>()), Sz)
                    + Lc::term(real(sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Dxy, Q2m) | (Q2p, Dxy) => {
                Lc::primitive(I)
                    + Lc::term(real(spin_scale::<R>()), Sz)
                    + Lc::term(real(sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Q0, Sx) => {
                Lc::term(real(-sqrt_two::<R>() / four::<R>()), Sx)
                    + Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Qyz)
            }
            (Q0, Sy) => {
                Lc::term(real(-sqrt_two::<R>() / four::<R>()), Sy)
                    + Lc::term(
                        imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()),
                        Qzx,
                    )
            }
            (Q0, Q0) => Lc::primitive(I) + Lc::term(real(-sqrt_two::<R>() / two::<R>()), Q0),
            (Q0, Qxy) | (Qxy, Q0) => Lc::term(real(sqrt_two::<R>() / two::<R>()), Qxy),
            (Q0, Qyz) => {
                Lc::term(imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()), Sx)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Qyz)
            }
            (Q0, Qzx) => {
                Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Sy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Qzx)
            }
            (Q0, Sp) => {
                Lc::term(real(-sqrt_two::<R>() / four::<R>()), Sp)
                    + Lc::term(real(three::<R>() * sqrt_two::<R>() / four::<R>()), Q1p)
            }
            (Q0, Sm) => {
                Lc::term(real(-sqrt_two::<R>() / four::<R>()), Sm)
                    + Lc::term(real(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q1m)
            }
            (Q0, Q1p) => {
                Lc::term(real(three::<R>() * sqrt_two::<R>() / four::<R>()), Sp)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q1p)
            }
            (Q0, Q1m) => {
                Lc::term(real(-three::<R>() * sqrt_two::<R>() / four::<R>()), Sm)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q1m)
            }
            (Q0, Q2m) | (Q2m, Q0) => Lc::term(real(sqrt_two::<R>() / two::<R>()), Q2m),
            (Sm, Sm) | (Sm, Q1m) | (Q2m, Sz) => Lc::term(real(spin_scale::<R>()), Q2m),
            (Qxy, Sz) => Lc::term(imaginary(spin_scale::<R>()), Dxy),
            (Qxy, Dxy) => Lc::term(imaginary(-spin_scale::<R>()), Sz),
            (Qxy, Qzx) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qyz)
            }
            (Qxy, Sp) | (Qxy, Q1p) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Qxy, Sm) | (Q2p, Qyz) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Qxy, Q1m) | (Q2p, Sy) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Qxy, Q2p) => {
                Lc::term(imaginary(R::one()), I)
                    + Lc::term(imaginary(-spin_scale::<R>()), Sz)
                    + Lc::term(imaginary(sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Qxy, Q2m) => {
                Lc::term(imaginary(-R::one()), I)
                    + Lc::term(imaginary(-spin_scale::<R>()), Sz)
                    + Lc::term(imaginary(-sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Qyz, Sx) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
            }
            (Qyz, Q0) => {
                Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Sx)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Qyz)
            }
            (Qyz, Qxy) | (Qzx, Dxy) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Qyz, Qzx) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qyz, Sp) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qyz, Sm) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qyz, Q1p) => {
                Lc::term(imaginary(R::one()), I)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qyz, Q1m) => {
                Lc::term(imaginary(-R::one()), I)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qyz, Q2m) | (Sp, Qxy) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Qzx, Sy) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
            }
            (Qzx, Q0) => {
                Lc::term(imaginary(-three::<R>() * sqrt_two::<R>() / four::<R>()), Sy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Qzx)
            }
            (Qzx, Qxy) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qyz)
            }
            (Qzx, Qyz) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qzx, Sp) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qzx, Sm) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qzx, Q1p) => {
                Lc::primitive(I)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qzx, Q1m) => {
                Lc::primitive(I)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Qzx, Q2m) | (Q1p, Dxy) => {
                Lc::term(real(-spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sp, Sx) => {
                Lc::primitive(I)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sp, Sy) => {
                Lc::term(imaginary(R::one()), I)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sp, Sz) => {
                Lc::term(real(-spin_scale::<R>() / two::<R>()), Sp)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Q1p)
            }
            (Sp, Q0) => {
                Lc::term(real(-sqrt_two::<R>() / four::<R>()), Sp)
                    + Lc::term(real(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q1p)
            }
            (Sp, Qyz) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sp, Qzx) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sp, Sm) | (Q1p, Q1m) => {
                Lc::term(real(two::<R>()), I)
                    + Lc::term(real(spin_scale::<R>()), Sz)
                    + Lc::term(real(-sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Sp, Q1p) | (Q1p, Q1p) | (Q2p, Sz) => Lc::term(real(-spin_scale::<R>()), Q2p),
            (Sp, Q1m) | (Q1p, Sm) => {
                Lc::term(real(spin_scale::<R>()), Sz)
                    + Lc::term(real(three::<R>() * sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Sp, Q2p)
            | (Sm, Q2m)
            | (Q1p, Q2p)
            | (Q1m, Q2m)
            | (Q2p, Sp)
            | (Q2p, Q1p)
            | (Q2p, Q2p)
            | (Q2m, Sm)
            | (Q2m, Q1m)
            | (Q2m, Q2m) => Lc::zero(),
            (Sp, Q2m) => {
                Lc::term(real(spin_scale::<R>()), Sm) + Lc::term(real(spin_scale::<R>()), Q1m)
            }
            (Sm, Sx) => {
                Lc::primitive(I)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sm, Sy) => {
                Lc::term(imaginary(-R::one()), I)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sm, Sz) | (Q1m, Sz) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sm)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Q1m)
            }
            (Sm, Q0) => {
                Lc::term(real(-sqrt_two::<R>() / four::<R>()), Sm)
                    + Lc::term(real(three::<R>() * sqrt_two::<R>() / four::<R>()), Q1m)
            }
            (Sm, Qxy) | (Q1m, Qxy) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Sm, Qyz) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sm, Qzx) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Sm, Sp) | (Q1m, Q1p) => {
                Lc::term(real(two::<R>()), I)
                    + Lc::term(real(-spin_scale::<R>()), Sz)
                    + Lc::term(real(-sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Sm, Q1p) | (Q1m, Sp) => {
                Lc::term(real(spin_scale::<R>()), Sz)
                    + Lc::term(real(-three::<R>() * sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Sm, Q2p) | (Q1m, Q2p) => {
                Lc::term(real(spin_scale::<R>()), Sp) + Lc::term(real(-spin_scale::<R>()), Q1p)
            }
            (Q1p, Sx) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Q1p, Sy) => {
                Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Q1p, Sz) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sp)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Q1p)
            }
            (Q1p, Q0) => {
                Lc::term(real(-three::<R>() * sqrt_two::<R>() / four::<R>()), Sp)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q1p)
            }
            (Q1p, Qyz) => {
                Lc::term(imaginary(R::one()), I)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Q1p, Qzx) => {
                Lc::primitive(I)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Q1p, Q2m) => {
                Lc::term(real(-spin_scale::<R>()), Sm) + Lc::term(real(-spin_scale::<R>()), Q1m)
            }
            (Q1m, Sx) => {
                Lc::term(real(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Q1m, Sy) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(three::<R>() * sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Q1m, Q0) => {
                Lc::term(real(three::<R>() * sqrt_two::<R>() / four::<R>()), Sm)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q1m)
            }
            (Q1m, Qyz) => {
                Lc::term(imaginary(-R::one()), I)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(imaginary(sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Q1m, Qzx) => {
                Lc::primitive(I)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Sz)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Dxy)
                    + Lc::term(real(-sqrt_two::<R>() / four::<R>()), Q0)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qxy)
            }
            (Q2p, Qxy) => {
                Lc::term(imaginary(R::one()), I)
                    + Lc::term(imaginary(spin_scale::<R>()), Sz)
                    + Lc::term(imaginary(sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Q2p, Sm) => {
                Lc::term(real(spin_scale::<R>()), Sp) + Lc::term(real(spin_scale::<R>()), Q1p)
            }
            (Q2p, Q1m) => {
                Lc::term(real(-spin_scale::<R>()), Sp) + Lc::term(real(-spin_scale::<R>()), Q1p)
            }
            (Q2p, Q2m) => {
                Lc::term(real(two::<R>()), I)
                    + Lc::term(real(two::<R>() * spin_scale::<R>()), Sz)
                    + Lc::term(real(sqrt_two::<R>()), Q0)
            }
            (Q2m, Sy) | (Q2m, Qyz) => {
                Lc::term(imaginary(-spin_scale::<R>() / two::<R>()), Sx)
                    + Lc::term(real(-spin_scale::<R>() / two::<R>()), Sy)
                    + Lc::term(real(spin_scale::<R>() / two::<R>()), Qyz)
                    + Lc::term(imaginary(spin_scale::<R>() / two::<R>()), Qzx)
            }
            (Q2m, Qxy) => {
                Lc::term(imaginary(-R::one()), I)
                    + Lc::term(imaginary(spin_scale::<R>()), Sz)
                    + Lc::term(imaginary(-sqrt_two::<R>() / two::<R>()), Q0)
            }
            (Q2m, Sp) | (Q2m, Q1p) => {
                Lc::term(real(spin_scale::<R>()), Sm) + Lc::term(real(-spin_scale::<R>()), Q1m)
            }
            (Q2m, Q2p) => {
                Lc::term(real(two::<R>()), I)
                    + Lc::term(real(-two::<R>() * spin_scale::<R>()), Sz)
                    + Lc::term(real(sqrt_two::<R>()), Q0)
            }
        }
    }
}

impl<C: One> Adjoint<C> for SpinOneOperator {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, C> {
        use SpinOneOperator::*;
        LinearCombination::primitive(match self {
            Sp => Sm,
            Sm => Sp,
            Q1p => Q1m,
            Q1m => Q1p,
            Q2p => Q2m,
            Q2m => Q2p,
            I | Sx | Sy | Sz | Dxy | Q0 | Qxy | Qyz | Qzx => *self,
        })
    }
}
impl<C: Zero + One> Trace<C> for SpinOneOperator {
    fn trace(&self) -> C {
        if *self == Self::I {
            C::one() + C::one() + C::one()
        } else {
            C::zero()
        }
    }
}
impl Identity for SpinOneOperator {
    fn identity() -> Self {
        Self::I
    }
}
impl Dimension for SpinOneOperator {
    fn dim() -> Result<core::num::NonZeroUsize, DimensionOverflow> {
        Ok(core::num::NonZeroUsize::new(9).unwrap())
    }
}
impl qmb_core::HomogeneousParity for SpinOneOperator {
    fn parity(&self) -> qmb_core::Parity {
        qmb_core::Parity::Even
    }
}
impl qmb_core::ParityOperator for SpinOneOperator {
    fn parity_operator() -> Self {
        Self::I
    }
}
