//! Spin-1 dipoles, quadrupoles, and direct multiplication rules.

use crate::coefficient::{
    four_thirds, half, imaginary, sixth, third, three_halves, two, two_thirds,
};
use core::ops::{Add, Div, Neg};
use num_traits::{One, Zero};
use qmb_core::{
    Adjoint, Dimension, DimensionOverflow, Identity, ImaginaryUnit, LinearCombination, Product,
    Trace,
};

/// Fifteen primitives spanning the nine-dimensional spin-1 operator space.
///
/// With `hbar = 1`, `SPlus = Sx + i Sy` and `SMinus = Sx - i Sy`.
/// The quadrupoles have no additional normalization or spherical-tensor phases.
/// Primitive equality compares variants, not operator values after expansion.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpinOneOperator {
    /// The identity.
    I,
    /// The x component of spin.
    Sx,
    /// The y component of spin.
    Sy,
    /// The z component, with eigenvalues `1, 0, -1`.
    Sz,
    /// The raising operator `Sx + i Sy`.
    SPlus,
    /// The lowering operator `Sx - i Sy`.
    SMinus,
    /// The Cartesian quadrupole `Sx² - Sy²`.
    Dxy,
    /// The traceless diagonal quadrupole `3 Sz² - 2 I`.
    Q0,
    /// The anticommutator `Sx Sy + Sy Sx`.
    Qxy,
    /// The anticommutator `Sy Sz + Sz Sy`.
    Qyz,
    /// The anticommutator `Sz Sx + Sx Sz`.
    Qzx,
    /// The quadrupole `{Sz, SPlus} = Qzx + i Qyz`.
    Q1Plus,
    /// The quadrupole `{Sz, SMinus} = Qzx - i Qyz`.
    Q1Minus,
    /// The double raising operator `SPlus² = Dxy + i Qxy`.
    Q2Plus,
    /// The double lowering operator `SMinus² = Dxy - i Qxy`.
    Q2Minus,
}

impl<C> Product<C> for SpinOneOperator
where
    C: One + Add<Output = C> + Div<Output = C> + Neg<Output = C> + ImaginaryUnit,
{
    type Term = Self;
    fn product(&self, rhs: &Self) -> LinearCombination<Self, C> {
        use LinearCombination as Lc;
        use SpinOneOperator::*;
        // Each branch emits only mathematically nonzero terms, in basis order.
        // The match is exhaustive; coefficients are constructed directly in C.
        match (*self, *rhs) {
            (I, op) | (op, I) => Lc::primitive(op),
            (Sx, Sx) | (Qyz, Qyz) => {
                Lc::term(two_thirds(), I) + Lc::term(half(), Dxy) + Lc::term(-sixth::<C>(), Q0)
            }
            (Sx, Sy) => Lc::term(imaginary(half()), Sz) + Lc::term(half(), Qxy),
            (Sx, Sz) => Lc::term(imaginary(-half::<C>()), Sy) + Lc::term(half(), Qzx),
            (Sx, Dxy) | (Sy, Qxy) | (Qzx, Sz) => {
                Lc::term(half(), Sx) + Lc::term(imaginary(-half::<C>()), Qyz)
            }
            (Sx, Q0) => Lc::term(-half::<C>(), Sx) + Lc::term(imaginary(-three_halves::<C>()), Qyz),
            (Sx, Qxy) | (Qyz, Sz) => Lc::term(half(), Sy) + Lc::term(imaginary(half()), Qzx),
            (Sx, Qyz) => Lc::term(imaginary(half()), Dxy) + Lc::term(imaginary(half()), Q0),
            (Sx, Qzx) | (Qyz, Sy) => Lc::term(half(), Sz) + Lc::term(imaginary(-half::<C>()), Qxy),
            (Sx, SPlus) => {
                Lc::term(two_thirds(), I)
                    + Lc::term(-half::<C>(), Sz)
                    + Lc::term(half(), Dxy)
                    + Lc::term(-sixth::<C>(), Q0)
                    + Lc::term(imaginary(half()), Qxy)
            }
            (Sx, SMinus) => {
                Lc::term(two_thirds(), I)
                    + Lc::term(half(), Sz)
                    + Lc::term(half(), Dxy)
                    + Lc::term(-sixth::<C>(), Q0)
                    + Lc::term(imaginary(-half::<C>()), Qxy)
            }
            (Sx, Q1Plus) => {
                Lc::term(half(), Sz)
                    + Lc::term(-half::<C>(), Dxy)
                    + Lc::term(-half::<C>(), Q0)
                    + Lc::term(imaginary(-half::<C>()), Qxy)
            }
            (Sx, Q1Minus) => {
                Lc::term(half(), Sz)
                    + Lc::term(half(), Dxy)
                    + Lc::term(half(), Q0)
                    + Lc::term(imaginary(-half::<C>()), Qxy)
            }
            (Sx, Q2Plus) | (Qzx, Q2Plus) | (SMinus, Dxy) | (Q1Minus, Dxy) => {
                Lc::term(half(), Sx)
                    + Lc::term(imaginary(half()), Sy)
                    + Lc::term(imaginary(-half::<C>()), Qyz)
                    + Lc::term(-half::<C>(), Qzx)
            }
            (Sx, Q2Minus) | (SPlus, Dxy) => {
                Lc::term(half(), Sx)
                    + Lc::term(imaginary(-half::<C>()), Sy)
                    + Lc::term(imaginary(-half::<C>()), Qyz)
                    + Lc::term(half(), Qzx)
            }
            (Sy, Sx) => Lc::term(imaginary(-half::<C>()), Sz) + Lc::term(half(), Qxy),
            (Sy, Sy) | (Qzx, Qzx) => {
                Lc::term(two_thirds(), I)
                    + Lc::term(-half::<C>(), Dxy)
                    + Lc::term(-sixth::<C>(), Q0)
            }
            (Sy, Sz) | (Qyz, Dxy) => Lc::term(imaginary(half()), Sx) + Lc::term(half(), Qyz),
            (Sy, Dxy) => Lc::term(-half::<C>(), Sy) + Lc::term(imaginary(-half::<C>()), Qzx),
            (Sy, Q0) => Lc::term(-half::<C>(), Sy) + Lc::term(imaginary(three_halves()), Qzx),
            (Sy, Qyz) | (Qzx, Sx) => Lc::term(half(), Sz) + Lc::term(imaginary(half()), Qxy),
            (Sy, Qzx) => Lc::term(imaginary(half()), Dxy) + Lc::term(imaginary(-half::<C>()), Q0),
            (Sy, SPlus) => {
                Lc::term(imaginary(two_thirds()), I)
                    + Lc::term(imaginary(-half::<C>()), Sz)
                    + Lc::term(imaginary(-half::<C>()), Dxy)
                    + Lc::term(imaginary(-sixth::<C>()), Q0)
                    + Lc::term(half(), Qxy)
            }
            (Sy, SMinus) => {
                Lc::term(imaginary(-two_thirds::<C>()), I)
                    + Lc::term(imaginary(-half::<C>()), Sz)
                    + Lc::term(imaginary(half()), Dxy)
                    + Lc::term(imaginary(sixth()), Q0)
                    + Lc::term(half(), Qxy)
            }
            (Sy, Q1Plus) => {
                Lc::term(imaginary(half()), Sz)
                    + Lc::term(imaginary(half()), Dxy)
                    + Lc::term(imaginary(-half::<C>()), Q0)
                    + Lc::term(-half::<C>(), Qxy)
            }
            (Sy, Q1Minus) => {
                Lc::term(imaginary(-half::<C>()), Sz)
                    + Lc::term(imaginary(half()), Dxy)
                    + Lc::term(imaginary(-half::<C>()), Q0)
                    + Lc::term(half(), Qxy)
            }
            (Sy, Q2Plus) | (Qyz, Q2Plus) => {
                Lc::term(imaginary(half()), Sx)
                    + Lc::term(-half::<C>(), Sy)
                    + Lc::term(half(), Qyz)
                    + Lc::term(imaginary(-half::<C>()), Qzx)
            }
            (Sy, Q2Minus) | (Q1Plus, Qxy) => {
                Lc::term(imaginary(-half::<C>()), Sx)
                    + Lc::term(-half::<C>(), Sy)
                    + Lc::term(-half::<C>(), Qyz)
                    + Lc::term(imaginary(-half::<C>()), Qzx)
            }
            (Sz, Sx) => Lc::term(imaginary(half()), Sy) + Lc::term(half(), Qzx),
            (Sz, Sy) | (Dxy, Qyz) => Lc::term(imaginary(-half::<C>()), Sx) + Lc::term(half(), Qyz),
            (Sz, Sz) | (Dxy, Dxy) | (Qxy, Qxy) => Lc::term(two_thirds(), I) + Lc::term(third(), Q0),
            (Sz, Dxy) => Lc::term(imaginary(C::one()), Qxy),
            (Sz, Q0) | (Q0, Sz) => Lc::primitive(Sz),
            (Sz, Qxy) => Lc::term(imaginary(-C::one()), Dxy),
            (Sz, Qyz) | (Qxy, Sx) => Lc::term(half(), Sy) + Lc::term(imaginary(-half::<C>()), Qzx),
            (Sz, Qzx) | (Dxy, Sx) | (Qxy, Sy) => {
                Lc::term(half(), Sx) + Lc::term(imaginary(half()), Qyz)
            }
            (Sz, SPlus) | (Sz, Q1Plus) => Lc::term(half(), SPlus) + Lc::term(half(), Q1Plus),
            (Sz, SMinus) => Lc::term(-half::<C>(), SMinus) + Lc::term(half(), Q1Minus),
            (Sz, Q1Minus) => Lc::term(half(), SMinus) + Lc::term(-half::<C>(), Q1Minus),
            (Sz, Q2Plus) | (Q0, Q2Plus) | (SPlus, SPlus) | (Q1Plus, SPlus) | (Q2Plus, Q0) => {
                Lc::primitive(Q2Plus)
            }
            (Sz, Q2Minus) | (Q1Minus, SMinus) | (Q1Minus, Q1Minus) => Lc::term(-C::one(), Q2Minus),
            (Dxy, Sy) => Lc::term(-half::<C>(), Sy) + Lc::term(imaginary(half()), Qzx),
            (Dxy, Sz) => Lc::term(imaginary(-C::one()), Qxy),
            (Dxy, Q0) | (Q0, Dxy) => Lc::primitive(Dxy),
            (Dxy, Qxy) => Lc::term(imaginary(C::one()), Sz),
            (Dxy, Qzx) | (Qxy, Qyz) => {
                Lc::term(imaginary(-half::<C>()), Sy) + Lc::term(-half::<C>(), Qzx)
            }
            (Dxy, SPlus) | (Dxy, Q1Plus) | (Q2Minus, Sx) | (Q2Minus, Qzx) => {
                Lc::term(half(), Sx)
                    + Lc::term(imaginary(-half::<C>()), Sy)
                    + Lc::term(imaginary(half()), Qyz)
                    + Lc::term(-half::<C>(), Qzx)
            }
            (Dxy, SMinus) | (Q2Plus, Sx) => {
                Lc::term(half(), Sx)
                    + Lc::term(imaginary(half()), Sy)
                    + Lc::term(imaginary(half()), Qyz)
                    + Lc::term(half(), Qzx)
            }
            (Dxy, Q1Minus) | (Q2Plus, Qzx) => {
                Lc::term(-half::<C>(), Sx)
                    + Lc::term(imaginary(-half::<C>()), Sy)
                    + Lc::term(imaginary(-half::<C>()), Qyz)
                    + Lc::term(-half::<C>(), Qzx)
            }
            (Dxy, Q2Plus) | (Q2Minus, Dxy) => {
                Lc::term(two_thirds(), I) + Lc::term(-C::one(), Sz) + Lc::term(third(), Q0)
            }
            (Dxy, Q2Minus) | (Q2Plus, Dxy) => {
                Lc::term(two_thirds(), I) + Lc::primitive(Sz) + Lc::term(third(), Q0)
            }
            (Q0, Sx) => Lc::term(-half::<C>(), Sx) + Lc::term(imaginary(three_halves()), Qyz),
            (Q0, Sy) => Lc::term(-half::<C>(), Sy) + Lc::term(imaginary(-three_halves::<C>()), Qzx),
            (Q0, Q0) => Lc::term(two(), I) + Lc::term(-C::one(), Q0),
            (Q0, Qxy) | (Qxy, Q0) => Lc::primitive(Qxy),
            (Q0, Qyz) => {
                Lc::term(imaginary(-three_halves::<C>()), Sx) + Lc::term(-half::<C>(), Qyz)
            }
            (Q0, Qzx) => Lc::term(imaginary(three_halves()), Sy) + Lc::term(-half::<C>(), Qzx),
            (Q0, SPlus) => Lc::term(-half::<C>(), SPlus) + Lc::term(three_halves(), Q1Plus),
            (Q0, SMinus) => {
                Lc::term(-half::<C>(), SMinus) + Lc::term(-three_halves::<C>(), Q1Minus)
            }
            (Q0, Q1Plus) => Lc::term(three_halves(), SPlus) + Lc::term(-half::<C>(), Q1Plus),
            (Q0, Q1Minus) => {
                Lc::term(-three_halves::<C>(), SMinus) + Lc::term(-half::<C>(), Q1Minus)
            }
            (Q0, Q2Minus)
            | (SMinus, SMinus)
            | (SMinus, Q1Minus)
            | (Q2Minus, Sz)
            | (Q2Minus, Q0) => Lc::primitive(Q2Minus),
            (Qxy, Sz) => Lc::term(imaginary(C::one()), Dxy),
            (Qxy, Dxy) => Lc::term(imaginary(-C::one()), Sz),
            (Qxy, Qzx) => Lc::term(imaginary(half()), Sx) + Lc::term(-half::<C>(), Qyz),
            (Qxy, SPlus) | (Qxy, Q1Plus) => {
                Lc::term(imaginary(half()), Sx)
                    + Lc::term(half(), Sy)
                    + Lc::term(-half::<C>(), Qyz)
                    + Lc::term(imaginary(-half::<C>()), Qzx)
            }
            (Qxy, SMinus) | (Q2Plus, Qyz) => {
                Lc::term(imaginary(-half::<C>()), Sx)
                    + Lc::term(half(), Sy)
                    + Lc::term(half(), Qyz)
                    + Lc::term(imaginary(-half::<C>()), Qzx)
            }
            (Qxy, Q1Minus) | (Q2Plus, Sy) => {
                Lc::term(imaginary(half()), Sx)
                    + Lc::term(-half::<C>(), Sy)
                    + Lc::term(-half::<C>(), Qyz)
                    + Lc::term(imaginary(half()), Qzx)
            }
            (Qxy, Q2Plus) => {
                Lc::term(imaginary(two_thirds()), I)
                    + Lc::term(imaginary(-C::one()), Sz)
                    + Lc::term(imaginary(third()), Q0)
            }
            (Qxy, Q2Minus) => {
                Lc::term(imaginary(-two_thirds::<C>()), I)
                    + Lc::term(imaginary(-C::one()), Sz)
                    + Lc::term(imaginary(-third::<C>()), Q0)
            }
            (Qyz, Sx) => {
                Lc::term(imaginary(-half::<C>()), Dxy) + Lc::term(imaginary(-half::<C>()), Q0)
            }
            (Qyz, Q0) => Lc::term(imaginary(three_halves()), Sx) + Lc::term(-half::<C>(), Qyz),
            (Qyz, Qxy) | (Qzx, Dxy) => {
                Lc::term(imaginary(half()), Sy) + Lc::term(-half::<C>(), Qzx)
            }
            (Qyz, Qzx) => Lc::term(imaginary(-half::<C>()), Sz) + Lc::term(-half::<C>(), Qxy),
            (Qyz, SPlus) => {
                Lc::term(imaginary(half()), Sz)
                    + Lc::term(imaginary(-half::<C>()), Dxy)
                    + Lc::term(imaginary(-half::<C>()), Q0)
                    + Lc::term(half(), Qxy)
            }
            (Qyz, SMinus) => {
                Lc::term(imaginary(-half::<C>()), Sz)
                    + Lc::term(imaginary(-half::<C>()), Dxy)
                    + Lc::term(imaginary(-half::<C>()), Q0)
                    + Lc::term(-half::<C>(), Qxy)
            }
            (Qyz, Q1Plus) => {
                Lc::term(imaginary(two_thirds()), I)
                    + Lc::term(imaginary(-half::<C>()), Sz)
                    + Lc::term(imaginary(half()), Dxy)
                    + Lc::term(imaginary(-sixth::<C>()), Q0)
                    + Lc::term(-half::<C>(), Qxy)
            }
            (Qyz, Q1Minus) => {
                Lc::term(imaginary(-two_thirds::<C>()), I)
                    + Lc::term(imaginary(-half::<C>()), Sz)
                    + Lc::term(imaginary(-half::<C>()), Dxy)
                    + Lc::term(imaginary(sixth()), Q0)
                    + Lc::term(-half::<C>(), Qxy)
            }
            (Qyz, Q2Minus) | (SPlus, Qxy) => {
                Lc::term(imaginary(half()), Sx)
                    + Lc::term(half(), Sy)
                    + Lc::term(half(), Qyz)
                    + Lc::term(imaginary(half()), Qzx)
            }
            (Qzx, Sy) => Lc::term(imaginary(-half::<C>()), Dxy) + Lc::term(imaginary(half()), Q0),
            (Qzx, Q0) => {
                Lc::term(imaginary(-three_halves::<C>()), Sy) + Lc::term(-half::<C>(), Qzx)
            }
            (Qzx, Qxy) => Lc::term(imaginary(-half::<C>()), Sx) + Lc::term(-half::<C>(), Qyz),
            (Qzx, Qyz) => Lc::term(imaginary(half()), Sz) + Lc::term(-half::<C>(), Qxy),
            (Qzx, SPlus) => {
                Lc::term(half(), Sz)
                    + Lc::term(half(), Dxy)
                    + Lc::term(-half::<C>(), Q0)
                    + Lc::term(imaginary(half()), Qxy)
            }
            (Qzx, SMinus) => {
                Lc::term(half(), Sz)
                    + Lc::term(-half::<C>(), Dxy)
                    + Lc::term(half(), Q0)
                    + Lc::term(imaginary(half()), Qxy)
            }
            (Qzx, Q1Plus) => {
                Lc::term(two_thirds(), I)
                    + Lc::term(-half::<C>(), Sz)
                    + Lc::term(-half::<C>(), Dxy)
                    + Lc::term(-sixth::<C>(), Q0)
                    + Lc::term(imaginary(-half::<C>()), Qxy)
            }
            (Qzx, Q1Minus) => {
                Lc::term(two_thirds(), I)
                    + Lc::term(half(), Sz)
                    + Lc::term(-half::<C>(), Dxy)
                    + Lc::term(-sixth::<C>(), Q0)
                    + Lc::term(imaginary(half()), Qxy)
            }
            (Qzx, Q2Minus) | (Q1Plus, Dxy) => {
                Lc::term(-half::<C>(), Sx)
                    + Lc::term(imaginary(half()), Sy)
                    + Lc::term(imaginary(half()), Qyz)
                    + Lc::term(-half::<C>(), Qzx)
            }
            (SPlus, Sx) => {
                Lc::term(two_thirds(), I)
                    + Lc::term(half(), Sz)
                    + Lc::term(half(), Dxy)
                    + Lc::term(-sixth::<C>(), Q0)
                    + Lc::term(imaginary(half()), Qxy)
            }
            (SPlus, Sy) => {
                Lc::term(imaginary(two_thirds()), I)
                    + Lc::term(imaginary(half()), Sz)
                    + Lc::term(imaginary(-half::<C>()), Dxy)
                    + Lc::term(imaginary(-sixth::<C>()), Q0)
                    + Lc::term(half(), Qxy)
            }
            (SPlus, Sz) => Lc::term(-half::<C>(), SPlus) + Lc::term(half(), Q1Plus),
            (SPlus, Q0) => Lc::term(-half::<C>(), SPlus) + Lc::term(-three_halves::<C>(), Q1Plus),
            (SPlus, Qyz) => {
                Lc::term(imaginary(half()), Sz)
                    + Lc::term(imaginary(half()), Dxy)
                    + Lc::term(imaginary(half()), Q0)
                    + Lc::term(-half::<C>(), Qxy)
            }
            (SPlus, Qzx) => {
                Lc::term(half(), Sz)
                    + Lc::term(-half::<C>(), Dxy)
                    + Lc::term(half(), Q0)
                    + Lc::term(imaginary(-half::<C>()), Qxy)
            }
            (SPlus, SMinus) | (Q1Plus, Q1Minus) => {
                Lc::term(four_thirds(), I) + Lc::primitive(Sz) + Lc::term(-third::<C>(), Q0)
            }
            (SPlus, Q1Plus) | (Q1Plus, Q1Plus) | (Q2Plus, Sz) => Lc::term(-C::one(), Q2Plus),
            (SPlus, Q1Minus) | (Q1Plus, SMinus) => Lc::primitive(Sz) + Lc::primitive(Q0),
            (SPlus, Q2Plus)
            | (SMinus, Q2Minus)
            | (Q1Plus, Q2Plus)
            | (Q1Minus, Q2Minus)
            | (Q2Plus, SPlus)
            | (Q2Plus, Q1Plus)
            | (Q2Plus, Q2Plus)
            | (Q2Minus, SMinus)
            | (Q2Minus, Q1Minus)
            | (Q2Minus, Q2Minus) => Lc::zero(),
            (SPlus, Q2Minus) => Lc::primitive(SMinus) + Lc::primitive(Q1Minus),
            (SMinus, Sx) => {
                Lc::term(two_thirds(), I)
                    + Lc::term(-half::<C>(), Sz)
                    + Lc::term(half(), Dxy)
                    + Lc::term(-sixth::<C>(), Q0)
                    + Lc::term(imaginary(-half::<C>()), Qxy)
            }
            (SMinus, Sy) => {
                Lc::term(imaginary(-two_thirds::<C>()), I)
                    + Lc::term(imaginary(half()), Sz)
                    + Lc::term(imaginary(half()), Dxy)
                    + Lc::term(imaginary(sixth()), Q0)
                    + Lc::term(half(), Qxy)
            }
            (SMinus, Sz) | (Q1Minus, Sz) => Lc::term(half(), SMinus) + Lc::term(half(), Q1Minus),
            (SMinus, Q0) => Lc::term(-half::<C>(), SMinus) + Lc::term(three_halves(), Q1Minus),
            (SMinus, Qxy) | (Q1Minus, Qxy) => {
                Lc::term(imaginary(-half::<C>()), Sx)
                    + Lc::term(half(), Sy)
                    + Lc::term(-half::<C>(), Qyz)
                    + Lc::term(imaginary(half()), Qzx)
            }
            (SMinus, Qyz) => {
                Lc::term(imaginary(-half::<C>()), Sz)
                    + Lc::term(imaginary(half()), Dxy)
                    + Lc::term(imaginary(half()), Q0)
                    + Lc::term(half(), Qxy)
            }
            (SMinus, Qzx) => {
                Lc::term(half(), Sz)
                    + Lc::term(half(), Dxy)
                    + Lc::term(-half::<C>(), Q0)
                    + Lc::term(imaginary(-half::<C>()), Qxy)
            }
            (SMinus, SPlus) | (Q1Minus, Q1Plus) => {
                Lc::term(four_thirds(), I) + Lc::term(-C::one(), Sz) + Lc::term(-third::<C>(), Q0)
            }
            (SMinus, Q1Plus) | (Q1Minus, SPlus) => Lc::primitive(Sz) + Lc::term(-C::one(), Q0),
            (SMinus, Q2Plus) | (Q1Minus, Q2Plus) => {
                Lc::primitive(SPlus) + Lc::term(-C::one(), Q1Plus)
            }
            (Q1Plus, Sx) => {
                Lc::term(half(), Sz)
                    + Lc::term(half(), Dxy)
                    + Lc::term(half(), Q0)
                    + Lc::term(imaginary(half()), Qxy)
            }
            (Q1Plus, Sy) => {
                Lc::term(imaginary(half()), Sz)
                    + Lc::term(imaginary(-half::<C>()), Dxy)
                    + Lc::term(imaginary(half()), Q0)
                    + Lc::term(half(), Qxy)
            }
            (Q1Plus, Sz) => Lc::term(half(), SPlus) + Lc::term(-half::<C>(), Q1Plus),
            (Q1Plus, Q0) => Lc::term(-three_halves::<C>(), SPlus) + Lc::term(-half::<C>(), Q1Plus),
            (Q1Plus, Qyz) => {
                Lc::term(imaginary(two_thirds()), I)
                    + Lc::term(imaginary(half()), Sz)
                    + Lc::term(imaginary(half()), Dxy)
                    + Lc::term(imaginary(-sixth::<C>()), Q0)
                    + Lc::term(-half::<C>(), Qxy)
            }
            (Q1Plus, Qzx) => {
                Lc::term(two_thirds(), I)
                    + Lc::term(half(), Sz)
                    + Lc::term(-half::<C>(), Dxy)
                    + Lc::term(-sixth::<C>(), Q0)
                    + Lc::term(imaginary(-half::<C>()), Qxy)
            }
            (Q1Plus, Q2Minus) => Lc::term(-C::one(), SMinus) + Lc::term(-C::one(), Q1Minus),
            (Q1Minus, Sx) => {
                Lc::term(half(), Sz)
                    + Lc::term(-half::<C>(), Dxy)
                    + Lc::term(-half::<C>(), Q0)
                    + Lc::term(imaginary(half()), Qxy)
            }
            (Q1Minus, Sy) => {
                Lc::term(imaginary(-half::<C>()), Sz)
                    + Lc::term(imaginary(-half::<C>()), Dxy)
                    + Lc::term(imaginary(half()), Q0)
                    + Lc::term(-half::<C>(), Qxy)
            }
            (Q1Minus, Q0) => Lc::term(three_halves(), SMinus) + Lc::term(-half::<C>(), Q1Minus),
            (Q1Minus, Qyz) => {
                Lc::term(imaginary(-two_thirds::<C>()), I)
                    + Lc::term(imaginary(half()), Sz)
                    + Lc::term(imaginary(-half::<C>()), Dxy)
                    + Lc::term(imaginary(sixth()), Q0)
                    + Lc::term(-half::<C>(), Qxy)
            }
            (Q1Minus, Qzx) => {
                Lc::term(two_thirds(), I)
                    + Lc::term(-half::<C>(), Sz)
                    + Lc::term(-half::<C>(), Dxy)
                    + Lc::term(-sixth::<C>(), Q0)
                    + Lc::term(imaginary(half()), Qxy)
            }
            (Q2Plus, Qxy) => {
                Lc::term(imaginary(two_thirds()), I)
                    + Lc::term(imaginary(C::one()), Sz)
                    + Lc::term(imaginary(third()), Q0)
            }
            (Q2Plus, SMinus) => Lc::primitive(SPlus) + Lc::primitive(Q1Plus),
            (Q2Plus, Q1Minus) => Lc::term(-C::one(), SPlus) + Lc::term(-C::one(), Q1Plus),
            (Q2Plus, Q2Minus) => {
                Lc::term(four_thirds(), I) + Lc::term(two(), Sz) + Lc::term(two_thirds(), Q0)
            }
            (Q2Minus, Sy) | (Q2Minus, Qyz) => {
                Lc::term(imaginary(-half::<C>()), Sx)
                    + Lc::term(-half::<C>(), Sy)
                    + Lc::term(half(), Qyz)
                    + Lc::term(imaginary(half()), Qzx)
            }
            (Q2Minus, Qxy) => {
                Lc::term(imaginary(-two_thirds::<C>()), I)
                    + Lc::term(imaginary(C::one()), Sz)
                    + Lc::term(imaginary(-third::<C>()), Q0)
            }
            (Q2Minus, SPlus) | (Q2Minus, Q1Plus) => {
                Lc::primitive(SMinus) + Lc::term(-C::one(), Q1Minus)
            }
            (Q2Minus, Q2Plus) => {
                Lc::term(four_thirds(), I) + Lc::term(-two::<C>(), Sz) + Lc::term(two_thirds(), Q0)
            }
        }
    }
}

impl<C: One> Adjoint<C> for SpinOneOperator {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, C> {
        use SpinOneOperator::*;
        LinearCombination::primitive(match self {
            SPlus => SMinus,
            SMinus => SPlus,
            Q1Plus => Q1Minus,
            Q1Minus => Q1Plus,
            Q2Plus => Q2Minus,
            Q2Minus => Q2Plus,
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
