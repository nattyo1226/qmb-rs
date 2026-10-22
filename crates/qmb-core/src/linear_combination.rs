//! Ordered, finite linear combinations returned by algebra rules.

use crate::{
    Complex,
    expr::parts::{IntoOperands, ScaleExpr, SumExpr},
};
use core::{ops::Add, slice};
use num_traits::{Num, One, Zero};

/// A coefficient and an operand; distinct from an operation's associated `Term`.
pub type Term<P, T> = ScaleExpr<P, T>;

/// An ordered sum preserving duplicates, zero coefficients, and input order.
///
/// Empty and single-term containers do not allocate. Storage and element mapping
/// require no arithmetic, cloning, or comparison capabilities.
#[derive(Clone, Debug, PartialEq)]
pub struct LinearCombination<P, T> {
    terms: SumExpr<ScaleExpr<P, T>>,
}

impl<P, T> LinearCombination<P, T> {
    /// Constructs the empty sum.
    pub fn zero() -> Self {
        Self {
            terms: SumExpr::Empty,
        }
    }

    /// Constructs one term, retaining even a zero coefficient.
    pub fn term(coefficient: Complex<T>, operand: P) -> Self {
        Self {
            terms: SumExpr::One(ScaleExpr {
                coefficient,
                operand,
            }),
        }
    }

    /// Collects terms in order, without allocation for at most one term.
    pub fn from_terms(terms: impl IntoIterator<Item = Term<P, T>>) -> Self {
        Self {
            terms: terms.into_iter().collect(),
        }
    }

    /// Borrows the stored terms in order.
    pub fn terms(&self) -> slice::Iter<'_, Term<P, T>> {
        self.terms.iter()
    }

    /// Moves the terms out, taking over any existing allocation.
    pub fn into_terms(self) -> IntoTerms<P, T> {
        self.terms.into_iter()
    }

    /// Returns the number of stored terms, including zero coefficients.
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Returns whether the term sequence is empty, not mathematical zero.
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    /// Transforms operands once in order, moving coefficients unchanged.
    pub fn map<Q>(self, mut f: impl FnMut(P) -> Q) -> LinearCombination<Q, T> {
        LinearCombination {
            terms: self.terms.map(|term| ScaleExpr {
                coefficient: term.coefficient,
                operand: f(term.operand),
            }),
        }
    }
}

impl<P, T: Zero + One> LinearCombination<P, T> {
    /// Constructs one operand with coefficient one.
    pub fn primitive(operand: P) -> Self {
        Self::term(Complex::new(T::one(), T::zero()), operand)
    }
}

impl<P, T: Clone + Num> LinearCombination<P, T> {
    /// Multiplies coefficients, without cloning operands or removing zero terms.
    pub fn scaled(self, coefficient: Complex<T>) -> Self {
        Self {
            terms: self.terms.map(|term| ScaleExpr {
                coefficient: term.coefficient * coefficient.clone(),
                operand: term.operand,
            }),
        }
    }
}

impl<P, T> Default for LinearCombination<P, T> {
    fn default() -> Self {
        Self::zero()
    }
}
impl<P, T> Add for LinearCombination<P, T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            terms: self.terms.concat(rhs.terms),
        }
    }
}
impl<P, T> FromIterator<Term<P, T>> for LinearCombination<P, T> {
    fn from_iter<I: IntoIterator<Item = Term<P, T>>>(iter: I) -> Self {
        Self::from_terms(iter)
    }
}
impl<P, T> IntoIterator for LinearCombination<P, T> {
    type Item = Term<P, T>;
    type IntoIter = IntoTerms<P, T>;
    fn into_iter(self) -> Self::IntoIter {
        self.into_terms()
    }
}
impl<'a, P, T> IntoIterator for &'a LinearCombination<P, T> {
    type Item = &'a Term<P, T>;
    type IntoIter = slice::Iter<'a, Term<P, T>>;
    fn into_iter(self) -> Self::IntoIter {
        self.terms()
    }
}

/// Owning, exact-size, double-ended iteration over linear-combination terms.
pub type IntoTerms<P, T> = IntoOperands<Term<P, T>>;
