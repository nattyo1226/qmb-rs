//! Ordered, finite linear combinations returned by algebra rules.

use crate::{
    coefficient::multiply,
    expr::parts::{IntoOperands, ScaleExpr, SumExpr},
};
use core::{
    ops::{Add, Mul},
    slice,
};
use num_traits::One;

/// A coefficient and an operand; distinct from an operation's associated `Term`.
pub type Term<P, C> = ScaleExpr<P, C>;

/// An ordered sum preserving duplicates, zero coefficients, and input order.
///
/// Empty and single-term containers do not allocate. Storage and element mapping
/// require no arithmetic, cloning, or comparison capabilities.
#[derive(Clone, Debug, PartialEq)]
pub struct LinearCombination<P, C> {
    terms: SumExpr<ScaleExpr<P, C>>,
}

impl<P, C> LinearCombination<P, C> {
    /// Constructs the empty sum.
    pub fn zero() -> Self {
        Self {
            terms: SumExpr::Empty,
        }
    }

    /// Constructs one term, retaining even a zero coefficient.
    pub fn term(coefficient: C, operand: P) -> Self {
        Self {
            terms: SumExpr::One(ScaleExpr {
                coefficient,
                operand,
            }),
        }
    }

    /// Collects terms in order, without allocation for at most one term.
    pub fn from_terms(terms: impl IntoIterator<Item = Term<P, C>>) -> Self {
        Self {
            terms: terms.into_iter().collect(),
        }
    }

    /// Borrows the stored terms in order.
    pub fn terms(&self) -> slice::Iter<'_, Term<P, C>> {
        self.terms.iter()
    }

    /// Moves the terms out, taking over any existing allocation.
    pub fn into_terms(self) -> IntoTerms<P, C> {
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
    pub fn map<Q>(self, mut f: impl FnMut(P) -> Q) -> LinearCombination<Q, C> {
        LinearCombination {
            terms: self.terms.map(|term| ScaleExpr {
                coefficient: term.coefficient,
                operand: f(term.operand),
            }),
        }
    }

    /// Constructs one operand with coefficient one.
    pub fn primitive(operand: P) -> Self
    where
        C: One,
    {
        Self::term(C::one(), operand)
    }

    /// Multiplies coefficients, without cloning operands or removing zero terms.
    /// Converts the input once even for an empty sum; the input need not be Clone.
    pub fn scaled(self, coefficient: impl Into<C>) -> Self
    where
        C: Clone + Mul<Output = C>,
    {
        let coefficient = coefficient.into();
        Self {
            terms: self.terms.map(|term| ScaleExpr {
                coefficient: multiply(term.coefficient, coefficient.clone()),
                operand: term.operand,
            }),
        }
    }
}

impl<P, C> Default for LinearCombination<P, C> {
    fn default() -> Self {
        Self::zero()
    }
}
impl<P, C> Add for LinearCombination<P, C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            terms: self.terms.concat(rhs.terms),
        }
    }
}
impl<P, C> FromIterator<Term<P, C>> for LinearCombination<P, C> {
    fn from_iter<I: IntoIterator<Item = Term<P, C>>>(iter: I) -> Self {
        Self::from_terms(iter)
    }
}
impl<P, C> IntoIterator for LinearCombination<P, C> {
    type Item = Term<P, C>;
    type IntoIter = IntoTerms<P, C>;
    fn into_iter(self) -> Self::IntoIter {
        self.into_terms()
    }
}
impl<'a, P, C> IntoIterator for &'a LinearCombination<P, C> {
    type Item = &'a Term<P, C>;
    type IntoIter = slice::Iter<'a, Term<P, C>>;
    fn into_iter(self) -> Self::IntoIter {
        self.terms()
    }
}

/// Owning, exact-size, double-ended iteration over linear-combination terms.
pub type IntoTerms<P, C> = IntoOperands<Term<P, C>>;
