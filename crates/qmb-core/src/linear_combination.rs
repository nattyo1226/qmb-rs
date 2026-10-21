//! Ordered, finite linear combinations returned by local rules.

use alloc::{vec, vec::Vec};
use core::{iter::FusedIterator, ops::Add, slice};

use num_complex::Complex;
use num_traits::{Num, One, Zero};

/// One coefficient and primitive, without any numeric bounds on their storage.
#[derive(Clone, Debug, PartialEq)]
pub struct Term<P, T> {
    /// The coefficient of this term, including an explicitly stored zero.
    pub coefficient: Complex<T>,
    /// The primitive multiplied by the coefficient.
    pub primitive: P,
}

/// An ordered sum of terms, preserving duplicates and zero coefficients.
///
/// Empty and single-term combinations need no heap allocation. Equality compares
/// the stored sequence, not mathematical values. No sorting, term merging, or
/// zero-coefficient removal is performed.
#[derive(Clone, Debug, PartialEq)]
pub struct LinearCombination<P, T> {
    repr: Repr<P, T>,
}

#[derive(Clone, Debug, PartialEq)]
enum Repr<P, T> {
    Zero,
    Term(Term<P, T>),
    Sum(Vec<Term<P, T>>),
}

impl<P, T> LinearCombination<P, T> {
    /// Returns the empty sum without requiring numeric operations.
    pub fn zero() -> Self {
        Self { repr: Repr::Zero }
    }

    /// Constructs one term, retaining even a zero coefficient.
    pub fn term(coefficient: Complex<T>, primitive: P) -> Self {
        Self {
            repr: Repr::Term(Term {
                coefficient,
                primitive,
            }),
        }
    }

    /// Collects terms in input order without simplifying them.
    ///
    /// No allocation occurs when the iterator yields at most one term.
    pub fn from_terms(terms: impl IntoIterator<Item = Term<P, T>>) -> Self {
        let mut terms = terms.into_iter();
        let Some(first) = terms.next() else {
            return Self::zero();
        };
        let Some(second) = terms.next() else {
            return Self {
                repr: Repr::Term(first),
            };
        };
        let mut sum = vec![first, second];
        sum.extend(terms);
        Self {
            repr: Repr::Sum(sum),
        }
    }

    /// Borrows terms in their stored order.
    pub fn terms(&self) -> slice::Iter<'_, Term<P, T>> {
        match &self.repr {
            Repr::Zero => [].iter(),
            Repr::Term(term) => slice::from_ref(term).iter(),
            Repr::Sum(terms) => terms.iter(),
        }
    }

    /// Takes ownership of the terms, without cloning primitives or coefficients.
    pub fn into_terms(self) -> IntoTerms<P, T> {
        let repr = match self.repr {
            Repr::Zero => IntoTermsRepr::Single(None.into_iter()),
            Repr::Term(term) => IntoTermsRepr::Single(Some(term).into_iter()),
            Repr::Sum(terms) => IntoTermsRepr::Sum(terms.into_iter()),
        };
        IntoTerms { repr }
    }

    /// Returns the stored term count, including zero coefficients.
    pub fn len(&self) -> usize {
        self.terms().len()
    }

    /// Returns whether the term sequence is empty, not whether its value is zero.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<P, T: Zero + One> LinearCombination<P, T> {
    /// Constructs one primitive with coefficient one.
    pub fn primitive(primitive: P) -> Self {
        Self::term(Complex::new(T::one(), T::zero()), primitive)
    }
}

impl<P, T: Clone + Num> LinearCombination<P, T> {
    /// Multiplies every coefficient, preserving term order and zero terms.
    ///
    /// This consumes the combination and never requires `P: Clone`.
    pub fn scaled(self, coefficient: Complex<T>) -> Self {
        match self.repr {
            Repr::Zero => Self::zero(),
            Repr::Term(term) => Self::term(term.coefficient * coefficient, term.primitive),
            Repr::Sum(mut terms) => {
                for term in &mut terms {
                    term.coefficient = term.coefficient.clone() * coefficient.clone();
                }
                Self {
                    repr: Repr::Sum(terms),
                }
            }
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

    /// Concatenates the left and right sequences, without arithmetic or cloning.
    fn add(self, rhs: Self) -> Self {
        let repr = match (self.repr, rhs.repr) {
            (Repr::Zero, repr) | (repr, Repr::Zero) => repr,
            (Repr::Term(left), Repr::Term(right)) => Repr::Sum(vec![left, right]),
            (Repr::Sum(mut left), Repr::Term(right)) => {
                left.push(right);
                Repr::Sum(left)
            }
            (Repr::Term(left), Repr::Sum(mut right)) => {
                right.insert(0, left);
                Repr::Sum(right)
            }
            (Repr::Sum(mut left), Repr::Sum(right)) => {
                left.extend(right);
                Repr::Sum(left)
            }
        };
        Self { repr }
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

/// An owning iterator over a finite linear combination.
///
/// Empty and single-term iterators do not allocate. Larger iterators take over
/// the combination's existing storage.
pub struct IntoTerms<P, T> {
    repr: IntoTermsRepr<P, T>,
}

enum IntoTermsRepr<P, T> {
    Single(core::option::IntoIter<Term<P, T>>),
    Sum(alloc::vec::IntoIter<Term<P, T>>),
}

impl<P, T> Iterator for IntoTerms<P, T> {
    type Item = Term<P, T>;

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.repr {
            IntoTermsRepr::Single(iter) => iter.next(),
            IntoTermsRepr::Sum(iter) => iter.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.len();
        (len, Some(len))
    }
}

impl<P, T> DoubleEndedIterator for IntoTerms<P, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match &mut self.repr {
            IntoTermsRepr::Single(iter) => iter.next_back(),
            IntoTermsRepr::Sum(iter) => iter.next_back(),
        }
    }
}

impl<P, T> ExactSizeIterator for IntoTerms<P, T> {
    fn len(&self) -> usize {
        match &self.repr {
            IntoTermsRepr::Single(iter) => iter.len(),
            IntoTermsRepr::Sum(iter) => iter.len(),
        }
    }
}

impl<P, T> FusedIterator for IntoTerms<P, T> {}
