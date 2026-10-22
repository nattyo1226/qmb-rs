//! Reusable syntax containers, independent of local structures and evaluation.

use alloc::{vec, vec::Vec};
use core::{iter::FusedIterator, slice};

use crate::Complex;

/// A coefficient and its operand, without numeric bounds on storage.
#[derive(Clone, Debug, PartialEq)]
pub struct ScaleExpr<E, T> {
    /// The coefficient, including explicitly stored zero and one.
    pub coefficient: Complex<T>,
    /// The scaled operand.
    pub operand: E,
}

/// An ordered sequence with inline storage for zero or one operand.
///
/// All variants represent their stored sequence, including an empty `Many`.
#[derive(Clone, Debug, PartialEq, Default)]
pub enum SumExpr<E> {
    /// The empty sequence.
    #[default]
    Empty,
    /// One operand, stored inline.
    One(E),
    /// An owned sequence of operands.
    Many(Vec<E>),
}

impl<E> SumExpr<E> {
    /// Borrows operands in their stored order.
    pub fn iter(&self) -> slice::Iter<'_, E> {
        self.as_slice().iter()
    }

    /// Borrows the contiguous stored sequence.
    pub fn as_slice(&self) -> &[E] {
        match self {
            Self::Empty => &[],
            Self::One(value) => slice::from_ref(value),
            Self::Many(values) => values,
        }
    }

    /// Returns the number of operands.
    pub fn len(&self) -> usize {
        self.as_slice().len()
    }

    /// Returns whether there are no operands.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Transforms each operand once, in order, without cloning.
    pub fn map<F>(self, f: impl FnMut(E) -> F) -> SumExpr<F> {
        self.into_iter().map(f).collect()
    }

    /// Concatenates sequences without arithmetic or cloning.
    pub fn concat(self, rhs: Self) -> Self {
        match (self, rhs) {
            (Self::Empty, rhs) => rhs,
            (lhs, Self::Empty) => lhs,
            (Self::One(lhs), Self::One(rhs)) => Self::Many(vec![lhs, rhs]),
            (Self::Many(mut lhs), Self::One(rhs)) => {
                lhs.push(rhs);
                Self::Many(lhs)
            }
            (Self::One(lhs), Self::Many(mut rhs)) => {
                rhs.insert(0, lhs);
                Self::Many(rhs)
            }
            (Self::Many(mut lhs), Self::Many(rhs)) => {
                lhs.extend(rhs);
                Self::Many(lhs)
            }
        }
    }
}

impl<E> FromIterator<E> for SumExpr<E> {
    fn from_iter<I: IntoIterator<Item = E>>(iter: I) -> Self {
        let mut iter = iter.into_iter();
        let Some(first) = iter.next() else {
            return Self::Empty;
        };
        let Some(second) = iter.next() else {
            return Self::One(first);
        };
        let mut values = vec![first, second];
        values.extend(iter);
        Self::Many(values)
    }
}

/// An owning iterator that preserves inline empty and single-element storage.
pub enum IntoOperands<E> {
    /// Zero or one operand.
    Inline(core::option::IntoIter<E>),
    /// Operands from an existing allocation.
    Many(alloc::vec::IntoIter<E>),
}

impl<E> IntoIterator for SumExpr<E> {
    type Item = E;
    type IntoIter = IntoOperands<E>;

    fn into_iter(self) -> Self::IntoIter {
        match self {
            Self::Empty => IntoOperands::Inline(None.into_iter()),
            Self::One(value) => IntoOperands::Inline(Some(value).into_iter()),
            Self::Many(values) => IntoOperands::Many(values.into_iter()),
        }
    }
}

impl<'a, E> IntoIterator for &'a SumExpr<E> {
    type Item = &'a E;
    type IntoIter = slice::Iter<'a, E>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<E> Iterator for IntoOperands<E> {
    type Item = E;

    fn next(&mut self) -> Option<E> {
        match self {
            Self::Inline(iter) => iter.next(),
            Self::Many(iter) => iter.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.len(), Some(self.len()))
    }
}

impl<E> DoubleEndedIterator for IntoOperands<E> {
    fn next_back(&mut self) -> Option<E> {
        match self {
            Self::Inline(iter) => iter.next_back(),
            Self::Many(iter) => iter.next_back(),
        }
    }
}

impl<E> ExactSizeIterator for IntoOperands<E> {
    fn len(&self) -> usize {
        match self {
            Self::Inline(iter) => iter.len(),
            Self::Many(iter) => iter.len(),
        }
    }
}

impl<E> FusedIterator for IntoOperands<E> {}

/// Ordered binary multiplication; the right operand acts first.
#[derive(Clone, Debug, PartialEq)]
pub struct ProductExpr<E> {
    /// Left operand.
    pub lhs: E,
    /// Right operand.
    pub rhs: E,
}

/// The left action of an operator on a state.
#[derive(Clone, Debug, PartialEq)]
pub struct ActExpr<S, O> {
    /// State receiving the action.
    pub state: S,
    /// Operator acting on the state.
    pub operator: O,
}

/// An unevaluated adjoint.
#[derive(Clone, Debug, PartialEq)]
pub struct AdjointExpr<E> {
    /// Operand whose adjoint is taken.
    pub operand: E,
}
