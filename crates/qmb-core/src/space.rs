//! Immutable ordinary tensor spaces, generic over an equality-comparable index.

use alloc::{sync::Arc, vec, vec::Vec};
use core::{iter::FusedIterator, marker::PhantomData, slice};

use crate::{DimensionOverflow, LocalStructure, expr::parts::SumExpr};

/// Errors constructing a tensor space.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum SpaceError {
    /// The same basic index occurs more than once, including inside children.
    #[error("duplicate factor index")]
    DuplicateIndex,
}

enum Node<I, L: LocalStructure> {
    Factor(I),
    Tensor(SumExpr<HilbertSpace<I, L>>),
}

/// A shared factor or ordered, possibly nested ordinary tensor space.
///
/// `L` represents one fixed local definition. Index identity is supplied by
/// `I: Eq`; generation, hashing and ordering are not required by this crate.
/// Cloning this handle never clones `I` or `L`.
///
/// ```
/// use core::num::NonZeroUsize;
/// use qmb_core::{HilbertSpace, LocalStructure};
/// struct Qubit;
/// impl LocalStructure for Qubit {
///     type State = ();
///     type Operator = ();
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(2).unwrap();
/// }
/// let a = HilbertSpace::<_, Qubit>::factor("a");
/// let b = HilbertSpace::factor("b");
/// let h = a.tensor_product(&b)?;
/// assert_eq!(h.dimension()?, 4);
/// assert!(a.is_subsystem_of(&h));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// Structures remain distinct even if their dimensions and primitive types agree:
///
/// ```compile_fail
/// use core::num::NonZeroUsize;
/// use qmb_core::{HilbertSpace, LocalStructure};
/// struct A;
/// struct B;
/// impl LocalStructure for A {
///     type State = ();
///     type Operator = ();
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
/// }
/// impl LocalStructure for B {
///     type State = ();
///     type Operator = ();
///     const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
/// }
/// let a = HilbertSpace::<u8, A>::factor(0);
/// let b = HilbertSpace::<u8, B>::factor(0);
/// let _ = a.tensor_product(&b);
/// ```
pub struct HilbertSpace<I, L: LocalStructure> {
    root: Option<Arc<Node<I, L>>>,
    metadata: PhantomData<fn() -> L>,
}

impl<I, L: LocalStructure> Clone for HilbertSpace<I, L> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
            metadata: PhantomData,
        }
    }
}

impl<I, L: LocalStructure> HilbertSpace<I, L> {
    fn new(node: Node<I, L>) -> Self {
        Self {
            root: Some(Arc::new(node)),
            metadata: PhantomData,
        }
    }

    fn node(&self) -> &Node<I, L> {
        self.root.as_deref().expect("live space handle")
    }

    /// Creates a leaf from an already issued index, without generating an index.
    pub fn factor(index: I) -> Self {
        Self::new(Node::Factor(index))
    }

    /// Creates the empty tensor space, of dimension one.
    pub fn unit() -> Self {
        Self::new(Node::Tensor(SumExpr::Empty))
    }

    /// Returns the index of a leaf, or `None` for a tensor (even a singleton).
    pub fn factor_index(&self) -> Option<&I> {
        match self.node() {
            Node::Factor(i) => Some(i),
            Node::Tensor(_) => None,
        }
    }

    /// Borrows the immediate children, preserving the stored grouping.
    pub fn children(&self) -> slice::Iter<'_, Self> {
        match self.node() {
            Node::Factor(_) => [].iter(),
            Node::Tensor(s) => s.iter(),
        }
    }

    /// Lazily borrows all leaf indices in tensor order, using a work stack.
    pub fn indices(&self) -> Indices<'_, I, L> {
        Indices {
            pending: vec![self],
        }
    }

    /// Computes the dimension, failing only when this integer query overflows.
    pub fn dimension(&self) -> Result<usize, DimensionOverflow> {
        self.indices().try_fold(1usize, |d, _| {
            d.checked_mul(L::DIMENSION.get()).ok_or(DimensionOverflow)
        })
    }

    /// Whether two handles share the same immutable space node.
    pub fn shares_node_with(&self, rhs: &Self) -> bool {
        Arc::ptr_eq(
            self.root.as_ref().expect("live space"),
            rhs.root.as_ref().expect("live space"),
        )
    }
}

impl<I: Eq, L: LocalStructure> HilbertSpace<I, L> {
    /// Constructs a tensor, retaining children and grouping and rejecting duplicates.
    pub fn tensor(children: impl IntoIterator<Item = Self>) -> Result<Self, SpaceError> {
        let children: SumExpr<Self> = children.into_iter().collect();
        // Every existing child is valid. A sole child cannot introduce duplicates,
        // and preserving this case avoids quadratic construction of deep grouping.
        if let SumExpr::Many(parts) = &children {
            let mut seen: Vec<&I> = Vec::new();
            for index in parts.iter().flat_map(Self::indices) {
                if seen.contains(&index) {
                    return Err(SpaceError::DuplicateIndex);
                }
                seen.push(index);
            }
        }
        Ok(Self::new(Node::Tensor(children)))
    }

    /// Constructs a binary tensor with this space's factors on the left.
    pub fn tensor_product(&self, rhs: &Self) -> Result<Self, SpaceError> {
        Self::tensor([self.clone(), rhs.clone()])
    }

    /// Tests equality of ordered factors, independently of grouping and sharing.
    pub fn is_compatible_with(&self, rhs: &Self) -> bool {
        self.shares_node_with(rhs) || self.indices().eq(rhs.indices())
    }

    /// Tests whether this space's factors form an ordered subsequence of `ambient`.
    pub fn is_subsystem_of(&self, ambient: &Self) -> bool {
        let mut ambient = ambient.indices();
        self.indices()
            .all(|index| ambient.any(|other| index == other))
    }

    /// Tests whether this space contains an index, without inferring its structure.
    pub fn contains_index(&self, index: &I) -> bool {
        self.indices().any(|i| i == index)
    }
}

/// A borrowing, depth-first iterator over a space's leaf indices.
pub struct Indices<'a, I, L: LocalStructure> {
    pending: Vec<&'a HilbertSpace<I, L>>,
}

impl<'a, I, L: LocalStructure> Iterator for Indices<'a, I, L> {
    type Item = &'a I;
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(space) = self.pending.pop() {
            match space.node() {
                Node::Factor(index) => return Some(index),
                Node::Tensor(children) => self.pending.extend(children.iter().rev()),
            }
        }
        None
    }
}
impl<I, L: LocalStructure> FusedIterator for Indices<'_, I, L> {}

impl<I, L: LocalStructure> Drop for HilbertSpace<I, L> {
    fn drop(&mut self) {
        let mut first = self.root.take();
        let mut pending = Vec::new();
        while let Some(root) = first.take().or_else(|| pending.pop()) {
            if let Some(Node::Tensor(children)) = Arc::into_inner(root) {
                for mut child in children {
                    if let Some(root) = child.root.take() {
                        pending.push(root);
                    }
                }
            }
        }
    }
}
