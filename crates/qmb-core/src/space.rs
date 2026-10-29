//! Immutable ordinary and graded tensor spaces, generic over an equality-comparable index.

use alloc::{sync::Arc, vec, vec::Vec};
use core::{iter::FusedIterator, slice};

use crate::{Dimension, DimensionOverflow, expr::parts::SumExpr};

/// Composition rule for an ordered tensor product.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TensorKind {
    /// Ordinary tensor product, with commuting distinct factors.
    Ordinary,
    /// Z₂-graded tensor product, with Koszul signs for odd factors.
    Graded,
}

/// Errors constructing a tensor space.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum SpaceError {
    /// The same basic index occurs more than once, including inside children.
    #[error("duplicate factor index")]
    DuplicateIndex,
    /// Different effective composition rules remain after removing units.
    #[error("mixed tensor kinds")]
    MixedTensorKinds,
}

enum Node<I> {
    Factor(I),
    Tensor(TensorKind, SumExpr<HilbertSpace<I>>),
}

/// A shared factor or ordered, possibly nested tensor space.
///
/// Local primitive types and dimensions belong to expressions and their query
/// inputs. Index identity is supplied by I: Eq; cloning never clones an index.
///
/// ```
/// use core::num::NonZeroUsize;
/// use qmb_core::HilbertSpace;
/// let a = HilbertSpace::factor("a");
/// let b = HilbertSpace::factor("b");
/// let h = a.tensor(&b)?;
/// assert_eq!(h.dim_with(NonZeroUsize::new(2).unwrap())?, 4);
/// assert!(a.is_subsystem_of(&h));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct HilbertSpace<I> {
    root: Option<Arc<Node<I>>>,
    // Saturated leaf count: 0 (unit), 1 (factor), 2 (multiple leaves).
    // Structural metadata only; no dimensions or evaluation results are cached.
    leaf_class: u8,
    effective_kind: Option<TensorKind>,
}

impl<I> Clone for HilbertSpace<I> {
    fn clone(&self) -> Self {
        Self {
            root: self.root.clone(),
            leaf_class: self.leaf_class,
            effective_kind: self.effective_kind,
        }
    }
}

impl<I> HilbertSpace<I> {
    fn new(node: Node<I>) -> Self {
        let (leaf_class, effective_kind) = match &node {
            Node::Factor(_) => (1, None),
            Node::Tensor(kind, children) => {
                let mut nonempty = children.iter().filter(|c| c.leaf_class != 0);
                match (nonempty.next(), nonempty.next()) {
                    (None, _) => (0, None),
                    (Some(child), None) => (child.leaf_class, child.effective_kind),
                    (Some(_), Some(_)) => (2, Some(*kind)),
                }
            }
        };
        Self {
            root: Some(Arc::new(node)),
            leaf_class,
            effective_kind,
        }
    }

    fn node(&self) -> &Node<I> {
        self.root.as_deref().expect("live space handle")
    }

    /// Creates a leaf from an already issued index, without generating an index.
    pub fn factor(index: I) -> Self {
        Self::new(Node::Factor(index))
    }

    /// Creates the empty tensor space, of dimension one.
    pub fn unit() -> Self {
        Self::new(Node::Tensor(TensorKind::Ordinary, SumExpr::Empty))
    }

    /// Returns the index of a leaf, or `None` for a tensor (even a singleton).
    pub fn factor_index(&self) -> Option<&I> {
        match self.node() {
            Node::Factor(i) => Some(i),
            Node::Tensor(_, _) => None,
        }
    }

    /// Returns the stored tensor kind, including empty and singleton wrappers.
    pub fn kind(&self) -> Option<TensorKind> {
        match self.node() {
            Node::Factor(_) => None,
            Node::Tensor(kind, _) => Some(*kind),
        }
    }

    /// Returns the composition rule after ignoring units and singleton wrappers.
    pub fn effective_kind(&self) -> Option<TensorKind> {
        self.effective_kind
    }

    /// Borrows the immediate children, preserving the stored grouping.
    pub fn children(&self) -> slice::Iter<'_, Self> {
        match self.node() {
            Node::Factor(_) => [].iter(),
            Node::Tensor(_, s) => s.iter(),
        }
    }

    /// Lazily borrows all leaf indices in tensor order, using a work stack.
    pub fn indices(&self) -> Indices<'_, I> {
        Indices {
            leaves: self.leaves(),
        }
    }

    /// Lazily borrows the actual leaves in tensor order, including a leaf input.
    /// Empty tensors yield no leaves; traversal uses an explicit work stack.
    pub fn leaves(&self) -> Leaves<'_, I> {
        Leaves {
            pending: vec![self],
        }
    }

    /// Computes the tensor-space dimension using primitive P's span per factor.
    ///
    /// Use a state primitive for a state-space dimension, and an operator
    /// primitive for an operator-space dimension. The empty tensor has dimension
    /// one even when P's local dimension cannot be represented.
    pub fn dim<P: Dimension>(&self) -> Result<usize, DimensionOverflow> {
        let mut indices = self.indices();
        if indices.next().is_none() {
            return Ok(1);
        }
        let local = P::dim()?.get();
        indices.try_fold(local, |d, _| d.checked_mul(local).ok_or(DimensionOverflow))
    }

    /// Computes the dimension for an explicitly chosen uniform local space.
    ///
    /// This query does not associate a definition or primitive type with indices.
    pub fn dim_with(&self, local: core::num::NonZeroUsize) -> Result<usize, DimensionOverflow> {
        self.indices().try_fold(1usize, |d, _| {
            d.checked_mul(local.get()).ok_or(DimensionOverflow)
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

impl<I: Eq> HilbertSpace<I> {
    /// Constructs a tensor, retaining children and grouping and rejecting duplicates.
    pub fn tensor_all(children: impl IntoIterator<Item = Self>) -> Result<Self, SpaceError> {
        Self::tensor_all_with(TensorKind::Ordinary, children)
    }

    /// Constructs a graded tensor, preserving grouping.
    pub fn tensor_all_g(children: impl IntoIterator<Item = Self>) -> Result<Self, SpaceError> {
        Self::tensor_all_with(TensorKind::Graded, children)
    }

    /// Constructs a tensor with an explicit composition rule.
    pub fn tensor_all_with(
        kind: TensorKind,
        children: impl IntoIterator<Item = Self>,
    ) -> Result<Self, SpaceError> {
        let children: SumExpr<Self> = children.into_iter().collect();
        if children
            .iter()
            .filter(|c| c.leaf_class != 0)
            .take(2)
            .count()
            == 2
            && children
                .iter()
                .any(|c| c.effective_kind.is_some_and(|k| k != kind))
        {
            return Err(SpaceError::MixedTensorKinds);
        }
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
        Ok(Self::new(Node::Tensor(kind, children)))
    }

    /// Constructs a binary tensor with this space's factors on the left.
    pub fn tensor(&self, rhs: &Self) -> Result<Self, SpaceError> {
        self.tensor_with(TensorKind::Ordinary, rhs)
    }
    /// Constructs a binary graded tensor.
    pub fn tensor_g(&self, rhs: &Self) -> Result<Self, SpaceError> {
        self.tensor_with(TensorKind::Graded, rhs)
    }
    /// Constructs a binary tensor with an explicit rule.
    pub fn tensor_with(&self, kind: TensorKind, rhs: &Self) -> Result<Self, SpaceError> {
        Self::tensor_all_with(kind, [self.clone(), rhs.clone()])
    }

    /// Tests equality of ordered factors, independently of grouping and sharing.
    pub fn is_compatible_with(&self, rhs: &Self) -> bool {
        self.shares_node_with(rhs)
            || (self.effective_kind == rhs.effective_kind && self.indices().eq(rhs.indices()))
    }

    /// Tests whether this space's factors form an ordered subsequence of `ambient`.
    pub fn is_subsystem_of(&self, ambient: &Self) -> bool {
        if self.leaf_class == 2 && self.effective_kind != ambient.effective_kind {
            return false;
        }
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
pub struct Indices<'a, I> {
    leaves: Leaves<'a, I>,
}

impl<'a, I> Iterator for Indices<'a, I> {
    type Item = &'a I;
    fn next(&mut self) -> Option<Self::Item> {
        self.leaves.next().and_then(HilbertSpace::factor_index)
    }
}
impl<I> FusedIterator for Indices<'_, I> {}

/// A borrowing, depth-first iterator over the actual leaf space handles.
pub struct Leaves<'a, I> {
    pending: Vec<&'a HilbertSpace<I>>,
}

impl<'a, I> Iterator for Leaves<'a, I> {
    type Item = &'a HilbertSpace<I>;
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(space) = self.pending.pop() {
            match space.node() {
                Node::Factor(_) => return Some(space),
                Node::Tensor(_, children) => self.pending.extend(children.iter().rev()),
            }
        }
        None
    }
}
impl<I> FusedIterator for Leaves<'_, I> {}

impl<I> Drop for HilbertSpace<I> {
    fn drop(&mut self) {
        let mut first = self.root.take();
        let mut pending = Vec::new();
        while let Some(root) = first.take().or_else(|| pending.pop()) {
            if let Some(Node::Tensor(_, children)) = Arc::into_inner(root) {
                for mut child in children {
                    if let Some(root) = child.root.take() {
                        pending.push(root);
                    }
                }
            }
        }
    }
}
