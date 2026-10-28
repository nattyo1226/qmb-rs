//! Pure chain geometry with borrowing connections to caller-owned spaces.

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

use qmb_core::HilbertSpace;

/// Failure to connect a chain to a site-ordered input.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum ChainError {
    /// The supplied number of sites differs from the chain's count.
    #[error("site count mismatch")]
    SiteCountMismatch,
}

/// Boundary condition for an undirected chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// Connect consecutive sites only.
    Open,
    /// Also connect the last site to the first, without loops or duplicate edges.
    Periodic,
}

/// Immutable site count and boundary condition, without indices or assignments.
///
/// ```
/// use qmb_geometry::{Boundary, Chain};
/// let chain = Chain::new(3, Boundary::Open);
/// let values = ["a", "b", "c"];
/// assert_eq!(chain.nn_site_pairs_with(&values)?.collect::<Vec<_>>(),
///            [(&"a", &"b"), (&"b", &"c")]);
/// # Ok::<(), qmb_geometry::ChainError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chain {
    site_count: usize,
    boundary: Boundary,
}

impl Chain {
    /// Constructs geometry without issuing indices; zero sites are allowed.
    pub fn new(site_count: usize, boundary: Boundary) -> Self {
        Self {
            site_count,
            boundary,
        }
    }

    /// Returns the number of sites, including any assigned empty sites.
    pub fn site_count(&self) -> usize {
        self.site_count
    }

    /// Returns the boundary condition.
    pub fn boundary(&self) -> Boundary {
        self.boundary
    }

    /// Enumerates normalized, unique site pairs in lexicographic order.
    /// Zero and one sites have no edges; two sites have one for either boundary.
    pub fn nn_site_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        let n = self.site_count;
        let closing = n > 2 && self.boundary == Boundary::Periodic;
        (0..n.saturating_sub(1)).flat_map(move |i| {
            core::iter::once((i, i + 1)).chain((closing && i == 0).then_some((0, n - 1)))
        })
    }

    /// Maps site pairs to borrowed elements, without cloning or inspecting them.
    /// Validates the exact input length before returning even an empty iterator.
    pub fn nn_site_pairs_with<'a, T>(
        &'a self,
        values: &'a [T],
    ) -> Result<impl Iterator<Item = (&'a T, &'a T)> + 'a, ChainError> {
        if values.len() != self.site_count {
            return Err(ChainError::SiteCountMismatch);
        }
        Ok(self.nn_site_pairs().map(|(i, j)| (&values[i], &values[j])))
    }

    /// Enumerates all leaf combinations across neighboring sites, left outermost.
    /// Empty sites retain their positions and contribute no leaf pairs.
    /// Index uniqueness is the responsibility of space or model construction.
    pub fn nn_pairs_with<'a, I>(
        &'a self,
        factors: &'a [HilbertSpace<I>],
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        Ok(self.nn_site_pairs_with(factors)?.flat_map(|(left, right)| {
            left.leaves()
                .flat_map(move |a| right.leaves().map(move |b| (a, b)))
        }))
    }

    /// Pairs corresponding leaf positions, stopping at the shorter site's end.
    pub fn nn_pairs_by_pos_with<'a, I>(
        &'a self,
        factors: &'a [HilbertSpace<I>],
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        Ok(self
            .nn_site_pairs_with(factors)?
            .flat_map(|(left, right)| left.leaves().zip(right.leaves())))
    }

    /// Reads the stored top-level partition as sites and borrows site pairs.
    /// A leaf is one site; an empty tensor is zero sites. Nested grouping stays intact.
    pub fn nn_site_pairs_on<'a, I>(
        &'a self,
        space: &'a HilbertSpace<I>,
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        self.nn_site_pairs_with(sites_on(space))
    }

    /// Reads top-level sites and delegates to the Cartesian leaf-pair connection.
    pub fn nn_pairs_on<'a, I>(
        &'a self,
        space: &'a HilbertSpace<I>,
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        self.nn_pairs_with(sites_on(space))
    }

    /// Reads top-level sites and delegates to the positional leaf-pair connection.
    pub fn nn_pairs_by_pos_on<'a, I>(
        &'a self,
        space: &'a HilbertSpace<I>,
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        self.nn_pairs_by_pos_with(sites_on(space))
    }
}

fn sites_on<I>(space: &HilbertSpace<I>) -> &[HilbertSpace<I>] {
    if space.factor_index().is_some() {
        core::slice::from_ref(space)
    } else {
        space.children().as_slice()
    }
}
