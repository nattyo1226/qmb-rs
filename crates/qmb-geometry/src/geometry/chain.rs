use crate::{Boundary, ChainError, connection, relation::nearest_neighbor};
use qmb_core::HilbertSpace;

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
        nearest_neighbor::site_pairs(self)
    }

    /// Maps site pairs to borrowed elements, without cloning or inspecting them.
    /// Validates the exact input length before returning even an empty iterator.
    pub fn nn_site_pairs_with<'a, T>(
        &'a self,
        values: &'a [T],
    ) -> Result<impl Iterator<Item = (&'a T, &'a T)> + 'a, ChainError> {
        connection::site_pairs_with(self.site_count(), self.nn_site_pairs(), values)
    }

    /// Enumerates all local-factor combinations across neighboring sites, left outermost.
    /// Within each site, factors follow their stored tensor order.
    /// Empty sites retain their positions and contribute no leaf pairs.
    /// Index uniqueness is the responsibility of space or model construction.
    pub fn nn_pairs_with<'a, I>(
        &'a self,
        factors: &'a [HilbertSpace<I>],
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        Ok(connection::leaf_pairs(self.nn_site_pairs_with(factors)?))
    }

    /// Pairs corresponding local factors in stored tensor order.
    /// Stops at the shorter site's end.
    pub fn nn_zip_with<'a, I>(
        &'a self,
        factors: &'a [HilbertSpace<I>],
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        Ok(connection::leaf_pairs_by_pos(
            self.nn_site_pairs_with(factors)?,
        ))
    }

    /// Reads the stored top-level partition as sites and borrows site pairs.
    /// A leaf is one site; an empty tensor is zero sites. Nested grouping stays intact.
    pub fn nn_site_pairs_in<'a, I>(
        &'a self,
        space: &'a HilbertSpace<I>,
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        self.nn_site_pairs_with(connection::sites_on(space))
    }

    /// Reads top-level sites and delegates to the Cartesian leaf-pair connection.
    pub fn nn_pairs_in<'a, I>(
        &'a self,
        space: &'a HilbertSpace<I>,
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        self.nn_pairs_with(connection::sites_on(space))
    }

    /// Reads top-level sites and delegates to the positional leaf-pair connection.
    pub fn nn_zip_in<'a, I>(
        &'a self,
        space: &'a HilbertSpace<I>,
    ) -> Result<impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a, ChainError>
    {
        self.nn_zip_with(connection::sites_on(space))
    }
}
