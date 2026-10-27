//! Immutable site index assignments and geometric relations.
//!
//! A chain owns ordered index lists and its boundary condition. Space conversion
//! is available separately; callers manage Hilbert spaces, grouping and local algebras.

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

extern crate alloc;

use alloc::vec::Vec;

pub mod space;
pub use space::ToSpace;

/// Invalid chain assignments or site references.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum ChainError {
    /// A site number is outside the chain.
    #[error("site number is out of bounds")]
    SiteOutOfBounds,
    /// The supplied number of sites differs from the requested count.
    #[error("site count mismatch")]
    SiteCountMismatch,
    /// An index occurs more than once, within or across sites.
    #[error("duplicate factor index")]
    DuplicateIndex,
}

/// Boundary condition for an undirected chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// Connect consecutive sites only.
    Open,
    /// Also connect the last site to the first, without loops or duplicate edges.
    Periodic,
}

/// An immutable chain owning ordered indices for each site.
///
/// Empty and variable-length sites are retained. Construction checks index
/// uniqueness with equality alone; borrowing sites or relations requires no
/// index cloning or generation. Cloning a chain requires `I: Clone` and copies
/// the same identities without issuing new indices.
///
/// ```
/// use qmb_geometry::{Boundary, Chain, NearestNeighborRelation};
/// let chain = Chain::from_sites([vec![9, 3], vec![8], vec![]], Boundary::Open)?;
/// assert_eq!(chain.site(0)?, [9, 3]);
/// assert_eq!(chain.nn_pairs().map(|(&a, &b)| (a, b)).collect::<Vec<_>>(),
///            [(9, 8), (3, 8)]);
/// # Ok::<(), qmb_geometry::ChainError>(())
/// ```
#[derive(Clone)]
pub struct Chain<I> {
    sites: Vec<Vec<I>>,
    boundary: Boundary,
}

impl<I> Chain<I> {
    /// Borrows all index lists in site order.
    pub fn sites(&self) -> &[Vec<I>] {
        &self.sites
    }

    /// Borrows a site's ordered indices, with an explicit bounds error.
    pub fn site(&self, site: usize) -> Result<&[I], ChainError> {
        self.sites
            .get(site)
            .map(Vec::as_slice)
            .ok_or(ChainError::SiteOutOfBounds)
    }

    /// Returns the boundary condition.
    pub fn boundary(&self) -> Boundary {
        self.boundary
    }

    /// Enumerates normalized, unique site pairs in lexicographic order.
    pub fn nn_site_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        let n = self.sites.len();
        let closing = n > 2 && self.boundary == Boundary::Periodic;
        (0..n.saturating_sub(1)).flat_map(move |i| {
            core::iter::once((i, i + 1)).chain((closing && i == 0).then_some((0, n - 1)))
        })
    }
}

impl<I: Eq> Chain<I> {
    /// Moves existing indices into a chain, rejecting duplicates within or across sites.
    pub fn from_sites<S: IntoIterator<Item = I>>(
        sites: impl IntoIterator<Item = S>,
        boundary: Boundary,
    ) -> Result<Self, ChainError> {
        let sites: Vec<Vec<I>> = sites
            .into_iter()
            .map(|site| site.into_iter().collect())
            .collect();
        let mut seen = Vec::new();
        for index in sites.iter().flatten() {
            if seen.contains(&index) {
                return Err(ChainError::DuplicateIndex);
            }
            seen.push(index);
        }
        Ok(Self { sites, boundary })
    }

    /// Constructs a chain and checks the explicitly declared site count.
    pub fn with_sites<S: IntoIterator<Item = I>>(
        site_count: usize,
        sites: impl IntoIterator<Item = S>,
        boundary: Boundary,
    ) -> Result<Self, ChainError> {
        let chain = Self::from_sites(sites, boundary)?;
        if chain.sites.len() != site_count {
            return Err(ChainError::SiteCountMismatch);
        }
        Ok(chain)
    }
}

impl<I: Default + Eq> Chain<I> {
    /// Issues indices in site order and then within each site.
    ///
    /// `I::default()` must issue a fresh independent identity each time. Zero
    /// degrees retain the requested number of empty sites; zero sites issue none.
    pub fn new(
        site_count: usize,
        degrees_per_site: usize,
        boundary: Boundary,
    ) -> Result<Self, ChainError> {
        Self::from_sites(
            (0..site_count).map(|_| (0..degrees_per_site).map(|_| I::default())),
            boundary,
        )
    }
}

/// A geometry capable of enumerating nearest-neighbor index pairs.
pub trait NearestNeighborRelation {
    /// Identity type of a basic factor.
    type Index;
    /// All combinations across each neighboring site pair, in site-pair order.
    fn nn_pairs(&self) -> impl Iterator<Item = (&Self::Index, &Self::Index)> + '_;
    /// Matching index positions across each site pair, stopping at the shorter side.
    fn nn_pairs_by_pos(&self) -> impl Iterator<Item = (&Self::Index, &Self::Index)> + '_;
}

/// Cartesian index pairs, left outer and right inner, borrowing both lists.
pub fn index_pairs<'a, I>(
    left: &'a [I],
    right: &'a [I],
) -> impl Iterator<Item = (&'a I, &'a I)> + 'a {
    left.iter()
        .flat_map(move |a| right.iter().map(move |b| (a, b)))
}

/// Zips index positions, stopping at the shorter list.
pub fn index_pairs_by_pos<'a, I>(
    left: &'a [I],
    right: &'a [I],
) -> impl Iterator<Item = (&'a I, &'a I)> + 'a {
    left.iter().zip(right.iter())
}

impl<I> NearestNeighborRelation for Chain<I> {
    type Index = I;
    fn nn_pairs(&self) -> impl Iterator<Item = (&I, &I)> + '_ {
        self.nn_site_pairs()
            .flat_map(|(i, j)| index_pairs(&self.sites[i], &self.sites[j]))
    }
    fn nn_pairs_by_pos(&self) -> impl Iterator<Item = (&I, &I)> + '_ {
        self.nn_site_pairs()
            .flat_map(|(i, j)| index_pairs_by_pos(&self.sites[i], &self.sites[j]))
    }
}
