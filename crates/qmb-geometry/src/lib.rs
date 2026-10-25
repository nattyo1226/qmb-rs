//! Immutable site geometries, independent of any concrete local algebra or index generator.
//!
//! Sites are positions in a vector. Grouping and leaf identity are carried by
//! existing Hilbert spaces; cloning a geometry shares their nodes.

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

extern crate alloc;

use alloc::vec::Vec;
use qmb_core::{HilbertSpace, LocalStructure, SpaceError};

/// Invalid site assignments or spaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum GeometryError {
    /// A site number is outside the geometry.
    #[error("site number is out of bounds")]
    SiteOutOfBounds,
    /// The supplied number of site spaces differs from the requested count.
    #[error("site count mismatch")]
    SiteCountMismatch,
    /// The sites do not partition the supplied full space.
    #[error("sites do not partition the full space")]
    IncompatibleSpace,
    /// A factor index occurs in more than one site.
    #[error("{0}")]
    Space(#[from] SpaceError),
}

/// Immutable site order and its full Hilbert space.
pub struct Geometry<I, L: LocalStructure> {
    sites: Vec<HilbertSpace<I, L>>,
    space: HilbertSpace<I, L>,
}

impl<I, L: LocalStructure> Clone for Geometry<I, L> {
    fn clone(&self) -> Self {
        Self {
            sites: self.sites.clone(),
            space: self.space.clone(),
        }
    }
}

impl<I, L: LocalStructure> Geometry<I, L> {
    /// Site-ordered factor or tensor spaces, including empty sites.
    pub fn sites(&self) -> &[HilbertSpace<I, L>] {
        &self.sites
    }

    /// Borrows the stored full space without issuing indices or rebuilding it.
    pub fn space(&self) -> &HilbertSpace<I, L> {
        &self.space
    }

    /// Looks up a site with an explicit out-of-bounds error.
    pub fn site(&self, site: usize) -> Result<&HilbertSpace<I, L>, GeometryError> {
        self.sites.get(site).ok_or(GeometryError::SiteOutOfBounds)
    }
}

impl<I: Eq, L: LocalStructure> Geometry<I, L> {
    /// Tensors existing site spaces in site order, preserving their grouping.
    pub fn from_sites(
        sites: impl IntoIterator<Item = HilbertSpace<I, L>>,
    ) -> Result<Self, GeometryError> {
        let sites: Vec<_> = sites.into_iter().collect();
        let space = HilbertSpace::tensor(sites.iter().cloned())?;
        Ok(Self { sites, space })
    }

    /// Uses a supplied full space without changing its factor order or grouping.
    ///
    /// Sites must partition its indices, and each site's internal leaf order
    /// must be a subsequence of the full order. Site order itself is independent.
    pub fn with_space(
        sites: impl IntoIterator<Item = HilbertSpace<I, L>>,
        space: &HilbertSpace<I, L>,
    ) -> Result<Self, GeometryError> {
        let mut geometry = Self::from_sites(sites)?;
        if geometry.space.indices().count() != space.indices().count()
            || !geometry
                .sites
                .iter()
                .all(|site| site.is_subsystem_of(space))
        {
            return Err(GeometryError::IncompatibleSpace);
        }
        geometry.space = space.clone();
        Ok(geometry)
    }
}

impl<I: Default + Eq, L: LocalStructure> Geometry<I, L> {
    /// Issues one index per local degree, in site order and then leaf order.
    ///
    /// `I::default()` must generate a fresh independent index on each call.
    /// One degree produces a factor, multiple degrees an ordered tensor.
    /// Zero degrees retain the requested number of empty tensor sites.
    pub fn new(site_count: usize, degrees_per_site: usize) -> Result<Self, GeometryError> {
        let sites = (0..site_count)
            .map(|_| match degrees_per_site {
                1 => Ok(HilbertSpace::factor(I::default())),
                _ => HilbertSpace::tensor(
                    (0..degrees_per_site).map(|_| HilbertSpace::factor(I::default())),
                ),
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_sites(sites)
    }
}

/// Boundary condition for an undirected chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// Connect consecutive sites only.
    Open,
    /// Also connect the last site to the first, without loops or duplicate edges.
    Periodic,
}

/// A chain with an immutable site assignment.
pub struct Chain<I, L: LocalStructure> {
    geometry: Geometry<I, L>,
    boundary: Boundary,
}

impl<I, L: LocalStructure> Clone for Chain<I, L> {
    fn clone(&self) -> Self {
        Self {
            geometry: self.geometry.clone(),
            boundary: self.boundary,
        }
    }
}

impl<I, L: LocalStructure> Chain<I, L> {
    /// Assigns chain adjacency to an existing geometry.
    pub fn from_geometry(geometry: Geometry<I, L>, boundary: Boundary) -> Self {
        Self { geometry, boundary }
    }
    /// Borrows the site assignment.
    pub fn geometry(&self) -> &Geometry<I, L> {
        &self.geometry
    }
    /// Borrows site spaces in site order.
    pub fn sites(&self) -> &[HilbertSpace<I, L>] {
        self.geometry.sites()
    }
    /// Borrows the stored full space.
    pub fn space(&self) -> &HilbertSpace<I, L> {
        self.geometry.space()
    }
    /// Returns the boundary condition.
    pub fn boundary(&self) -> Boundary {
        self.boundary
    }
    /// Enumerates normalized, unique site pairs in lexicographic order.
    pub fn nn_site_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        let n = self.sites().len();
        let closing = n > 2 && self.boundary == Boundary::Periodic;
        (0..n.saturating_sub(1)).flat_map(move |i| {
            core::iter::once((i, i + 1)).chain((closing && i == 0).then_some((0, n - 1)))
        })
    }
}

impl<I: Eq, L: LocalStructure> Chain<I, L> {
    /// Constructs a chain from existing site spaces.
    pub fn from_sites(
        sites: impl IntoIterator<Item = HilbertSpace<I, L>>,
        boundary: Boundary,
    ) -> Result<Self, GeometryError> {
        Ok(Self::from_geometry(Geometry::from_sites(sites)?, boundary))
    }
    /// Constructs a chain, rejecting a mismatch with the declared site count.
    pub fn with_sites(
        site_count: usize,
        sites: impl IntoIterator<Item = HilbertSpace<I, L>>,
        boundary: Boundary,
    ) -> Result<Self, GeometryError> {
        let geometry = Geometry::from_sites(sites)?;
        if geometry.sites().len() != site_count {
            return Err(GeometryError::SiteCountMismatch);
        }
        Ok(Self::from_geometry(geometry, boundary))
    }
}

impl<I: Default + Eq, L: LocalStructure> Chain<I, L> {
    /// Constructs uniform sites; `I::default()` must issue fresh indices.
    pub fn new(
        site_count: usize,
        degrees_per_site: usize,
        boundary: Boundary,
    ) -> Result<Self, GeometryError> {
        Ok(Self::from_geometry(
            Geometry::new(site_count, degrees_per_site)?,
            boundary,
        ))
    }
}

/// A geometry capable of enumerating nearest-neighbor index pairs.
pub trait NearestNeighborRelation {
    /// Identity type of a basic factor.
    type Index;
    /// All combinations across each neighboring site pair, in site-pair order.
    fn nn_pairs(&self) -> impl Iterator<Item = (&Self::Index, &Self::Index)> + '_;
    /// Matching leaf positions across each site pair, stopping at the shorter side.
    fn nn_pairs_by_pos(&self) -> impl Iterator<Item = (&Self::Index, &Self::Index)> + '_;
}

/// Cartesian leaf-index pairs, left outer and right inner, borrowing both spaces.
pub fn index_pairs<'a, I, L: LocalStructure>(
    left: &'a HilbertSpace<I, L>,
    right: &'a HilbertSpace<I, L>,
) -> impl Iterator<Item = (&'a I, &'a I)> + 'a {
    left.indices()
        .flat_map(move |a| right.indices().map(move |b| (a, b)))
}

/// Zips leaf-index positions, stopping at the shorter side.
pub fn index_pairs_by_pos<'a, I, L: LocalStructure>(
    left: &'a HilbertSpace<I, L>,
    right: &'a HilbertSpace<I, L>,
) -> impl Iterator<Item = (&'a I, &'a I)> + 'a {
    left.indices().zip(right.indices())
}

impl<I, L: LocalStructure> NearestNeighborRelation for Chain<I, L> {
    type Index = I;
    fn nn_pairs(&self) -> impl Iterator<Item = (&I, &I)> + '_ {
        self.nn_site_pairs()
            .flat_map(|(i, j)| index_pairs(&self.sites()[i], &self.sites()[j]))
    }
    fn nn_pairs_by_pos(&self) -> impl Iterator<Item = (&I, &I)> + '_ {
        self.nn_site_pairs()
            .flat_map(|(i, j)| index_pairs_by_pos(&self.sites()[i], &self.sites()[j]))
    }
}
