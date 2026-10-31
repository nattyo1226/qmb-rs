use crate::{Boundary, Geometry};

/// Immutable site count and boundary condition, without indices or assignments.
///
/// Use [`crate::NearestNeighbors`] or [`crate::Consecutive`] to select sites,
/// and [`crate::SiteMapper`] to connect them to caller-owned data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chain {
    num_sites: usize,
    boundary: Boundary,
}

impl Chain {
    /// Constructs geometry without issuing indices; zero sites are allowed.
    pub fn new(num_sites: usize, boundary: Boundary) -> Self {
        Self {
            num_sites,
            boundary,
        }
    }

    /// Returns the number of sites, including any assigned empty sites.
    pub fn num_sites(&self) -> usize {
        self.num_sites
    }

    /// Returns the boundary condition.
    pub fn boundary(&self) -> Boundary {
        self.boundary
    }
}

impl Geometry for Chain {
    fn num_sites(&self) -> usize {
        self.num_sites
    }
}
