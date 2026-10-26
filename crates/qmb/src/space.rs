//! Connections from geometric index assignments to ordinary Hilbert spaces.

use qmb_core::{HilbertSpace, SpaceError};
use qmb_geometry::Chain;

/// Converts a borrowed object to an ordinary space without issuing indices.
///
/// Implementations preserve index identities and document their factor order.
/// The returned space need not share nodes with previous conversions. Index
/// cloning is a requirement of individual implementations, not of this trait.
///
/// ```
/// use qmb::prelude::*;
/// use qmb::{Boundary, Chain};
/// let chain = Chain::from_sites([vec![9, 3], vec![], vec![8]], Boundary::Open)?;
/// let space = chain.to_space()?;
/// assert_eq!(space.indices().copied().collect::<Vec<_>>(), [9, 3, 8]);
/// assert_eq!(space.children().count(), 3);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait ToSpace {
    /// Index type preserved by the conversion.
    type Index;

    /// Builds or shares a space in the implementation's documented factor order.
    fn to_space(&self) -> Result<HilbertSpace<Self::Index>, SpaceError>;
}

impl<I: Clone + Eq> ToSpace for Chain<I> {
    type Index = I;

    /// Builds factors in site order and then index order within each site.
    ///
    /// Each site is a tensor except that a single index is represented directly
    /// by its factor. Empty sites are retained, and an empty chain gives the unit
    /// space. This conversion leaves the chain and its boundary unchanged.
    fn to_space(&self) -> Result<HilbertSpace<I>, SpaceError> {
        let sites = self.sites().iter().map(|indices| match indices.as_slice() {
            [index] => Ok(HilbertSpace::factor(index.clone())),
            _ => HilbertSpace::tensor(indices.iter().cloned().map(HilbertSpace::factor)),
        });
        HilbertSpace::tensor(sites.collect::<Result<Vec<_>, _>>()?)
    }
}
