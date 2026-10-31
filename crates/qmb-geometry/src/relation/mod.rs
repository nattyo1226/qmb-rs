//! Site relations for each supported geometry.

mod consecutive;
mod nearest_neighbor;
use crate::Geometry;
pub use consecutive::{Consecutive, ConsecutiveError};
pub use nearest_neighbor::NearestNeighbors;

/// A geometric rule selecting ordered groups of site numbers on `G`.
///
/// Multiple rules may apply to one geometry, and one rule may support multiple
/// geometries. Each implementation defines ordering, repetitions and boundaries.
/// Group types can retain metadata in addition to their site slice.
///
/// Groups may borrow stored data from the rule or geometry. Generated groups
/// must be owned values: an ordinary iterator cannot lend a reusable buffer.
pub trait Relation<G: Geometry> {
    /// Invalid rule settings or other errors specific to this implementation.
    type Error;

    /// A group whose site slice remains stable while the group is borrowed.
    type Group<'a>: AsRef<[usize]>
    where
        Self: 'a,
        G: 'a;

    /// Validates settings and lazily enumerates groups on the supplied geometry.
    fn groups<'a>(
        &'a self,
        geometry: &'a G,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error>;
}
