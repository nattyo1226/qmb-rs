//! Geometries and their structural information.

pub(crate) mod chain;

/// Common indexing contract for a geometry with a caller-owned site sequence.
///
/// Site numbers are `0..num_sites()`. Their meaning and count must remain stable
/// while a borrowing mapper is in use. Coordinates, directions and adjacency
/// belong to concrete geometries, rather than this minimal common interface.
pub trait Geometry {
    /// Returns the number of sites, including any assigned empty sites.
    fn num_sites(&self) -> usize;
}
