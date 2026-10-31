use crate::{Boundary, Chain, Relation};
use core::convert::Infallible;

/// Undirected nearest-neighbor edges, normalized and unique.
///
/// On a chain, returns `[i, j]` with `i < j` in lexicographic order. Zero and
/// one sites have no edges; two sites have one for either boundary condition.
/// This differs from ordered [`crate::Consecutive<2>`] on periodic chains.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NearestNeighbors;

impl Relation<Chain> for NearestNeighbors {
    type Error = Infallible;
    type Group<'a> = [usize; 2];
    fn groups<'a>(
        &'a self,
        geometry: &'a Chain,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
        let n = geometry.num_sites();
        let closing = n > 2 && geometry.boundary() == Boundary::Periodic;
        Ok((0..n.saturating_sub(1)).flat_map(move |i| {
            core::iter::once([i, i + 1]).chain((closing && i == 0).then(|| [0, n - 1]))
        }))
    }
}
