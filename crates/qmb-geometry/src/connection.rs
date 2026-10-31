//! Borrowing connections between geometric groups and caller-owned data.
use crate::{Geometry, Relation};
use alloc::{vec, vec::Vec};
use qmb_core::{HilbertSpace, space::Leaves};

/// Failure to connect site numbers to a caller-owned sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum MappingError {
    /// The supplied number of sites differs from the geometry's count.
    #[error("site count mismatch: expected {expected}, got {actual}")]
    SiteCountMismatch {
        /// Geometry's site count.
        expected: usize,
        /// Input slice length.
        actual: usize,
    },
    /// A custom relation returned a site outside the mapped sequence.
    #[error("site {site} is outside 0..{num_sites}")]
    SiteOutOfRange {
        /// Invalid site number.
        site: usize,
        /// Number of mapped sites.
        num_sites: usize,
    },
}

/// A geometry connected to a caller-owned site-ordered sequence.
///
/// Checks exact length even when a relation would produce no groups. Coordinate
/// meanings remain the caller's responsibility. No bounds are required on `T`.
pub struct SiteMapper<'g, 's, G, T> {
    geometry: &'g G,
    values: &'s [T],
}

impl<'g, 's, G: Geometry, T> SiteMapper<'g, 's, G, T> {
    /// Checks the site count and borrows both inputs.
    pub fn new(geometry: &'g G, values: &'s [T]) -> Result<Self, MappingError> {
        if geometry.num_sites() != values.len() {
            return Err(MappingError::SiteCountMismatch {
                expected: geometry.num_sites(),
                actual: values.len(),
            });
        }
        Ok(Self { geometry, values })
    }

    /// Returns the original geometry.
    pub fn geometry(&self) -> &'g G {
        self.geometry
    }

    /// Returns the original site-ordered slice.
    pub fn values(&self) -> &'s [T] {
        self.values
    }

    /// Applies a rule, returning its own error or a lazy iterator of mapping results.
    ///
    /// Each group is range-checked before mapping; previous results are not rolled back.
    pub fn groups<'a, R: Relation<G>>(
        &'a self,
        relation: &'a R,
    ) -> Result<impl Iterator<Item = MappedGroupResult<'s, T, R::Group<'a>>> + 'a, R::Error> {
        let groups = relation.groups(self.geometry)?;
        let values = self.values;
        Ok(groups.map(move |source| {
            if let Some(&site) = source.as_ref().iter().find(|&&site| site >= values.len()) {
                return Err(MappingError::SiteOutOfRange {
                    site,
                    num_sites: values.len(),
                });
            }
            Ok(MappedGroup { values, source })
        }))
    }
}

/// An original typed group connected to borrowed site values in group order.
pub struct MappedGroup<'s, T, Group> {
    values: &'s [T],
    source: Group,
}

/// Per-group result returned by [`SiteMapper::groups`].
pub type MappedGroupResult<'s, T, Group> = Result<MappedGroup<'s, T, Group>, MappingError>;

impl<'s, T, Group: AsRef<[usize]>> MappedGroup<'s, T, Group> {
    /// Returns the mapped value at a group position, rather than a site number.
    pub fn get(&self, position: usize) -> Option<&'s T> {
        self.source
            .as_ref()
            .get(position)
            .map(|&site| &self.values[site])
    }

    /// Returns the number of sites in this group.
    pub fn len(&self) -> usize {
        self.source.as_ref().len()
    }

    /// Returns whether this group contains no sites.
    pub fn is_empty(&self) -> bool {
        self.source.as_ref().is_empty()
    }

    /// Returns the original group, retaining relation-specific metadata.
    pub fn source(&self) -> &Group {
        &self.source
    }

    /// Borrows values in group order, including repetitions.
    pub fn iter(&self) -> impl Iterator<Item = &'s T> + '_ {
        let values = self.values;
        self.source.as_ref().iter().map(move |&site| &values[site])
    }
}

impl<'s, I, Group: AsRef<[usize]>> MappedGroup<'s, HilbertSpace<I>, Group> {
    /// Calls a custom combiner once with each site's leaves in tensor order.
    ///
    /// Output items are unrestricted. Neither Index nor Leaves needs Clone or Eq.
    pub fn combine_leaves<Out: Iterator>(
        &self,
        combine: impl FnOnce(Vec<Leaves<'s, I>>) -> Out,
    ) -> Out {
        combine(self.iter().map(HilbertSpace::leaves).collect())
    }

    /// Selects one leaf per site in all combinations, first site slowest.
    ///
    /// Caches leaf references for replay. An empty site gives no rows; a zero-site
    /// group gives one empty row. Results are generated lazily.
    pub fn leaves_cartesian(&self) -> impl Iterator<Item = Vec<&'s HilbertSpace<I>>> {
        self.combine_leaves(cartesian_leaves)
    }

    /// Zips leaf positions up to the shortest site in site order.
    ///
    /// An empty site gives no rows; a zero-site group gives one empty row.
    pub fn leaves_zip(&self) -> impl Iterator<Item = Vec<&'s HilbertSpace<I>>> {
        self.combine_leaves(zip_leaves)
    }
}

// Cache borrowed leaves for replay: Leaves itself need not implement Clone.
// No Eq/Clone bound on Index; the cached items are shared references.
fn cartesian_leaves<I>(
    mut streams: Vec<Leaves<'_, I>>,
) -> impl Iterator<Item = Vec<&HilbertSpace<I>>> {
    let mut cache: Vec<Vec<_>> = (0..streams.len()).map(|_| Vec::new()).collect();
    let mut positions = vec![0; streams.len()];
    let mut started = false;
    let mut done = false;
    core::iter::from_fn(move || {
        if done {
            return None;
        }
        if !started {
            for (stream, cached) in streams.iter_mut().zip(&mut cache) {
                let Some(first) = stream.next() else {
                    done = true;
                    return None;
                };
                cached.push(first);
            }
            started = true;
        } else {
            let mut advanced = false;
            for slot in (0..streams.len()).rev() {
                let next = positions[slot] + 1;
                if next < cache[slot].len() {
                    positions[slot] = next;
                    advanced = true;
                    break;
                }
                if let Some(leaf) = streams[slot].next() {
                    cache[slot].push(leaf);
                    positions[slot] = next;
                    advanced = true;
                    break;
                }
                positions[slot] = 0;
            }
            if !advanced {
                done = true;
                return None;
            }
        }
        Some(
            cache
                .iter()
                .zip(&positions)
                .map(|(leaves, &pos)| leaves[pos])
                .collect(),
        )
    })
}

fn zip_leaves<I>(mut streams: Vec<Leaves<'_, I>>) -> impl Iterator<Item = Vec<&HilbertSpace<I>>> {
    let mut done = false;
    core::iter::from_fn(move || {
        if done {
            return None;
        }
        if streams.is_empty() {
            done = true;
            return Some(Vec::new());
        }
        let row: Option<Vec<_>> = streams.iter_mut().map(Iterator::next).collect();
        if row.is_none() {
            done = true;
        }
        row
    })
}

impl<'g, 's, G: Geometry, I> SiteMapper<'g, 's, G, HilbertSpace<I>> {
    /// Reads a space's stored top-level partition as sites.
    ///
    /// A leaf is one site; the unit tensor is zero sites. A tensor containing
    /// one unit is one empty site. Nested site grouping is preserved.
    pub fn from_space(geometry: &'g G, space: &'s HilbertSpace<I>) -> Result<Self, MappingError> {
        let sites = if space.factor_index().is_some() {
            core::slice::from_ref(space)
        } else {
            space.children().as_slice()
        };
        Self::new(geometry, sites)
    }
}
