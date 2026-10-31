//! Bounded storage for variable-length site groups.

/// A site group exceeded its fixed storage capacity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum SiteGroupError {
    /// The requested group does not fit in its storage.
    #[error("site group length {len} exceeds capacity {capacity}")]
    CapacityExceeded {
        /// Maximum number of stored sites.
        capacity: usize,
        /// Requested number of sites.
        len: usize,
    },
}

/// A variable-length site group backed by an array of maximum length `CAP`.
///
/// Construction preserves order and repetitions. Each relation determines their
/// meaning; the group does not sort, deduplicate or validate a geometry's range.
/// No heap allocation is needed, and unused storage is never exposed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SiteGroup<const CAP: usize> {
    sites: [usize; CAP],
    len: usize,
}

impl<const CAP: usize> SiteGroup<CAP> {
    /// Copies the supplied sites, rejecting a length greater than `CAP`.
    pub fn new(sites: impl AsRef<[usize]>) -> Result<Self, SiteGroupError> {
        let sites = sites.as_ref();
        if sites.len() > CAP {
            return Err(SiteGroupError::CapacityExceeded {
                capacity: CAP,
                len: sites.len(),
            });
        }
        let mut result = Self::default();
        result.sites[..sites.len()].copy_from_slice(sites);
        result.len = sites.len();
        Ok(result)
    }

    /// Borrows the used portion of the group in its original order.
    pub fn as_slice(&self) -> &[usize] {
        &self.sites[..self.len]
    }

    /// Returns the number of used sites, rather than the storage capacity.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns whether the group contains no sites.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

impl<const CAP: usize> Default for SiteGroup<CAP> {
    fn default() -> Self {
        Self {
            sites: [0; CAP],
            len: 0,
        }
    }
}

impl<const CAP: usize> AsRef<[usize]> for SiteGroup<CAP> {
    fn as_ref(&self) -> &[usize] {
        self.as_slice()
    }
}
