//! Standard UUID indices. Randomness and platform dependencies stay outside core.
//!
//! ```
//! use qmb_basic::Index;
//! let index = Index::new();
//! let same = index;
//! assert_eq!(index, same);
//! assert_ne!(index, Index::new());
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt;

/// An opaque, randomly generated version 4 UUID identifying one basic factor.
///
/// Copies preserve identity. Only explicit generation issues a new index.
/// The UUID library obtains randomness from the operating system; callers do
/// not supply a seed. Generation follows `uuid::Uuid::new_v4`'s platform and
/// failure behavior (it can panic if randomness is unavailable).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Index(uuid::Uuid);

impl Index {
    /// Issues a fresh UUID without a local-structure argument or seed.
    pub fn new() -> Self {
        Self(uuid::Uuid::new_v4())
    }
}

impl Default for Index {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Index {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// A core space using the standard UUID index.
pub type HilbertSpace<L> = qmb_core::HilbertSpace<Index, L>;
