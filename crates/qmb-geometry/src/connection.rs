//! Connect site relations to caller-owned sequences and Hilbert spaces.

use qmb_core::HilbertSpace;

/// Failure to connect a chain to a site-ordered input.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum ChainError {
    /// The supplied number of sites differs from the chain's count.
    #[error("site count mismatch")]
    SiteCountMismatch,
}

pub(crate) fn site_pairs_with<'a, T>(
    site_count: usize,
    pairs: impl Iterator<Item = (usize, usize)> + 'a,
    values: &'a [T],
) -> Result<impl Iterator<Item = (&'a T, &'a T)> + 'a, ChainError> {
    if values.len() != site_count {
        return Err(ChainError::SiteCountMismatch);
    }
    Ok(pairs.map(|(i, j)| (&values[i], &values[j])))
}

pub(crate) fn leaf_pairs<'a, I: 'a>(
    pairs: impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a,
) -> impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a {
    pairs.flat_map(|(left, right)| {
        left.leaves()
            .flat_map(move |a| right.leaves().map(move |b| (a, b)))
    })
}

pub(crate) fn leaf_pairs_by_pos<'a, I: 'a>(
    pairs: impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a,
) -> impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a {
    pairs.flat_map(|(left, right)| left.leaves().zip(right.leaves()))
}

pub(crate) fn sites_on<I>(space: &HilbertSpace<I>) -> &[HilbertSpace<I>] {
    if space.factor_index().is_some() {
        core::slice::from_ref(space)
    } else {
        space.children().as_slice()
    }
}
