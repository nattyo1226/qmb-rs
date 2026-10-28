//! Nearest-neighbor site relations for chains.

use crate::{Boundary, Chain};

pub(crate) fn site_pairs(chain: &Chain) -> impl Iterator<Item = (usize, usize)> + '_ {
    let n = chain.site_count();
    let closing = n > 2 && chain.boundary() == Boundary::Periodic;
    (0..n.saturating_sub(1)).flat_map(move |i| {
        core::iter::once((i, i + 1)).chain((closing && i == 0).then_some((0, n - 1)))
    })
}
