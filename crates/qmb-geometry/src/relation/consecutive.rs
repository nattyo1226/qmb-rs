use crate::{Boundary, Chain, Relation};

/// Invalid settings for a consecutive-site rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum ConsecutiveError {
    /// A consecutive placement must contain at least one site.
    #[error("consecutive length must be positive")]
    ZeroLength,
}

/// Ordered consecutive placements of exactly `K` sites on a chain.
///
/// Open chains use starts `0..=N-K`; periodic chains use all `0..N` starts
/// and wrap in site order. `K > N` produces no groups, and `K == 0` is invalid.
/// Rotations of the same support are retained, including `[0,1]` and `[1,0]`
/// for two periodic sites. No site is repeated within a placement.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Consecutive<const K: usize>;

impl<const K: usize> Relation<Chain> for Consecutive<K> {
    type Error = ConsecutiveError;

    type Group<'a> = [usize; K];

    fn groups<'a>(
        &'a self,
        geometry: &'a Chain,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
        if K == 0 {
            return Err(ConsecutiveError::ZeroLength);
        }
        let n = geometry.num_sites();
        let count = if K > n {
            0
        } else if geometry.boundary() == Boundary::Periodic {
            n
        } else {
            n - K + 1
        };
        Ok((0..count).map(move |start| placement_at::<K>(n, start)))
    }
}

// Called only for 1 <= K <= n and start < n.
fn placement_at<const K: usize>(n: usize, start: usize) -> [usize; K] {
    core::array::from_fn(|offset| {
        // Avoid start + offset overflow, including at the largest site counts.
        if offset >= n - start {
            offset - (n - start)
        } else {
            start + offset
        }
    })
}

#[cfg(test)]
mod tests {
    use super::placement_at;
    #[test]
    fn wrapped_placements_at_usize_max_do_not_add_past_the_limit() {
        let n = usize::MAX;
        assert_eq!(placement_at::<3>(n, n - 1), [n - 1, 0, 1]);
        assert_eq!(placement_at::<4>(n, n - 2), [n - 2, n - 1, 0, 1]);
        assert_eq!(placement_at::<4>(n, n - 4), [n - 4, n - 3, n - 2, n - 1]);
        assert_eq!(placement_at::<1>(n, n - 1), [n - 1]);
    }
}
