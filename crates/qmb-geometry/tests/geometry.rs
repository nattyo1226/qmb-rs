use qmb_geometry::{Boundary::*, Chain, ChainError};

#[test]
fn boundaries_are_normalized_sorted_unique_and_loop_free() {
    for n in 0..7 {
        for boundary in [Open, Periodic] {
            let chain = Chain::new(n, boundary);
            let edges: Vec<_> = chain.nn_site_pairs().collect();
            let mut expected: Vec<_> = (0..n.saturating_sub(1)).map(|i| (i, i + 1)).collect();
            if boundary == Periodic && n > 2 {
                expected.push((0, n - 1));
                expected.sort();
            }
            assert_eq!(edges, expected);
            assert!(edges.iter().all(|&(i, j)| i < j && j < n));
            assert!(edges.windows(2).all(|w| w[0] < w[1]));
            assert_eq!(chain.site_count(), n);
            assert_eq!(chain.boundary(), boundary);
            // No Clone, Eq, Default, Hash or Ord requirement on site values.
            struct Opaque(usize);
            let values: Vec<_> = (0..n).map(Opaque).collect();
            let mapped: Vec<_> = chain.nn_site_pairs_with(&values).unwrap().collect();
            assert_eq!(mapped.len(), expected.len());
            for ((a, b), (i, j)) in mapped.into_iter().zip(expected) {
                assert_eq!((a.0, b.0), (i, j));
                assert!(std::ptr::eq(a, &values[i]));
                assert!(std::ptr::eq(b, &values[j]));
            }
        }
    }
}

#[test]
fn exact_length_is_validated_before_returning_even_without_pairs() {
    for n in 0..=4 {
        let chain = Chain::new(n, Periodic);
        for count in 0..=5 {
            let values = vec![(); count];
            assert_eq!(
                chain.nn_site_pairs_with(&values).err(),
                (n != count).then_some(ChainError::SiteCountMismatch)
            );
        }
    }
    let error = ChainError::SiteCountMismatch;
    assert_eq!(error.to_string(), "site count mismatch");
    assert!(std::error::Error::source(&error).is_none());
}
