//! Basis expansion and sparse coordinates share aggregation, not result trees.

use alloc::collections::BTreeMap;
use num_traits::Zero;

mod composite;
mod local;

// Remove exact zeros before intermediate maps feed later tensor products.
fn accumulate<K: Ord, V, C: Zero>(map: &mut BTreeMap<K, (C, V)>, key: K, c: C, value: V) {
    use alloc::collections::btree_map::Entry;
    match map.entry(key) {
        Entry::Vacant(entry) => {
            if !c.is_zero() {
                entry.insert((c, value));
            }
        }
        Entry::Occupied(mut entry) => {
            let existing = &mut entry.get_mut().0;
            *existing = core::mem::replace(existing, C::zero()) + c;
            if existing.is_zero() {
                entry.remove();
            }
        }
    }
}
