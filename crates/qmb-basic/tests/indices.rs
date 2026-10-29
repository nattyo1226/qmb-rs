use core::num::NonZeroUsize;
use qmb_basic::{HilbertSpace, Id128};

#[test]
fn uuid_generation_and_core_connection_preserve_existing_identity() {
    let a = Id128::new();
    let b = Id128::new();
    assert_ne!(a, b);
    assert_ne!(a, Id128::default());
    assert_ne!(Id128::default(), Id128::default());
    assert_eq!(core::mem::size_of::<Id128>(), 16);
    fn cloned<T: Clone>(value: &T) -> T {
        value.clone()
    }
    assert_eq!(cloned(&a), a);
    let copied = a;
    assert_eq!(copied, a);
    let mut set = std::collections::HashSet::new();
    set.insert(a);
    assert!(set.contains(&copied));
    let display = a.to_string();
    assert_eq!(display.len(), 36);
    assert_eq!(&display[14..15], "4");
    assert!(matches!(&display[19..20], "8" | "9" | "a" | "b"));
    let factor = HilbertSpace::factor(a);
    let copy = factor.clone();
    assert_eq!(copy.factor_index(), Some(&a));
    assert!(copy.is_compatible_with(&factor));
    let space = HilbertSpace::tensor_all([factor, HilbertSpace::factor(b)]).unwrap();
    assert_eq!(space.indices().copied().collect::<Vec<_>>(), vec![a, b]);
    assert_eq!(space.dim_with(NonZeroUsize::new(1).unwrap()), Ok(1));
}
