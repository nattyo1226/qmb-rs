use core::num::NonZeroUsize;
use qmb_basic::{HilbertSpace, Index};
use qmb_core::LocalStructure;

struct Unit;
impl LocalStructure for Unit {
    type State = ();
    type Operator = ();
    const DIMENSION: NonZeroUsize = NonZeroUsize::new(1).unwrap();
}

#[test]
fn uuid_generation_and_core_connection_preserve_existing_identity() {
    let a = Index::new();
    let b = Index::new();
    assert_ne!(a, b);
    let display = a.to_string();
    assert_eq!(display.len(), 36);
    assert_eq!(&display[14..15], "4");
    assert!(matches!(&display[19..20], "8" | "9" | "a" | "b"));
    let factor = HilbertSpace::<Unit>::factor(a);
    let copy = factor.clone();
    assert_eq!(copy.factor_index(), Some(&a));
    assert!(copy.is_compatible_with(&factor));
    let space = HilbertSpace::tensor([factor, HilbertSpace::factor(b)]).unwrap();
    assert_eq!(space.indices().copied().collect::<Vec<_>>(), vec![a, b]);
    assert_eq!(space.dimension(), Ok(1));
}
