use qmb_core::{HilbertSpace, SpaceError};
use qmb_geometry::{Boundary, Chain, ToSpace};

// Clone and Eq suffice; no Default, Copy, Hash or Ord is needed for conversion.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Id(usize);

#[test]
fn conversion_preserves_site_grouping_identity_order_and_empty_sites() {
    let chain = Chain::from_sites(
        [vec![Id(9), Id(3)], vec![], vec![Id(8)]],
        Boundary::Periodic,
    )
    .unwrap();
    let space = chain.to_space().unwrap();
    let again = chain.to_space().unwrap();
    assert_eq!(space.indices().map(|i| i.0).collect::<Vec<_>>(), [9, 3, 8]);
    assert!(space.is_compatible_with(&again));
    let sites: Vec<_> = space.children().collect();
    assert_eq!(sites.len(), 3);
    assert_eq!(sites[0].children().count(), 2);
    assert_eq!(sites[0].indices().map(|i| i.0).collect::<Vec<_>>(), [9, 3]);
    assert!(sites[1].indices().next().is_none());
    assert_eq!(sites[2].factor_index(), Some(&Id(8)));
    assert_eq!(chain.site(0).unwrap(), [Id(9), Id(3)]);
    assert!(chain.site(1).unwrap().is_empty());
    assert_eq!(chain.site(2).unwrap(), [Id(8)]);
    assert_eq!(chain.boundary(), Boundary::Periodic);
}

#[test]
fn empty_chain_and_empty_sites_give_unit_spaces_without_erasing_sites() {
    let empty = Chain::<Id>::from_sites(Vec::<Vec<Id>>::new(), Boundary::Open).unwrap();
    let space = empty.to_space().unwrap();
    assert!(space.is_compatible_with(&HilbertSpace::unit()));
    assert_eq!(space.children().count(), 0);
    let chain = Chain::<Id>::from_sites([vec![], vec![], vec![]], Boundary::Periodic).unwrap();
    let space = chain.to_space().unwrap();
    assert!(space.is_compatible_with(&HilbertSpace::unit()));
    assert_eq!(space.children().count(), 3);
    assert_eq!(chain.sites().len(), 3);
}

#[test]
fn converted_space_outlives_the_chain_and_keeps_borrowed_indices() {
    let names = [String::from("a"), String::from("b")];
    let space = {
        let chain = Chain::from_sites(names.iter().map(|s| [s.as_str()]), Boundary::Open).unwrap();
        chain.to_space().unwrap()
    };
    assert_eq!(space.indices().copied().collect::<Vec<_>>(), ["a", "b"]);
}

#[test]
fn trait_itself_requires_neither_index_clone_nor_equality() {
    struct Opaque;
    struct Existing(HilbertSpace<Opaque>);
    impl ToSpace for Existing {
        type Index = Opaque;
        fn to_space(&self) -> Result<HilbertSpace<Opaque>, SpaceError> {
            Ok(self.0.clone())
        }
    }
    let source = Existing(HilbertSpace::factor(Opaque));
    let space = source.to_space().unwrap();
    assert!(space.shares_node_with(&source.0));
}
