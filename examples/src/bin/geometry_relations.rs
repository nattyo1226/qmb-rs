//! Geometry relations using the library API and a custom open square grid.
//! The example-local grid demonstrates custom Geometry and Relation implementations.

use qmb_core::HilbertSpace;
use qmb_geometry::{
    Boundary, Chain, Consecutive, Geometry, MappedGroup, NearestNeighbors, Relation, SiteMapper,
};
#[cfg(test)]
use qmb_geometry::{ConsecutiveError, MappingError};

type Error = Box<dyn std::error::Error>;

#[derive(Debug, PartialEq, Eq)]
enum GridError {
    InvalidLength,
    InvalidGridSize,
}
impl std::fmt::Display for GridError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for GridError {}

// Deliberately arbitrary groups for mapper tests, not a proposed geometric relation.
#[cfg(test)]
struct TestGroups {
    groups: Vec<Vec<usize>>,
}

#[cfg(test)]
impl<G: Geometry> Relation<G> for TestGroups {
    type Error = GridError;

    type Group<'a>
        = &'a [usize]
    where
        Self: 'a,
        G: 'a;

    fn groups<'a>(
        &'a self,
        _geometry: &'a G,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
        Ok(self.groups.iter().map(Vec::as_slice))
    }
}

#[derive(Clone, Copy)]
enum Axis {
    X,
    Y,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Direction {
    East,
    North,
    West,
    South,
}

struct OpenSquare {
    width: usize,
    height: usize,
    count: usize,
}

impl OpenSquare {
    fn new(width: usize, height: usize) -> Result<Self, GridError> {
        let count = width
            .checked_mul(height)
            .ok_or(GridError::InvalidGridSize)?;
        Ok(Self {
            width,
            height,
            count,
        })
    }

    fn step(&self, site: usize, direction: Direction) -> Option<usize> {
        let (x, y) = (site % self.width, site / self.width);
        match direction {
            Direction::East if x + 1 < self.width => Some(site + 1),
            Direction::North if y + 1 < self.height => Some(site + self.width),
            Direction::West if x > 0 => Some(site - 1),
            Direction::South if y > 0 => Some(site - self.width),
            _ => None,
        }
    }
}

impl Geometry for OpenSquare {
    fn num_sites(&self) -> usize {
        self.count
    }
}

impl Relation<OpenSquare> for NearestNeighbors {
    type Error = GridError;

    type Group<'a> = [usize; 2];

    fn groups<'a>(
        &'a self,
        geometry: &'a OpenSquare,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
        Ok((0..geometry.num_sites()).flat_map(move |site| {
            [Direction::East, Direction::North]
                .into_iter()
                .filter_map(move |direction| {
                    geometry.step(site, direction).map(|next| [site, next])
                })
        }))
    }
}

struct AlongAxis<const K: usize> {
    axis: Axis,
}

impl<const K: usize> Relation<OpenSquare> for AlongAxis<K> {
    type Error = GridError;

    type Group<'a> = [usize; K];

    fn groups<'a>(
        &'a self,
        geometry: &'a OpenSquare,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
        if K == 0 {
            return Err(GridError::InvalidLength);
        }
        let direction = match self.axis {
            Axis::X => Direction::East,
            Axis::Y => Direction::North,
        };
        Ok((0..geometry.num_sites()).filter_map(move |start| {
            let mut sites = [start; K];
            for slot in 1..K {
                sites[slot] = geometry.step(sites[slot - 1], direction)?;
            }
            Some(sites)
        }))
    }
}

struct Neighborhood {
    sites: [usize; 5],
    directions: [Option<Direction>; 5],
    len: usize,
}

impl AsRef<[usize]> for Neighborhood {
    fn as_ref(&self) -> &[usize] {
        &self.sites[..self.len]
    }
}

struct ClosedNeighborhood;

impl Relation<OpenSquare> for ClosedNeighborhood {
    type Error = GridError;

    type Group<'a> = Neighborhood;

    fn groups<'a>(
        &'a self,
        geometry: &'a OpenSquare,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
        Ok((0..geometry.num_sites()).map(|center| {
            let mut group = Neighborhood {
                sites: [center; 5],
                directions: [None; 5],
                len: 1,
            };
            for direction in [
                Direction::East,
                Direction::North,
                Direction::West,
                Direction::South,
            ] {
                if let Some(site) = geometry.step(center, direction) {
                    group.sites[group.len] = site;
                    group.directions[group.len] = Some(direction);
                    group.len += 1;
                }
            }
            group
        }))
    }
}

// A local extension trait supplies capabilities specific to this custom group.
trait NeighborhoodView<'s, T: 's> {
    fn center(&self) -> &'s T;
    fn neighbors(&self) -> impl Iterator<Item = (Direction, &'s T)> + '_;
}
impl<'s, T> NeighborhoodView<'s, T> for MappedGroup<'s, T, Neighborhood> {
    fn center(&self) -> &'s T {
        self.get(0).expect("center")
    }
    fn neighbors(&self) -> impl Iterator<Item = (Direction, &'s T)> + '_ {
        self.source().directions[1..self.source().len]
            .iter()
            .zip(self.iter().skip(1))
            .map(|(direction, value)| (direction.expect("neighbor direction"), value))
    }
}

fn main() -> Result<(), Error> {
    let chain = Chain::new(5, Boundary::Open);
    let chain_spaces: Vec<_> = (0..5).map(HilbertSpace::factor).collect();
    let mapped = SiteMapper::new(&chain, &chain_spaces)?;
    println!(
        "Chain: {} nearest-neighbor edges",
        mapped.groups(&NearestNeighbors)?.count()
    );
    let mut count = 0;
    for group in mapped.groups(&Consecutive::<3>)? {
        let group = group?;
        assert_eq!(group.iter().count(), 3);
        assert_eq!(group.source().len(), 3);
        count += 1;
    }
    println!("Chain: {count} consecutive triples");

    let square = OpenSquare::new(3, 3)?;
    let square_spaces: Vec<_> = (0..9).map(HilbertSpace::factor).collect();
    let mapped = SiteMapper::new(&square, &square_spaces)?;
    println!(
        "Square: {} nearest-neighbor edges",
        mapped.groups(&NearestNeighbors)?.count()
    );
    for axis in [Axis::X, Axis::Y] {
        println!(
            "Square: {} straight triples",
            mapped.groups(&AlongAxis::<3> { axis })?.count()
        );
    }
    let rule = ClosedNeighborhood;
    let center_group = mapped
        .groups(&rule)?
        .nth(4)
        .expect("center")
        .expect("valid mapping");
    assert!(std::ptr::eq(center_group.center(), &square_spaces[4]));
    assert_eq!(center_group.neighbors().count(), 4);
    println!(
        "Square: center space plus {} neighbors",
        center_group.neighbors().count()
    );

    let spin_sites = [
        HilbertSpace::tensor_all([HilbertSpace::factor(10), HilbertSpace::factor(11)]).unwrap(),
        HilbertSpace::tensor_all([HilbertSpace::factor(20), HilbertSpace::factor(21)]).unwrap(),
    ];
    let spin_chain = Chain::new(2, Boundary::Open);
    let spin_map = SiteMapper::new(&spin_chain, &spin_sites)?;
    let rule = NearestNeighbors;
    let group = spin_map.groups(&rule)?.next().unwrap()?;
    assert_eq!(group.leaves_cartesian().count(), 4);
    assert_eq!(group.leaves_zip().count(), 2);
    assert_eq!(
        group
            .combine_leaves(|mut streams| streams.remove(0))
            .count(),
        2
    );
    println!("Two two-factor sites: 4 Cartesian pairs, 2 zipped pairs");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn custom_relation_error_is_preserved_without_display_or_debug_bounds() {
        struct UserError;
        struct CustomRelation {
            fail: bool,
        }
        impl Relation<Chain> for CustomRelation {
            type Error = UserError;
            type Group<'a> = [usize; 1];

            fn groups<'a>(
                &'a self,
                geometry: &'a Chain,
            ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error> {
                if self.fail {
                    return Err(UserError);
                }
                Ok((0..geometry.num_sites()).map(|site| [site]))
            }
        }
        let chain = Chain::new(3, Boundary::Open);
        let sites: Vec<_> = (0..3).map(HilbertSpace::factor).collect();
        let mapper = SiteMapper::new(&chain, &sites).unwrap();
        assert!(matches!(
            mapper.groups(&CustomRelation { fail: true }),
            Err(UserError)
        ));
        assert_eq!(
            positions(&mapper, &CustomRelation { fail: false }),
            [[0], [1], [2]]
        );
    }

    fn leaf_ids<'a>(rows: impl Iterator<Item = Vec<&'a HilbertSpace<u8>>>) -> Vec<Vec<u8>> {
        rows.map(|row| {
            row.into_iter()
                .map(|leaf| *leaf.factor_index().unwrap())
                .collect()
        })
        .collect()
    }

    #[test]
    fn leaf_combiners_keep_tensor_order_for_uneven_nested_sites() {
        let sites = [
            HilbertSpace::tensor_all([
                HilbertSpace::factor(9u8),
                HilbertSpace::tensor_all([HilbertSpace::unit(), HilbertSpace::factor(3)]).unwrap(),
            ])
            .unwrap(),
            HilbertSpace::tensor_all([
                HilbertSpace::factor(8),
                HilbertSpace::factor(2),
                HilbertSpace::factor(7),
            ])
            .unwrap(),
            HilbertSpace::tensor_all([HilbertSpace::factor(6), HilbertSpace::factor(5)]).unwrap(),
        ];
        let chain = Chain::new(3, Boundary::Open);
        let mapped = SiteMapper::new(&chain, &sites).unwrap();
        let rule = Consecutive::<3>;
        let group = mapped.groups(&rule).unwrap().next().unwrap().unwrap();
        assert_eq!(
            leaf_ids(group.leaves_cartesian()),
            [
                [9, 8, 6],
                [9, 8, 5],
                [9, 2, 6],
                [9, 2, 5],
                [9, 7, 6],
                [9, 7, 5],
                [3, 8, 6],
                [3, 8, 5],
                [3, 2, 6],
                [3, 2, 5],
                [3, 7, 6],
                [3, 7, 5],
            ]
        );
        assert_eq!(leaf_ids(group.leaves_zip()), [[9, 8, 6], [3, 2, 5]]);
        let custom = group.combine_leaves(|mut streams| {
            streams.reverse();
            let mut done = false;
            std::iter::from_fn(move || {
                if done {
                    return None;
                }
                let row: Option<Vec<_>> = streams.iter_mut().map(Iterator::next).collect();
                if row.is_none() {
                    done = true;
                }
                row
            })
        });
        assert_eq!(leaf_ids(custom), [[6, 8, 9], [5, 2, 3]]);
        let expected: Vec<_> = sites
            .iter()
            .map(|site| site.leaves().collect::<Vec<_>>())
            .collect();
        for row in group.leaves_cartesian() {
            for (slot, leaf) in row.into_iter().enumerate() {
                assert!(
                    expected[slot]
                        .iter()
                        .any(|original| std::ptr::eq(*original, leaf))
                );
            }
        }
        assert_eq!(*group.source(), [0, 1, 2]);
    }

    #[test]
    fn leaf_combiners_handle_no_sites_empty_site_single_site_and_exhaustion() {
        let chain = Chain::new(1, Boundary::Open);
        let sites = [HilbertSpace::factor(4u8)];
        let mapped = SiteMapper::new(&chain, &sites).unwrap();
        let rule = TestGroups {
            groups: vec![vec![], vec![0]],
        };
        let groups: Vec<_> = mapped.groups(&rule).unwrap().map(Result::unwrap).collect();
        assert_eq!(leaf_ids(groups[0].leaves_cartesian()), vec![vec![]]);
        assert_eq!(leaf_ids(groups[0].leaves_zip()), vec![vec![]]);
        assert_eq!(leaf_ids(groups[1].leaves_cartesian()), vec![vec![4]]);
        assert_eq!(leaf_ids(groups[1].leaves_zip()), vec![vec![4]]);
        let sites = [HilbertSpace::<u8>::unit(), HilbertSpace::factor(7)];
        let chain = Chain::new(2, Boundary::Open);
        let mapped = SiteMapper::new(&chain, &sites).unwrap();
        let rule = NearestNeighbors;
        let group = mapped.groups(&rule).unwrap().next().unwrap().unwrap();
        let mut product = group.leaves_cartesian();
        assert!(product.next().is_none());
        assert!(product.next().is_none());
        let mut zipped = group.leaves_zip();
        assert!(zipped.next().is_none());
        assert!(zipped.next().is_none());
    }

    #[test]
    fn leaf_combination_needs_no_index_traits_and_accepts_fn_once() {
        struct Opaque;
        let sites = [HilbertSpace::factor(Opaque), HilbertSpace::factor(Opaque)];
        let chain = Chain::new(2, Boundary::Open);
        let mapped = SiteMapper::new(&chain, &sites).unwrap();
        let rule = NearestNeighbors;
        let group = mapped.groups(&rule).unwrap().next().unwrap().unwrap();
        assert_eq!(group.leaves_cartesian().count(), 1);
        assert_eq!(group.leaves_zip().count(), 1);
        let token = String::from("consumed by the combiner");
        let mut custom = group.combine_leaves(move |mut streams| {
            drop(token);
            streams.remove(1).map(|leaf| vec![leaf])
        });
        let row = custom.next().unwrap();
        assert!(std::ptr::eq(row[0], &sites[1]));
        assert!(custom.next().is_none());
    }

    // No Group bound and no bound on Index: mapping does not evaluate or compare spaces.
    fn positions<G: Geometry, R: Relation<G>, I>(
        map: &SiteMapper<'_, '_, G, HilbertSpace<I>>,
        relation: &R,
    ) -> Vec<Vec<usize>> {
        map.groups(relation)
            .unwrap_or_else(|_| panic!("relation failed"))
            .map(|group| {
                let group = group.unwrap();
                for (&index, space) in group.source().as_ref().iter().zip(group.iter()) {
                    assert!(std::ptr::eq(space, &map.values()[index]));
                }
                group.source().as_ref().to_vec()
            })
            .collect()
    }

    #[test]
    fn one_geometry_and_space_map_accept_multiple_relations() {
        struct Opaque;
        let chain = Chain::new(4, Boundary::Open);
        let spaces: Vec<_> = (0..4).map(|_| HilbertSpace::factor(Opaque)).collect();
        let map = SiteMapper::new(&chain, &spaces).unwrap();
        assert_eq!(positions(&map, &NearestNeighbors), [[0, 1], [1, 2], [2, 3]]);
        assert_eq!(positions(&map, &Consecutive::<3>), [[0, 1, 2], [1, 2, 3]]);
        let explicit = TestGroups {
            groups: vec![vec![], vec![2], vec![3, 0, 1]],
        };
        assert_eq!(positions(&map, &explicit), explicit.groups);
    }

    #[test]
    fn same_relation_type_works_on_both_geometries_and_axes_are_distinct() {
        let chain = Chain::new(4, Boundary::Periodic);
        let values: Vec<_> = (0..4).map(HilbertSpace::factor).collect();
        assert_eq!(
            positions(
                &SiteMapper::new(&chain, &values).unwrap(),
                &NearestNeighbors
            ),
            [[0, 1], [0, 3], [1, 2], [2, 3]]
        );

        let square = OpenSquare::new(3, 3).unwrap();
        let values: Vec<_> = (0..9).map(HilbertSpace::factor).collect();
        let map = SiteMapper::new(&square, &values).unwrap();
        assert_eq!(
            positions(&map, &NearestNeighbors),
            [
                [0, 1],
                [0, 3],
                [1, 2],
                [1, 4],
                [2, 5],
                [3, 4],
                [3, 6],
                [4, 5],
                [4, 7],
                [5, 8],
                [6, 7],
                [7, 8],
            ]
        );
        assert_eq!(
            positions(&map, &AlongAxis::<3> { axis: Axis::X }),
            [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        );
        assert_eq!(
            positions(&map, &AlongAxis::<3> { axis: Axis::Y }),
            [[0, 3, 6], [1, 4, 7], [2, 5, 8]]
        );
    }

    #[test]
    fn neighborhoods_keep_center_directions_and_variable_lengths() {
        let square = OpenSquare::new(3, 3).unwrap();
        let values: Vec<_> = (0..9).map(HilbertSpace::factor).collect();
        let map = SiteMapper::new(&square, &values).unwrap();
        let rule = ClosedNeighborhood;
        let groups: Vec<_> = map.groups(&rule).unwrap().map(Result::unwrap).collect();
        assert_eq!(
            groups.iter().map(|g| g.iter().count()).collect::<Vec<_>>(),
            [3, 4, 3, 4, 5, 4, 3, 4, 3]
        );
        for (center, group) in groups.iter().enumerate() {
            assert!(std::ptr::eq(group.center(), &values[center]));
        }
        let corner: Vec<_> = groups[0]
            .neighbors()
            .map(|(direction, space)| (direction, *space.factor_index().unwrap()))
            .collect();
        assert_eq!(corner, [(Direction::East, 1), (Direction::North, 3)]);
        let top_right: Vec<_> = groups[8]
            .neighbors()
            .map(|(direction, space)| (direction, *space.factor_index().unwrap()))
            .collect();
        assert_eq!(top_right, [(Direction::West, 7), (Direction::South, 5)]);
    }

    #[test]
    fn validation_empty_inputs_and_periodic_order() {
        let chain = Chain::new(2, Boundary::Periodic);
        let values = [0, 1];
        assert!(matches!(
            SiteMapper::new(&chain, &values[..1]),
            Err(MappingError::SiteCountMismatch { .. })
        ));
        let map = SiteMapper::new(&chain, &values).unwrap();
        assert!(matches!(
            map.groups(&Consecutive::<0>),
            Err(ConsecutiveError::ZeroLength)
        ));
        let groups: Vec<_> = map
            .groups(&Consecutive::<2>)
            .unwrap()
            .map(|g| *g.unwrap().source())
            .collect();
        assert_eq!(groups, [[0, 1], [1, 0]]);
        assert_eq!(map.groups(&NearestNeighbors).unwrap().count(), 1);
        assert_eq!(map.groups(&Consecutive::<3>).unwrap().count(), 0);
        let bad = TestGroups {
            groups: vec![vec![2]],
        };
        assert!(matches!(
            map.groups(&bad).unwrap().next().unwrap(),
            Err(MappingError::SiteOutOfRange { .. })
        ));
        let empty = OpenSquare::new(0, 3).unwrap();
        let empty_map = SiteMapper::new(&empty, &[] as &[u8]).unwrap();
        assert_eq!(empty_map.groups(&ClosedNeighborhood).unwrap().count(), 0);
        assert_eq!(empty_map.groups(&NearestNeighbors).unwrap().count(), 0);
        assert_eq!(
            empty_map
                .groups(&AlongAxis::<3> { axis: Axis::X })
                .unwrap()
                .count(),
            0
        );
        assert!(matches!(
            OpenSquare::new(usize::MAX, 2),
            Err(GridError::InvalidGridSize)
        ));
    }

    #[test]
    fn mapped_spaces_feed_existing_expression_construction_without_evaluation_bounds() {
        use qmb_spin_half::{SpinHalf, SpinHalfOperator as Op};
        #[derive(PartialEq, Eq)]
        struct Index(usize);
        struct Scale;
        impl From<f64> for Scale {
            fn from(_: f64) -> Self {
                Self
            }
        }
        let chain = Chain::new(3, Boundary::Open);
        let spaces: Vec<_> = (0..3).map(|i| HilbertSpace::factor(Index(i))).collect();
        let ambient = HilbertSpace::tensor_all(spaces.iter().cloned()).unwrap();
        let map = SiteMapper::new(&chain, &spaces).unwrap();
        let relation = TestGroups {
            groups: vec![vec![], vec![1], vec![2, 0]],
        };
        let words: Vec<_> = map
            .groups(&relation)
            .unwrap()
            .map(|group| {
                let group = group.unwrap();
                let factors = group
                    .iter()
                    .map(|factor| SpinHalf::<Scale>::operator(Op::Sz).scaled(2.0).at(factor));
                SpinHalf::<Scale>::try_operator_product(factors)
                    .unwrap()
                    .bind(&ambient)
                    .unwrap()
            })
            .collect();
        assert_eq!(words.len(), 3);
        assert!(words.iter().all(|word| word.space().is_some()));
    }
}
