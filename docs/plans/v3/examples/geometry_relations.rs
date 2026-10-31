//! Adopted design prototype. No production geometry or API is added by this file.
//! The open square grid below exists only to exercise the interfaces.

use qmb_core::HilbertSpace;
use qmb_core::space::Leaves;
use qmb_geometry::{Boundary, Chain};

#[derive(Debug, PartialEq, Eq)]
enum Error {
    InvalidLength,
    InvalidGridSize,
    SiteCountMismatch,
    OutOfRange,
}

trait Geometry {
    /// Number of positions in a caller-owned, site-ordered slice.
    fn num_sites(&self) -> usize;
}

impl Geometry for Chain {
    fn num_sites(&self) -> usize {
        // Adapter to the current library; the adopted API will use num_sites.
        Chain::site_count(self)
    }
}

trait Relation<G: Geometry> {
    type Group<'a>: AsRef<[usize]>
    where
        Self: 'a,
        G: 'a;

    fn groups<'a>(
        &'a self,
        geometry: &'a G,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Error>;
}

struct SiteMapper<'g, 's, G, T> {
    geometry: &'g G,
    values: &'s [T],
}

impl<'g, 's, G: Geometry, T> SiteMapper<'g, 's, G, T> {
    fn new(geometry: &'g G, values: &'s [T]) -> Result<Self, Error> {
        if geometry.num_sites() != values.len() {
            return Err(Error::SiteCountMismatch);
        }
        Ok(Self { geometry, values })
    }

    fn groups<'a, R: Relation<G>>(
        &'a self,
        relation: &'a R,
    ) -> Result<impl Iterator<Item = Result<MappedGroup<'s, T, R::Group<'a>>, Error>> + 'a, Error>
    {
        let groups = relation.groups(self.geometry)?;
        let values = self.values;
        Ok(groups.map(move |source| {
            if source.as_ref().iter().any(|&site| site >= values.len()) {
                return Err(Error::OutOfRange);
            }
            Ok(MappedGroup { values, source })
        }))
    }
}

struct MappedGroup<'s, T, Group> {
    values: &'s [T],
    source: Group,
}

impl<'s, T, Group: AsRef<[usize]>> MappedGroup<'s, T, Group> {
    fn source(&self) -> &Group {
        &self.source
    }

    fn iter(&self) -> impl Iterator<Item = &'s T> + '_ {
        let values = self.values;
        self.source.as_ref().iter().map(move |&site| &values[site])
    }
}

impl<'s, I, Group: AsRef<[usize]>> MappedGroup<'s, HilbertSpace<I>, Group> {
    fn combine_leaves<Out: Iterator>(
        &self,
        combine: impl FnOnce(Vec<Leaves<'s, I>>) -> Out,
    ) -> Out {
        combine(self.iter().map(HilbertSpace::leaves).collect())
    }

    fn leaves_cartesian(&self) -> impl Iterator<Item = Vec<&'s HilbertSpace<I>>> {
        self.combine_leaves(cartesian_leaves)
    }

    fn leaves_zip(&self) -> impl Iterator<Item = Vec<&'s HilbertSpace<I>>> {
        self.combine_leaves(zip_leaves)
    }
}

// Cache borrowed leaves for replay: Leaves itself need not implement Clone.
// No Eq/Clone bound on Index; the cached items are shared references.
fn cartesian_leaves<I>(
    mut streams: Vec<Leaves<'_, I>>,
) -> impl Iterator<Item = Vec<&HilbertSpace<I>>> {
    let mut cache: Vec<Vec<_>> = (0..streams.len()).map(|_| Vec::new()).collect();
    let mut positions = vec![0; streams.len()];
    let mut started = false;
    let mut done = false;
    core::iter::from_fn(move || {
        if done {
            return None;
        }
        if !started {
            for (stream, cached) in streams.iter_mut().zip(&mut cache) {
                let Some(first) = stream.next() else {
                    done = true;
                    return None;
                };
                cached.push(first);
            }
            started = true;
        } else {
            let mut advanced = false;
            for slot in (0..streams.len()).rev() {
                let next = positions[slot] + 1;
                if next < cache[slot].len() {
                    positions[slot] = next;
                    advanced = true;
                    break;
                }
                if let Some(leaf) = streams[slot].next() {
                    cache[slot].push(leaf);
                    positions[slot] = next;
                    advanced = true;
                    break;
                }
                positions[slot] = 0;
            }
            if !advanced {
                done = true;
                return None;
            }
        }
        Some(
            cache
                .iter()
                .zip(&positions)
                .map(|(leaves, &pos)| leaves[pos])
                .collect(),
        )
    })
}

fn zip_leaves<I>(mut streams: Vec<Leaves<'_, I>>) -> impl Iterator<Item = Vec<&HilbertSpace<I>>> {
    let mut done = false;
    core::iter::from_fn(move || {
        if done {
            return None;
        }
        if streams.is_empty() {
            done = true;
            return Some(Vec::new());
        }
        let row: Option<Vec<_>> = streams.iter_mut().map(Iterator::next).collect();
        if row.is_none() {
            done = true;
        }
        row
    })
}

struct NearestNeighbors;

impl Relation<Chain> for NearestNeighbors {
    type Group<'a> = [usize; 2];

    fn groups<'a>(
        &'a self,
        geometry: &'a Chain,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Error> {
        Ok(geometry.nn_site_pairs().map(|(a, b)| [a, b]))
    }
}

struct Consecutive<const K: usize>;

impl<const K: usize> Relation<Chain> for Consecutive<K> {
    type Group<'a> = [usize; K];

    fn groups<'a>(
        &'a self,
        geometry: &'a Chain,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Error> {
        if K == 0 {
            return Err(Error::InvalidLength);
        }
        let n = geometry.num_sites();
        let count = if K > n {
            0
        } else if geometry.boundary() == Boundary::Periodic {
            n
        } else {
            n - K + 1
        };
        Ok((0..count).map(move |start| {
            core::array::from_fn(|offset| {
                // Both start and offset are below n. Avoid start + offset overflow.
                if offset >= n - start {
                    offset - (n - start)
                } else {
                    start + offset
                }
            })
        }))
    }
}

struct ExplicitRelations {
    groups: Vec<Vec<usize>>,
}

impl<G: Geometry> Relation<G> for ExplicitRelations {
    type Group<'a>
        = &'a [usize]
    where
        Self: 'a,
        G: 'a;

    fn groups<'a>(
        &'a self,
        _geometry: &'a G,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Error> {
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
    fn new(width: usize, height: usize) -> Result<Self, Error> {
        let count = width.checked_mul(height).ok_or(Error::InvalidGridSize)?;
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
    type Group<'a> = [usize; 2];

    fn groups<'a>(
        &'a self,
        geometry: &'a OpenSquare,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Error> {
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
    type Group<'a> = [usize; K];

    fn groups<'a>(
        &'a self,
        geometry: &'a OpenSquare,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Error> {
        if K == 0 {
            return Err(Error::InvalidLength);
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
    type Group<'a> = Neighborhood;

    fn groups<'a>(
        &'a self,
        geometry: &'a OpenSquare,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Error> {
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

// Optional, concrete group capabilities. Generic consumers still need only Relation<G>.
impl<'s, T> MappedGroup<'s, T, Neighborhood> {
    fn center(&self) -> &'s T {
        &self.values[self.source.sites[0]]
    }

    fn neighbors(&self) -> impl Iterator<Item = (Direction, &'s T)> + '_ {
        let values = self.values;
        self.source.directions[1..self.source.len]
            .iter()
            .zip(&self.source.sites[1..self.source.len])
            .map(move |(direction, &site)| (direction.expect("neighbor direction"), &values[site]))
    }
}

fn main() -> Result<(), Error> {
    let chain = Chain::new(5, Boundary::Open);
    let chain_spaces: Vec<_> = (0..5).map(HilbertSpace::factor).collect();
    let mapped = SiteMapper::new(&chain, &chain_spaces)?;
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

    let explicit = ExplicitRelations {
        groups: vec![vec![], vec![0], vec![0, 2, 4]],
    };
    assert_eq!(mapped.groups(&explicit)?.count(), 3);

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
    assert_eq!(group.combine_leaves(zip_leaves).count(), 2);
    println!("Two two-factor sites: 4 Cartesian pairs, 2 zipped pairs");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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
            zip_leaves(streams)
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
        let rule = ExplicitRelations {
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
            .unwrap()
            .map(|group| {
                let group = group.unwrap();
                for (&index, space) in group.source().as_ref().iter().zip(group.iter()) {
                    assert!(std::ptr::eq(space, &map.values[index]));
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
        let explicit = ExplicitRelations {
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
            Err(Error::SiteCountMismatch)
        ));
        let map = SiteMapper::new(&chain, &values).unwrap();
        assert!(matches!(
            map.groups(&Consecutive::<0>),
            Err(Error::InvalidLength)
        ));
        let groups: Vec<_> = map
            .groups(&Consecutive::<2>)
            .unwrap()
            .map(|g| *g.unwrap().source())
            .collect();
        assert_eq!(groups, [[0, 1], [1, 0]]);
        assert_eq!(map.groups(&NearestNeighbors).unwrap().count(), 1);
        assert_eq!(map.groups(&Consecutive::<3>).unwrap().count(), 0);
        let bad = ExplicitRelations {
            groups: vec![vec![2]],
        };
        assert!(matches!(
            map.groups(&bad).unwrap().next().unwrap(),
            Err(Error::OutOfRange)
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
            Err(Error::InvalidGridSize)
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
        let relation = ExplicitRelations {
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
