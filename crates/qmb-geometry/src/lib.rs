//! Geometric site relations and borrowing connections to caller-owned data.
//!
//! Geometry retains site counts and boundary conditions; [`Relation`] selects
//! ordered groups of site numbers. [`SiteMapper`] connects those groups to
//! a caller-owned slice without cloning its elements. For Hilbert spaces,
//! [`MappedGroup::leaves_cartesian`] and [`MappedGroup::leaves_zip`] combine
//! local factors while preserving site and tensor order.
//!
//! ```
//! use qmb_geometry::{Boundary, Chain, Consecutive, SiteMapper};
//! let chain = Chain::new(4, Boundary::Open);
//! let values = ["a", "b", "c", "d"];
//! let mapper = SiteMapper::new(&chain, &values)?;
//! let triples = mapper.groups(&Consecutive::<3>)?
//!     .map(|group| group.map(|group| group.iter().copied().collect::<Vec<_>>()))
//!     .collect::<Result<Vec<_>, _>>()?;
//! assert_eq!(triples, [["a", "b", "c"], ["b", "c", "d"]]);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

extern crate alloc;

mod boundary;
mod connection;
mod geometry;
mod group;
mod relation;

pub use boundary::Boundary;
pub use connection::{MappedGroup, MappedGroupResult, MappingError, SiteMapper};
pub use geometry::Geometry;
pub use geometry::chain::Chain;
pub use group::{SiteGroup, SiteGroupError};
pub use relation::{Consecutive, ConsecutiveError, NearestNeighbors, Relation};
