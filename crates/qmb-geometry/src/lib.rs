//! Pure chain geometry with borrowing connections to caller-owned spaces.

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod boundary;
mod connection;
mod geometry;
mod relation;

pub use boundary::Boundary;
pub use connection::ChainError;
pub use geometry::chain::Chain;
