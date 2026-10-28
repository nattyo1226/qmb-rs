//! Shared syntax parts and immutable local and composite expressions.

pub mod parts;

mod arithmetic;
pub mod composite;
pub(crate) mod local;
pub(crate) mod traversal;

pub use composite::{ExpressionError, Operator, State};
pub use local::{LocalOperator, LocalState};

#[cfg(test)]
mod tests;

#[cfg(test)]
mod composite_tests;
