//! Shared syntax parts and immutable local and composite expressions.

pub mod parts;

pub mod composite;
mod composite_eval;
mod evaluate;
mod inner_product;
mod local;
mod operators;

pub use composite::{ExpressionError, Operator, State};
pub use local::{LocalOperator, LocalState};

#[cfg(test)]
mod tests;

#[cfg(test)]
mod composite_tests;
