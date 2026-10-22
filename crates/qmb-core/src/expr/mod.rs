//! Shared syntax parts and immutable local expressions.

pub mod parts;

mod evaluate;
mod local;
mod operators;

pub use local::{LocalOperator, LocalState};

#[cfg(test)]
mod tests;
