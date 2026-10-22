//! Coefficient-independent associations for a fixed local structure.

use crate::LinearCombination;
use core::num::NonZeroUsize;

/// A static association of state and operator primitives with a finite dimension.
///
/// Metadata requires no numeric, cloning, or algebra capabilities. Operations
/// request those capabilities where they are used. Concrete rules must preserve
/// linear relations, associativity, identity, action, inner-product/adjoint
/// compatibility, and the ordinary unnormalized trace.
pub trait LocalStructure {
    /// State primitive type.
    type State;
    /// Operator primitive type.
    type Operator;
    /// Positive dimension of the local state space.
    const DIMENSION: NonZeroUsize;
}

/// A local structure providing its identity at coefficient precision `T`.
pub trait Identity<T>: LocalStructure {
    /// Returns the identity as a finite combination of operator primitives.
    fn identity() -> LinearCombination<Self::Operator, T>;
}
