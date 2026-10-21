//! Static associations between a local structure and its primitives.

use core::num::NonZeroUsize;

use crate::{ActOn, Adjoint, InnerProduct, LinearCombination, Product, Trace};

/// A fixed, finite-dimensional local structure at coefficient precision `T`.
///
/// The structure has no runtime definition parameters. For a given implementing
/// type, `State`, `Operator`, and `DIMENSION` must remain the same across all
/// supported `T`. This is an implementer contract, not enforced by the type
/// system. Primitive types need not contain a numeric type parameter.
///
/// Concrete rules must respect linearity, product associativity, the action of
/// products, the identity laws, the inner product and adjoint, and ordinary
/// trace. Primitives may be linearly dependent; rules must preserve their linear
/// relations. No blanket cloning, comparison, or thread-safety bound is imposed.
pub trait LocalStructure<T> {
    /// The state primitive type for this structure.
    type State: InnerProduct<T>;
    /// The operator primitive type and its required local capabilities.
    type Operator: Product<T> + ActOn<Self::State, T> + Adjoint<T> + Trace<T>;

    /// The positive dimension of the local state space.
    const DIMENSION: NonZeroUsize;

    /// Returns the identity, which need not be a single primitive.
    fn identity() -> LinearCombination<Self::Operator, T>;
}
