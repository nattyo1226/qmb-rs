//! Spin models built on a chain and caller-owned site spaces.

use crate::{SpinHalf, SpinHalfOperator, SpinHalfState};
use alloc::vec::Vec;
use core::ops::Neg;
use qmb_core::{ExpressionError, HilbertSpace, OperatorOf, SpaceError, State, StateOf};
use qmb_geometry::{Chain, MappingError, NearestNeighbors, SiteMapper};

/// Failure to construct a spin model on a site-ordered space assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum UtilityError {
    /// The input site count differs from the chain's count.
    #[error("{0}")]
    Mapping(#[from] MappingError),
    /// An index occurs more than once across the supplied site spaces.
    #[error("{0}")]
    Space(#[from] SpaceError),
    /// Placement or binding failed.
    #[error("{0}")]
    Expression(#[from] ExpressionError),
}

// Validate independently of edge enumeration, sharing the existing site trees.
fn model_space<I: Eq>(
    chain: &Chain,
    factors: &[HilbertSpace<I>],
) -> Result<HilbertSpace<I>, UtilityError> {
    SiteMapper::new(chain, factors)?;
    Ok(HilbertSpace::tensor_all(factors.iter().cloned())?)
}

/// Model constructors for the selected SpinHalf system.
///
/// No values are stored; select the coefficient type through `SpinHalf<C>`.
pub struct SpinHalfModels<Q>(core::marker::PhantomData<fn() -> Q>);

impl<C> SpinHalfModels<SpinHalf<C>> {
    /// Builds `-sum_(i,j) sum_(a in i,b in j) (Jx Sx_a Sx_b + Jy Sy_a Sy_b + Jz Sz_a Sz_b)`.
    ///
    /// Each input space is one site, with any number of leaves. All leaf pairs
    /// across nearest-neighbor sites interact; no within-site terms are added.
    /// Empty sites retain their positions without reconnecting neighbors.
    /// The exact site count and index uniqueness are checked even without edges.
    /// The result is bound to a tensor sharing the supplied site trees in order.
    /// Zero sites return a zero Hamiltonian bound to the unit space.
    pub fn xyz<I: Eq, R>(
        chain: &Chain,
        factors: &[HilbertSpace<I>],
        jx: R,
        jy: R,
        jz: R,
    ) -> Result<OperatorOf<SpinHalf<C>, I>, UtilityError>
    where
        C: Clone + From<R> + Neg<Output = C>,
    {
        let space = model_space(chain, factors)?;
        let axes = [
            (SpinHalfOperator::Sx, -C::from(jx)),
            (SpinHalfOperator::Sy, -C::from(jy)),
            (SpinHalfOperator::Sz, -C::from(jz)),
        ];
        let mut terms = Vec::new();
        let mapper = SiteMapper::new(chain, factors)?;
        for group in mapper
            .groups(&NearestNeighbors)
            .unwrap_or_else(|never| match never {})
        {
            let group = group?;
            for pair in group.leaves_cartesian() {
                let (left, right) = (pair[0], pair[1]);
                for (primitive, coupling) in &axes {
                    let local = SpinHalf::<C>::operator(*primitive);
                    terms.push((local.at(left)? * local.at(right)?)?.scaled(coupling.clone()));
                }
            }
        }
        Ok(OperatorOf::<SpinHalf<C>, I>::sum(terms)?.bind(&space)?)
    }

    /// Builds isotropic XYZ, with positive coupling ferromagnetic.
    /// Uses the same validation, all-leaf pairing and space binding as [`Self::xyz`].
    pub fn heisenberg<I: Eq, R: Clone>(
        chain: &Chain,
        factors: &[HilbertSpace<I>],
        coupling: R,
    ) -> Result<OperatorOf<SpinHalf<C>, I>, UtilityError>
    where
        C: Clone + From<R> + Neg<Output = C>,
    {
        Self::xyz(chain, factors, coupling.clone(), coupling.clone(), coupling)
    }

    /// Builds a product state alternating orientation by site, starting at `first`.
    ///
    /// All leaves of a site have the same orientation, and empty sites count
    /// toward parity. Leaves appear in site order and then local leaf order.
    /// Zero sites or all-empty sites give the unit state. Odd periodic chains
    /// are allowed. No numeric or cloning capabilities are required of `C`.
    pub fn neel<I: Eq>(
        chain: &Chain,
        factors: &[HilbertSpace<I>],
        first: SpinHalfState,
    ) -> Result<StateOf<SpinHalf<C>, I>, UtilityError> {
        let space = model_space(chain, factors)?;
        let mut states = Vec::new();
        for (i, site) in factors.iter().enumerate() {
            let orientation = if i % 2 == 0 {
                first
            } else {
                match first {
                    SpinHalfState::Up => SpinHalfState::Down,
                    SpinHalfState::Down => SpinHalfState::Up,
                }
            };
            for leaf in site.leaves() {
                states.push(SpinHalf::<C>::state(orientation).at(leaf)?);
            }
        }
        Ok(State::tensor_in(&space, states)?)
    }
}
