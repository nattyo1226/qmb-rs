//! Spin model and product-state builders on immutable chain assignments.

use crate::{SpinHalfOperator, SpinHalfState};
use alloc::vec::Vec;
use core::ops::Neg;
use qmb_core::{ExpressionError, HilbertSpace, LocalOperator, LocalState, Operator, State};
use qmb_geometry::{Chain, NearestNeighborRelation};

/// Failure to build a spin model on a site assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum UtilityError {
    /// A site has zero or multiple spin factors.
    #[error("site {0} must have exactly one spin factor")]
    InvalidSite(usize),
    /// Placement or binding failed.
    #[error("{0}")]
    Expression(#[from] ExpressionError),
}

// Borrow the leaf itself, preserving existing identity and Arc sharing, including
// a single factor wrapped in nested or unit-containing tensors.
fn factors<I>(chain: &Chain<I>) -> Result<Vec<&HilbertSpace<I>>, UtilityError> {
    chain
        .sites()
        .iter()
        .enumerate()
        .map(|(i, site)| {
            if site.indices().count() != 1 {
                return Err(UtilityError::InvalidSite(i));
            }
            let mut pending = alloc::vec![site];
            while let Some(space) = pending.pop() {
                if space.factor_index().is_some() {
                    return Ok(space);
                }
                pending.extend(space.children());
            }
            unreachable!("one indexed factor")
        })
        .collect()
}

/// Builds `-sum_(i,j) (Jx Sx_i Sx_j + Jy Sy_i Sy_j + Jz Sz_i Sz_j)`.
///
/// Couplings are real. Each site must contain exactly one factor. The result is
/// bound to the chain's full space, and no indices are issued.
pub fn xyz<I: Eq, R, C: Clone + From<R> + Neg<Output = C>>(
    chain: &Chain<I>,
    jx: R,
    jy: R,
    jz: R,
) -> Result<Operator<I, SpinHalfOperator, C>, UtilityError> {
    let factors = factors(chain)?;
    let axes = [
        (SpinHalfOperator::Sx, -C::from(jx)),
        (SpinHalfOperator::Sy, -C::from(jy)),
        (SpinHalfOperator::Sz, -C::from(jz)),
    ];
    let mut terms = Vec::new();
    for (left, right) in chain.nn_pairs() {
        let left = factors
            .iter()
            .find(|f| f.factor_index() == Some(left))
            .expect("chain index");
        let right = factors
            .iter()
            .find(|f| f.factor_index() == Some(right))
            .expect("chain index");
        for (primitive, coupling) in &axes {
            let local = LocalOperator::primitive(*primitive);
            terms.push(
                local
                    .at(left)?
                    .product(&local.at(right)?)?
                    .scaled(coupling.clone()),
            );
        }
    }
    Ok(Operator::sum(terms)?.bind(chain.space())?)
}

/// Builds the isotropic XYZ model `-J sum_(i,j) (Sx_i Sx_j + Sy_i Sy_j + Sz_i Sz_j)`.
///
/// Delegates to [`xyz`] with `Jx = Jy = Jz = coupling`; positive coupling is
/// ferromagnetic. The same site validation and full-space binding apply.
pub fn heisenberg<I: Eq, R: Clone, C: Clone + From<R> + Neg<Output = C>>(
    chain: &Chain<I>,
    coupling: R,
) -> Result<Operator<I, SpinHalfOperator, C>, UtilityError> {
    xyz(chain, coupling.clone(), coupling.clone(), coupling)
}

/// Builds an alternating product state, starting with the specified orientation.
///
/// Odd periodic chains are allowed. An empty chain yields the unit state.
/// The full-space factor order is preserved even if it differs from site order.
pub fn neel<I: Eq, C>(
    chain: &Chain<I>,
    first: SpinHalfState,
) -> Result<State<I, SpinHalfState, SpinHalfOperator, C>, UtilityError> {
    let factors = factors(chain)?;
    let states = chain
        .space()
        .indices()
        .map(|index| {
            let i = factors
                .iter()
                .position(|f| f.factor_index() == Some(index))
                .expect("partitioned space");
            let orientation = if i % 2 == 0 {
                first
            } else {
                match first {
                    SpinHalfState::Up => SpinHalfState::Down,
                    SpinHalfState::Down => SpinHalfState::Up,
                }
            };
            LocalState::primitive(orientation).at(factors[i])
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(State::tensor(states)?)
}
