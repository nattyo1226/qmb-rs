//! Spin model and product-state builders on immutable chain assignments.

use crate::{SpinHalfOperator, SpinHalfState};
use alloc::vec::Vec;
use core::ops::Neg;
use qmb_core::{ExpressionError, HilbertSpace, LocalOperator, LocalState, Operator, State};
use qmb_geometry::Chain;

/// Failure to build a spin model on a site assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum UtilityError {
    /// A site has zero or multiple spin factors.
    #[error("site {0} must have exactly one spin factor")]
    InvalidSite(usize),
    /// Site indices do not partition the supplied full space.
    #[error("site indices do not partition the full space")]
    IncompatibleSpace,
    /// Placement or binding failed.
    #[error("{0}")]
    Expression(#[from] ExpressionError),
}

// Validate the model-specific assignment and borrow the caller's actual leaves.
// Neither index cloning nor any relation between tensor grouping and site order
// is needed. The explicit space must cover exactly the chain's spin indices.
fn factors<'a, I: Eq>(
    chain: &Chain<I>,
    space: &'a HilbertSpace<I>,
) -> Result<Vec<&'a HilbertSpace<I>>, UtilityError> {
    for (i, site) in chain.sites().iter().enumerate() {
        if site.len() != 1 {
            return Err(UtilityError::InvalidSite(i));
        }
    }
    let mut leaves = Vec::new();
    let mut pending = alloc::vec![space];
    while let Some(space) = pending.pop() {
        if space.factor_index().is_some() {
            leaves.push(space);
        } else {
            pending.extend(space.children().rev());
        }
    }
    if leaves.len() != chain.sites().len() {
        return Err(UtilityError::IncompatibleSpace);
    }
    chain
        .sites()
        .iter()
        .map(|site| {
            leaves
                .iter()
                .copied()
                .find(|factor| factor.factor_index() == Some(&site[0]))
                .ok_or(UtilityError::IncompatibleSpace)
        })
        .collect()
}

/// Builds `-sum_(i,j) (Jx Sx_i Sx_j + Jy Sy_i Sy_j + Jz Sz_i Sz_j)`.
///
/// Couplings are real. Each site must contain exactly one factor. The result is
/// bound to the supplied full space, which must contain exactly the chain's
/// indices. Existing grouping and factor order are preserved; no indices are issued.
pub fn xyz<I: Eq, R, C: Clone + From<R> + Neg<Output = C>>(
    chain: &Chain<I>,
    space: &HilbertSpace<I>,
    jx: R,
    jy: R,
    jz: R,
) -> Result<Operator<I, SpinHalfOperator, C>, UtilityError> {
    let factors = factors(chain, space)?;
    let axes = [
        (SpinHalfOperator::Sx, -C::from(jx)),
        (SpinHalfOperator::Sy, -C::from(jy)),
        (SpinHalfOperator::Sz, -C::from(jz)),
    ];
    let mut terms = Vec::new();
    for (left, right) in chain.nn_site_pairs() {
        let left = factors[left];
        let right = factors[right];
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
    Ok(Operator::sum(terms)?.bind(space)?)
}

/// Builds the isotropic XYZ model `-J sum_(i,j) (Sx_i Sx_j + Sy_i Sy_j + Sz_i Sz_j)`.
///
/// Delegates to [`xyz`] with `Jx = Jy = Jz = coupling`; positive coupling is
/// ferromagnetic. The same site validation and full-space binding apply.
pub fn heisenberg<I: Eq, R: Clone, C: Clone + From<R> + Neg<Output = C>>(
    chain: &Chain<I>,
    space: &HilbertSpace<I>,
    coupling: R,
) -> Result<Operator<I, SpinHalfOperator, C>, UtilityError> {
    xyz(chain, space, coupling.clone(), coupling.clone(), coupling)
}

/// Builds an alternating product state, starting with the specified orientation.
///
/// Odd periodic chains are allowed. An empty chain yields the unit state.
/// The full-space factor order is preserved even if it differs from site order.
pub fn neel<I: Eq, C>(
    chain: &Chain<I>,
    space: &HilbertSpace<I>,
    first: SpinHalfState,
) -> Result<State<I, SpinHalfState, SpinHalfOperator, C>, UtilityError> {
    let factors = factors(chain, space)?;
    let states = space
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
