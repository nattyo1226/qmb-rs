//! Spin-conserving Hubbard chains on caller-owned two-mode site spaces.

use crate::{Fermion, FermionOperator as O};
use alloc::vec::Vec;
use core::ops::Neg;
use qmb_core::{ExpressionError, HilbertSpace, OperatorOf, QuantumSystem, SpaceError, TensorKind};
use qmb_geometry::{Chain, ChainError};

/// Failure to construct a Hubbard model on a site-ordered space assignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum UtilityError {
    /// The supplied site count does not match the chain.
    #[error("{0}")]
    Chain(#[from] ChainError),
    /// Indices are duplicated, or tensor kinds are incompatible.
    #[error("{0}")]
    Space(#[from] SpaceError),
    /// Placing or binding the expression failed.
    #[error("{0}")]
    Expression(#[from] ExpressionError),
    /// A site is not an effectively graded space with exactly two leaves.
    #[error("site {site} must contain exactly two graded modes in up, down order")]
    InvalidSite {
        /// Position in the supplied site sequence.
        site: usize,
    },
}

impl<C> Fermion<C> {
    /// Builds `-t Σ_(i,j),σ (f†_iσ f_jσ + f†_jσ f_iσ) + U Σ_i n_i↑ n_i↓`.
    ///
    /// Each site must have exactly two leaves in `(up, down)` order and an
    /// effective graded tensor kind. The result is bound to a graded tensor
    /// sharing those sites in input order. The chain supplies each undirected
    /// nearest-neighbor edge once, including for two-site periodic chains.
    ///
    /// Parameters must represent real values for this Hermite Hamiltonian.
    /// Construction requires no local evaluation, fractions, or imaginary unit.
    /// Site count, shape, kind and global index uniqueness are validated even
    /// for zero coefficients or chains without edges. Zero sites return zero
    /// on the graded unit space; one site retains all on-site terms.
    ///
    /// ```
    /// use num_complex::Complex;
    /// use qmb_core::HilbertSpace;
    /// use qmb_geometry::{Boundary, Chain};
    /// use qmb_fermion::Fermion;
    /// type F = Fermion<Complex<f64>>;
    /// let site = HilbertSpace::tensor_g([HilbertSpace::factor("up"), HilbertSpace::factor("down")])?;
    /// let h = F::hubbard(&Chain::new(1, Boundary::Open), &[site], 1., 4.)?;
    /// assert_eq!(h.trace()?, Complex::new(4., 0.));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn hubbard<I: Eq, R>(
        chain: &Chain,
        sites: &[HilbertSpace<I>],
        t: R,
        u: R,
    ) -> Result<OperatorOf<Self, I>, UtilityError>
    where
        C: Clone + From<R> + Neg<Output = C>,
    {
        Self::hubbard_with_interaction(chain, sites, C::from(t), C::from(u), O::Number)
    }

    /// Builds the Hubbard model with interaction `U (n_i↑-I/2)(n_i↓-I/2)`.
    ///
    /// Only the interaction changes from [`Self::hubbard`]. No constant energy
    /// is discarded. For L sites and total particle number N,
    /// `H_sym(t,U) = H(t,U) - (U/2) N + UL/4 I`.
    ///
    /// Uses `U Z_i↑ Z_i↓`, since `Z=I/2-n`, to construct the shifted interaction
    /// without requiring numeric fractions. Validation, site order, boundaries,
    /// parameter conventions and errors are the same as [`Self::hubbard`].
    pub fn symmetric_hubbard<I: Eq, R>(
        chain: &Chain,
        sites: &[HilbertSpace<I>],
        t: R,
        u: R,
    ) -> Result<OperatorOf<Self, I>, UtilityError>
    where
        C: Clone + From<R> + Neg<Output = C>,
    {
        Self::hubbard_with_interaction(chain, sites, C::from(t), C::from(u), O::MajoranaZ)
    }

    fn hubbard_with_interaction<I: Eq>(
        chain: &Chain,
        sites: &[HilbertSpace<I>],
        t: C,
        u: C,
        interaction: O,
    ) -> Result<OperatorOf<Self, I>, UtilityError>
    where
        C: Clone + Neg<Output = C>,
    {
        // Validate without depending on edges, coefficients, or evaluator traits.
        let _ = chain.nn_site_pairs_with(sites)?;
        let mut modes = Vec::with_capacity(sites.len());
        for (site, space) in sites.iter().enumerate() {
            let mut leaves = space.leaves();
            let pair = match (leaves.next(), leaves.next(), leaves.next()) {
                (Some(up), Some(down), None)
                    if space.effective_kind() == Some(TensorKind::Graded) =>
                {
                    [up, down]
                }
                _ => return Err(UtilityError::InvalidSite { site }),
            };
            modes.push(pair);
        }
        let space = HilbertSpace::tensor_g(sites.iter().cloned())?;
        let mut terms = Vec::new();
        let minus_t = -t;
        for [up, down] in &modes {
            terms.push(
                Self::operator_product_at([(*up, interaction), (*down, interaction)])?
                    .scaled(u.clone()),
            );
        }
        for (left, right) in chain.nn_site_pairs_with(&modes)? {
            for spin in 0..2 {
                for (from, to) in [(left[spin], right[spin]), (right[spin], left[spin])] {
                    terms.push(
                        Self::operator_product_at([(to, O::Create), (from, O::Annihilate)])?
                            .scaled(minus_t.clone()),
                    );
                }
            }
        }
        Ok(Self::operator_sum(terms)?.bind(&space)?)
    }
}
