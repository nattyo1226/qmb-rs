//! Spin-conserving Hubbard chains on caller-owned two-mode site spaces.

use crate::coefficient::quarter;
use crate::{Fermion, FermionOperator as O};
use alloc::vec::Vec;
use core::ops::{Add, Div, Neg};
use num_traits::One;
use qmb_core::{ExpressionError, HilbertSpace, OperatorOf, SpaceError, TensorKind};
use qmb_geometry::{Chain, MappingError, NearestNeighbors, SiteMapper};

/// Failure to construct a Hubbard model on a site-ordered space assignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum UtilityError {
    /// The supplied site count does not match the chain.
    #[error("{0}")]
    Mapping(#[from] MappingError),
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

/// Model constructors for the selected Fermion system.
///
/// No values are stored; select the coefficient type through `Fermion<C>`.
pub struct FermionModels<Q>(core::marker::PhantomData<fn() -> Q>);

impl<C> FermionModels<Fermion<C>> {
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
    /// let site = HilbertSpace::tensor_all_g([HilbertSpace::factor("up"), HilbertSpace::factor("down")])?;
    /// let h = qmb_fermion::models::FermionModels::<F>::hubbard(&Chain::new(1, Boundary::Open), &[site], 1., 4.)?;
    /// assert_eq!(h.trace()?, Complex::new(4., 0.));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn hubbard<I: Eq, R>(
        chain: &Chain,
        sites: &[HilbertSpace<I>],
        t: R,
        u: R,
    ) -> Result<OperatorOf<Fermion<C>, I>, UtilityError>
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
    /// Uses `(U/4) Mz_i↑ Mz_i↓`, since `Mz=I-2n`. The coefficient type must
    /// represent mathematical quarters. Validation, site order, boundaries,
    /// parameter conventions and errors are the same as [`Self::hubbard`].
    pub fn symmetric_hubbard<I: Eq, R>(
        chain: &Chain,
        sites: &[HilbertSpace<I>],
        t: R,
        u: R,
    ) -> Result<OperatorOf<Fermion<C>, I>, UtilityError>
    where
        C: Clone + From<R> + Neg<Output = C> + One + Add<Output = C> + Div<Output = C>,
    {
        Self::hubbard_with_interaction(chain, sites, C::from(t), C::from(u) * quarter::<C>(), O::Mz)
    }

    fn hubbard_with_interaction<I: Eq>(
        chain: &Chain,
        sites: &[HilbertSpace<I>],
        t: C,
        u: C,
        interaction: O,
    ) -> Result<OperatorOf<Fermion<C>, I>, UtilityError>
    where
        C: Clone + Neg<Output = C>,
    {
        // Validate without depending on edges, coefficients, or evaluator traits.
        SiteMapper::new(chain, sites)?;
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
        let space = HilbertSpace::tensor_all_g(sites.iter().cloned())?;
        let mut terms = Vec::new();
        let minus_t = -t;
        for [up, down] in &modes {
            terms.push(
                Fermion::<C>::operator_product_at([(*up, interaction), (*down, interaction)])?
                    .scaled(u.clone()),
            );
        }
        let mapper = SiteMapper::new(chain, &modes)?;
        for group in mapper
            .groups(&NearestNeighbors)
            .unwrap_or_else(|never| match never {})
        {
            let group = group?;
            let left = group.get(0).expect("nearest-neighbor pair");
            let right = group.get(1).expect("nearest-neighbor pair");
            for spin in 0..2 {
                for (from, to) in [(left[spin], right[spin]), (right[spin], left[spin])] {
                    terms.push(
                        Fermion::<C>::operator_product_at([
                            (to, O::Create),
                            (from, O::Annihilate),
                        ])?
                        .scaled(minus_t.clone()),
                    );
                }
            }
        }
        Ok(Fermion::<C>::operator_sum(terms)?.bind(&space)?)
    }
}
