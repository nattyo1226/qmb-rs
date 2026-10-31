//! Spin-1 chain Hamiltonians on caller-owned site spaces.

use crate::{SpinOne, SpinOneOperator as O};
use core::ops::Neg;
use qmb_core::{ExpressionError, HilbertSpace, OperatorOf, SpaceError};
use qmb_geometry::{Chain, MappingError, NearestNeighbors, SiteMapper};

/// Failure to construct a spin-1 model on a site-ordered space assignment.
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

// Validate independently of edges and coefficients, sharing the site trees.
fn model_space<I: Eq>(
    chain: &Chain,
    sites: &[HilbertSpace<I>],
) -> Result<HilbertSpace<I>, UtilityError> {
    SiteMapper::new(chain, sites)?;
    Ok(HilbertSpace::tensor_all(sites.iter().cloned())?)
}

/// Model constructors for the selected SpinOne system.
///
/// No values are stored; select the coefficient type through `SpinOne<C>`.
pub struct SpinOneModels<Q>(core::marker::PhantomData<fn() -> Q>);

impl<C> SpinOneModels<SpinOne<C>> {
    /// Builds `-sum_(i,j) sum_(a in i,b in j) (Jx Sx_a Sx_b + Jy Sy_a Sy_b + Jz Sz_a Sz_b)`.
    ///
    /// Each input space is one site, with any number of leaves. All leaf pairs
    /// across nearest-neighbor sites interact; no within-site terms are added.
    /// Empty sites retain their positions without reconnecting neighbors.
    /// Each undirected chain edge occurs once, including for two periodic sites.
    ///
    /// Site count and index uniqueness are checked even without edges or with
    /// zero couplings. The result is bound to a tensor sharing the supplied
    /// site trees in order. Zero sites give zero on the unit space.
    /// Parameters must represent real values for a Hermitian Hamiltonian.
    /// Construction does not evaluate local rules or perform scalar arithmetic
    /// beyond converting and negating the couplings.
    pub fn xyz<I: Eq, R>(
        chain: &Chain,
        sites: &[HilbertSpace<I>],
        jx: R,
        jy: R,
        jz: R,
    ) -> Result<OperatorOf<SpinOne<C>, I>, UtilityError>
    where
        C: Clone + From<R> + Neg<Output = C>,
    {
        let space = model_space(chain, sites)?;
        let axes = [
            (O::Sx, -C::from(jx)),
            (O::Sy, -C::from(jy)),
            (O::Sz, -C::from(jz)),
        ];
        let mut terms = Vec::new();
        let mapper = SiteMapper::new(chain, sites)?;
        for group in mapper
            .groups(&NearestNeighbors)
            .unwrap_or_else(|never| match never {})
        {
            let group = group?;
            for pair in group.leaves_cartesian() {
                let (left, right) = (pair[0], pair[1]);
                for (axis, coupling) in &axes {
                    terms.push(
                        SpinOne::<C>::operator_product_at([(left, *axis), (right, *axis)])?
                            .scaled(coupling.clone()),
                    );
                }
            }
        }
        Ok(SpinOne::<C>::operator_sum(terms)?.bind(&space)?)
    }

    /// Builds `-J sum_(i,j) sum_(a in i,b in j) S_a dot S_b`.
    ///
    /// This is [`Self::xyz`] with all three couplings equal; positive J is
    /// ferromagnetic. Validation, all-leaf pairing and space binding are the same.
    pub fn heisenberg<I: Eq, R: Clone>(
        chain: &Chain,
        sites: &[HilbertSpace<I>],
        j: R,
    ) -> Result<OperatorOf<SpinOne<C>, I>, UtilityError>
    where
        C: Clone + From<R> + Neg<Output = C>,
    {
        Self::xyz(chain, sites, j.clone(), j.clone(), j)
    }

    /// Builds `-sum_(i,j) sum_(a in i,b in j) [J A_ab + K A_ab^2]`, where `A_ab = S_a dot S_b`.
    ///
    /// The square is taken separately for each leaf pair, with all Cartesian
    /// cross terms retained. No constant energy shift is added or discarded.
    /// At K = 0 this equals [`Self::heisenberg`]; J > 0 is ferromagnetic.
    /// The angular convention `sum [cos(theta) A + sin(theta) A^2]` corresponds
    /// to `j = -cos(theta), k = -sin(theta)`.
    ///
    /// Site validation, boundaries, empty sites and binding follow [`Self::xyz`].
    /// Couplings must represent real values for a Hermitian Hamiltonian.
    /// The unscaled bilinear expression is shared in its delayed square, so
    /// construction requires no coefficient multiplication, fractions or Float.
    ///
    /// ```
    /// use num_complex::Complex;
    /// use qmb_core::HilbertSpace;
    /// use qmb_geometry::{Boundary, Chain};
    /// use qmb_spin_one::{SpinOne, SpinOneState as Ket};
    /// type Spin = SpinOne<Complex<f64>>;
    /// let sites = [HilbertSpace::factor("a"), HilbertSpace::factor("b")];
    /// let h = qmb_spin_one::models::SpinOneModels::<Spin>::blbq(&Chain::new(2, Boundary::Open), &sites, 1.0, 0.5)?;
    /// let psi = Spin::product_state(h.space().unwrap(), [Ket::Plus, Ket::Plus])?;
    /// assert!((psi.expectation(&h)?.re + 1.5).abs() < 1e-12);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn blbq<I: Eq, R>(
        chain: &Chain,
        sites: &[HilbertSpace<I>],
        j: R,
        k: R,
    ) -> Result<OperatorOf<SpinOne<C>, I>, UtilityError>
    where
        C: Clone + From<R> + Neg<Output = C>,
    {
        let space = model_space(chain, sites)?;
        let minus_j = -C::from(j);
        let minus_k = -C::from(k);
        let mut terms = Vec::new();
        let mapper = SiteMapper::new(chain, sites)?;
        for group in mapper
            .groups(&NearestNeighbors)
            .unwrap_or_else(|never| match never {})
        {
            let group = group?;
            for pair in group.leaves_cartesian() {
                let (left, right) = (pair[0], pair[1]);
                let bilinear =
                    SpinOne::<C>::try_operator_sum([O::Sx, O::Sy, O::Sz].map(|axis| {
                        SpinOne::<C>::operator_product_at([(left, axis), (right, axis)])
                    }))?;
                terms.push(bilinear.scaled(minus_j.clone()));
                terms.push((bilinear.clone() * &bilinear)?.scaled(minus_k.clone()));
            }
        }
        Ok(SpinOne::<C>::operator_sum(terms)?.bind(&space)?)
    }
}
