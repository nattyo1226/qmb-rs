//! Expansion preserves operation-specific capabilities and shared factor indices.

use num_traits::{One, Zero};
use qmb_core::{
    Act, Adjoint, Basis, BasisError, Conjugate, Dimension, DimensionOverflow, HilbertSpace,
    LinearCombination, LocalState, State,
};
use std::{
    cell::Cell,
    num::NonZeroUsize,
    ops::{Add, Mul},
};

// No Clone, Eq, ordering or hashing on either primitive. In particular, the
// operator implements Adjoint but deliberately has no Product or Identity.
struct Ket(bool);
struct Action;
impl qmb_core::HomogeneousParity for Action {
    fn parity(&self) -> qmb_core::Parity {
        qmb_core::Parity::Even
    }
}
impl qmb_core::ParityOperator for Action {
    fn parity_operator() -> Self {
        Self
    }
}
impl Dimension for Ket {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(2).unwrap())
    }
}
impl<C: From<i16>> Act<C> for Ket {
    type Operator = Action;
    type Term = Self;
    fn act(&self, _: &Action) -> LinearCombination<Self, C> {
        LinearCombination::term(1i16.into(), Ket(self.0))
    }
}
impl<C: From<i16>> Adjoint<C> for Action {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, C> {
        LinearCombination::term(1i16.into(), Action)
    }
}

struct CountingBasis {
    expansions: Cell<usize>,
    enumerations: Cell<usize>,
}
impl CountingBasis {
    fn new() -> Self {
        Self {
            expansions: Cell::new(0),
            enumerations: Cell::new(0),
        }
    }
}
impl<C: From<i16>> Basis<Ket, C> for CountingBasis {
    fn expand(&self, p: &Ket) -> LinearCombination<Ket, C> {
        self.expansions.set(self.expansions.get() + 1);
        LinearCombination::term(1i16.into(), Ket(p.0))
    }
    fn index(&self, p: &Ket) -> Option<usize> {
        Some(usize::from(p.0))
    }
    fn elements(&self) -> impl Iterator<Item = Ket> + '_ {
        self.enumerations.set(self.enumerations.get() + 1);
        [Ket(false), Ket(true)].into_iter()
    }
}

// No One, Float, Num, division, subtraction, negation, or approximate comparison.
#[derive(Clone, Debug, PartialEq)]
struct Coefficient(i32);
impl From<i16> for Coefficient {
    fn from(n: i16) -> Self {
        Self(i32::from(n))
    }
}
impl Add for Coefficient {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}
impl Mul for Coefficient {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0)
    }
}
impl Zero for Coefficient {
    fn zero() -> Self {
        Self(0)
    }
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}
impl Conjugate for Coefficient {
    type Output = Self;
    fn conjugate(&self) -> Self {
        self.clone()
    }
}

#[test]
fn local_expansion_needs_no_primitive_clone_product_float_or_coefficient_one() {
    type S = LocalState<Ket, Action, Coefficient>;
    let basis = CountingBasis::new();
    let state = S::primitive(Ket(true));
    let delayed = (state.clone() + state.scaled(Coefficient(2)))
        .act(&qmb_core::LocalOperator::primitive(Action));
    let expanded = delayed.expand_in(&basis).unwrap();
    let coords = expanded.coordinates_in(&basis).unwrap();
    assert_eq!(coords.len(), 1);
    assert_eq!(
        coords.entries().next().unwrap(),
        (&[1][..], &Coefficient(3))
    );
    assert_eq!(
        state
            .coordinates_in(&basis)
            .unwrap()
            .entries()
            .next()
            .unwrap()
            .1,
        &Coefficient(1)
    );
    assert_eq!(basis.enumerations.get(), 0);
}

// Used only by sparse expansion regressions to count coefficient products.
#[derive(Clone, Debug, PartialEq)]
struct Counted(f64);
thread_local! {
    static MULTIPLICATIONS: Cell<usize> = const { Cell::new(0) };
}

fn reset_multiplications() {
    MULTIPLICATIONS.with(|count| count.set(0));
}
fn multiplication_count() -> usize {
    MULTIPLICATIONS.with(Cell::get)
}
fn count_multiplication() {
    MULTIPLICATIONS.with(|count| count.set(count.get() + 1));
}

impl From<i16> for Counted {
    fn from(n: i16) -> Self {
        Self(f64::from(n))
    }
}
impl Add for Counted {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}
impl Mul for Counted {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        count_multiplication();
        Self(self.0 * rhs.0)
    }
}
impl Zero for Counted {
    fn zero() -> Self {
        Self(0.0)
    }
    fn is_zero(&self) -> bool {
        self.0 == 0.0
    }
}
impl One for Counted {
    fn one() -> Self {
        Self(1.0)
    }
}
impl Conjugate for Counted {
    type Output = Self;
    fn conjugate(&self) -> Self {
        self.clone()
    }
}

#[test]
fn cancellations_stay_sparse_through_tensor_expansion() {
    const FACTORS: usize = 12;
    type S = LocalState<Ket, Action, Counted>;
    let basis = CountingBasis::new();
    let space = HilbertSpace::tensor((0..FACTORS).map(HilbertSpace::factor)).unwrap();
    let state = State::tensor_in(
        &space,
        (0..FACTORS).map(|index| {
            let factor = HilbertSpace::factor(index);
            let zero = S::primitive(Ket(false)).at(&factor).unwrap();
            let one = S::primitive(Ket(true)).at(&factor).unwrap();
            State::sum(&factor, [zero.clone(), zero.scaled(Counted(-1.0)), one]).unwrap()
        }),
    )
    .unwrap();

    // A generous bound catches exponential work without relying on wall time
    // or a particular linear-time implementation.
    reset_multiplications();
    let expanded = state.expand_in(&basis).unwrap();
    assert!(multiplication_count() < 1_000);
    assert!(expanded.space().shares_node_with(state.space()));

    reset_multiplications();
    let coordinates = state.coordinates_in(&basis).unwrap();
    assert!(multiplication_count() < 1_000);
    assert_eq!(coordinates.len(), 1);
    let (digits, coefficient) = coordinates.entries().next().unwrap();
    assert_eq!(digits.len(), FACTORS);
    assert!(digits.iter().all(|digit| *digit == 1));
    assert_eq!(*coefficient, Counted(1.0));
    let expanded_coordinates = expanded.coordinates_in(&basis).unwrap();
    assert_eq!(
        expanded_coordinates.entries().collect::<Vec<_>>(),
        coordinates.entries().collect::<Vec<_>>()
    );

    let factor = HilbertSpace::factor(0usize);
    let zero = S::primitive(Ket(false)).at(&factor).unwrap();
    let one = S::primitive(Ket(true)).at(&factor).unwrap();
    let mut scaled = State::sum(&factor, [zero, one])
        .unwrap()
        .scaled(Counted(0.0));
    for _ in 0..100 {
        scaled = scaled.scaled(Counted(1.0));
    }
    reset_multiplications();
    let coordinates = scaled.coordinates_in(&basis).unwrap();
    assert!(multiplication_count() < 100);
    assert!(coordinates.is_empty());
    assert!(coordinates.space().shares_node_with(scaled.space()));
}

#[test]
fn underflowed_tensor_coefficients_do_not_keep_expanding() {
    const FACTORS: usize = 12;
    type S = LocalState<Ket, Action, Counted>;
    let basis = CountingBasis::new();
    let space = HilbertSpace::tensor((0..FACTORS).map(HilbertSpace::factor)).unwrap();
    let state = State::tensor_in(
        &space,
        (0..FACTORS).map(|index| {
            let factor = HilbertSpace::factor(index);
            let small = S::primitive(Ket(false)).at(&factor).unwrap();
            let one = S::primitive(Ket(true)).at(&factor).unwrap();
            State::sum(&factor, [small.scaled(Counted(1e-200)), one]).unwrap()
        }),
    )
    .unwrap();
    reset_multiplications();
    let coordinates = state.coordinates_in(&basis).unwrap();
    assert!(multiplication_count() < 1_000);
    assert_eq!(coordinates.len(), FACTORS + 1);
    for (digits, coefficient) in coordinates.entries() {
        let small_factors = digits.iter().filter(|digit| **digit == 0).count();
        assert!(small_factors <= 1);
        assert_eq!(
            *coefficient,
            Counted(if small_factors == 0 { 1.0 } else { 1e-200 })
        );
    }
}

#[test]
fn composite_expansion_retains_tiny_and_nan_coefficients() {
    type S = LocalState<Ket, Action, f64>;
    let factor = HilbertSpace::factor(0u8);
    let zero = S::primitive(Ket(false)).at(&factor).unwrap();
    let one = S::primitive(Ket(true)).at(&factor).unwrap();
    let zero_branch = State::sum(&factor, [zero.clone(), one.clone()])
        .unwrap()
        .scaled(0.0);
    let state = State::sum(
        &factor,
        [zero_branch, zero.scaled(1e-300), one.scaled(f64::NAN)],
    )
    .unwrap();
    let basis = CountingBasis::new();
    let expanded = state.expand_in(&basis).unwrap();
    for state in [&state, &expanded] {
        let coordinates = state.coordinates_in(&basis).unwrap();
        let entries: Vec<_> = coordinates.entries().collect();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0], (&[0][..], &1e-300));
        assert_eq!(entries[1].0, &[1]);
        assert!(entries[1].1.is_nan());
    }
}

#[test]
fn basis_enumeration_needs_no_index_comparison() {
    struct NonComparableIndex;
    let basis = CountingBasis::new();
    let single = HilbertSpace::factor(NonComparableIndex);
    let mut states = single.basis_states::<Ket, Action, f64, _>(&basis).unwrap();
    for _ in 0..2 {
        assert!(
            states
                .next()
                .unwrap()
                .unwrap()
                .space()
                .shares_node_with(&single)
        );
    }
    assert!(states.next().is_none());
    let mut operators = single.basis_operators::<Ket, f64, _>(&basis).unwrap();
    for _ in 0..2 {
        assert!(
            operators
                .next()
                .unwrap()
                .unwrap()
                .space()
                .unwrap()
                .shares_node_with(&single)
        );
    }
    assert!(operators.next().is_none());

    let unit = HilbertSpace::<NonComparableIndex>::unit();
    assert_eq!(
        unit.basis_states::<Ket, Action, f64, _>(&basis)
            .unwrap()
            .count(),
        1
    );
    assert_eq!(
        unit.basis_operators::<Ket, f64, _>(&basis).unwrap().count(),
        1
    );
}

#[derive(PartialEq, Eq)]
struct Index(u8); // No Clone, Hash, Ord or Debug.

#[test]
fn placed_expansion_shares_nonclone_indices_and_does_not_enumerate_all_basis_elements() {
    type S = LocalState<Ket, Action, f64>;
    let basis = CountingBasis::new();
    let a = HilbertSpace::factor(Index(11));
    let b = HilbertSpace::factor(Index(7));
    let state = State::tensor_in(
        &a.tensor_product(&b).unwrap(),
        [
            S::primitive(Ket(true)).at(&a).unwrap(),
            S::primitive(Ket(false)).at(&b).unwrap(),
        ],
    )
    .unwrap();
    let op = qmb_core::LocalOperator::primitive(Action).at(&b).unwrap();
    let expanded = state.act(&op).unwrap().expand_in(&basis).unwrap();
    assert!(expanded.space().shares_node_with(state.space()));
    let coords = expanded.coordinates_in(&basis).unwrap();
    assert!(coords.space().shares_node_with(state.space()));
    assert_eq!(coords.flat_entries().next().unwrap().unwrap(), (2, &1.0));
    assert_eq!(basis.enumerations.get(), 0);
    assert_eq!(
        state
            .space()
            .basis_states::<Ket, Action, f64, _>(&basis)
            .unwrap()
            .take(1)
            .count(),
        1
    );
}

struct Large;
impl Dimension for Large {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Err(DimensionOverflow)
    }
}
impl Act<f64> for Large {
    type Operator = Action;
    type Term = Self;
    fn act(&self, _: &Action) -> LinearCombination<Self, f64> {
        LinearCombination::term(1.0, Large)
    }
}
struct LargeBasis;
impl Basis<Large, f64> for LargeBasis {
    fn expand(&self, _: &Large) -> LinearCombination<Large, f64> {
        LinearCombination::term(1.0, Large)
    }
    fn index(&self, _: &Large) -> Option<usize> {
        Some(0)
    }
    fn elements(&self) -> impl Iterator<Item = Large> {
        std::iter::once(Large)
    }
}

#[test]
fn unit_space_does_not_request_an_unused_local_dimension() {
    let state = State::<u8, Large, Action, f64>::unit();
    let expanded = state.expand_in(&LargeBasis).unwrap();
    let coords = expanded.coordinates_in(&LargeBasis).unwrap();
    assert_eq!(coords.dimension().unwrap(), 1);
    assert_eq!(coords.flat_entries().next().unwrap().unwrap(), (0, &1.0));
    assert_eq!(
        state
            .space()
            .basis_states::<Large, Action, f64, _>(&LargeBasis)
            .unwrap()
            .count(),
        1
    );
    let factor = HilbertSpace::factor(0u8);
    assert!(matches!(
        LocalState::<Large, Action, f64>::primitive(Large)
            .at(&factor)
            .unwrap()
            .expand_in(&LargeBasis),
        Err(BasisError::Dimension(DimensionOverflow))
    ));
}
