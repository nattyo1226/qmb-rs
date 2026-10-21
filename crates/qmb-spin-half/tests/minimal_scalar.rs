//! Capability tests deliberately omit Clone, Copy, Num, Debug, and PartialEq.

use std::ops::{Add, Div, Mul, Neg};

use qmb_core::{
    Act, Adjoint, Basis, InnerProduct, LinearCombination as Lc, LocalStructure, OperatorBasis,
    Product, StateBasis, Trace,
    num_traits::{One, Zero},
};
use qmb_spin_half::{
    CartesianBasis, LadderBasis, SpinHalf, SpinHalfOperator as Op, SpinHalfState as State, SzBasis,
};

use Op::{I, SMinus, SPlus, Sx, Sy, Sz};
use State::{Down, Up};

const OPS: [Op; 6] = [I, Sx, Sy, Sz, SPlus, SMinus];

macro_rules! integer_scalar {
    ($name:ident, $integer:ty) => {
        struct $name($integer);
        impl Add for $name {
            type Output = Self;
            fn add(self, rhs: Self) -> Self {
                Self(self.0 + rhs.0)
            }
        }
        impl Mul for $name {
            type Output = Self;
            fn mul(self, rhs: Self) -> Self {
                Self(self.0 * rhs.0)
            }
        }
        impl Zero for $name {
            fn zero() -> Self {
                Self(0)
            }
            fn is_zero(&self) -> bool {
                self.0 == 0
            }
        }
        impl One for $name {
            fn one() -> Self {
                Self(1)
            }
        }
    };
}
integer_scalar!(Natural, u32);
integer_scalar!(Signed, i64);

impl Neg for Signed {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

#[test]
fn inner_adjoint_trace_and_state_basis_only_need_zero_and_one() {
    for a in [Up, Down] {
        for b in [Up, Down] {
            let inner = <State as InnerProduct<Natural>>::inner_product(&a, &b);
            assert_eq!(inner.re.0, if a == b { 1 } else { 0 });
            assert_eq!(inner.im.0, 0);
        }
    }
    for (op, expected) in OPS.into_iter().zip([I, Sx, Sy, Sz, SMinus, SPlus]) {
        let adjoint = <Op as Adjoint<Natural>>::adjoint(&op);
        assert_eq!(adjoint.len(), 1);
        let term = adjoint.into_iter().next().unwrap();
        assert_eq!(term.primitive, expected);
        assert_eq!(term.coefficient.re.0, 1);
        assert_eq!(term.coefficient.im.0, 0);
        let trace = <Op as Trace<Natural>>::trace(&op);
        assert_eq!(trace.re.0, if op == I { 2 } else { 0 });
        assert_eq!(trace.im.0, 0);
    }
    let elements: Vec<_> = <SzBasis as Basis<State, Natural>>::elements(&SzBasis).collect();
    assert_eq!(elements, [Up, Down]);
    for (i, state) in elements.into_iter().enumerate() {
        assert_eq!(
            <SzBasis as Basis<State, Natural>>::index(&SzBasis, &state),
            Some(i)
        );
        let expanded = <SzBasis as Basis<State, Natural>>::expand(&SzBasis, &state);
        let mut terms = expanded.into_terms();
        let term = terms.next().unwrap();
        assert_eq!(term.primitive, state);
        assert_eq!(term.coefficient.re.0, 1);
        assert_eq!(term.coefficient.im.0, 0);
        assert!(terms.next().is_none());
    }
}

#[test]
fn cartesian_basis_needs_negation_but_not_division() {
    let expected: [&[(i64, i64, Op)]; 6] = [
        &[(1, 0, I)],
        &[(1, 0, Sx)],
        &[(1, 0, Sy)],
        &[(1, 0, Sz)],
        &[(1, 0, Sx), (0, 1, Sy)],
        &[(1, 0, Sx), (0, -1, Sy)],
    ];
    let elements: Vec<_> =
        <CartesianBasis as Basis<Op, Signed>>::elements(&CartesianBasis).collect();
    assert_eq!(elements, [I, Sx, Sy, Sz]);
    for (op, expected) in OPS.into_iter().zip(expected) {
        let expanded = <CartesianBasis as Basis<Op, Signed>>::expand(&CartesianBasis, &op);
        let actual: Vec<_> = expanded
            .into_terms()
            .map(|t| (t.coefficient.re.0, t.coefficient.im.0, t.primitive))
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(
            <CartesianBasis as Basis<Op, Signed>>::index(&CartesianBasis, &op),
            elements.iter().position(|&p| p == op)
        );
    }
}

// Exact fractions with only the promised basic operations, and no cloning or
// parsing implementation. Values are small throughout these capability tests.
struct Rational {
    numerator: i64,
    denominator: i64,
}

impl Rational {
    fn new(mut numerator: i64, mut denominator: i64) -> Self {
        assert_ne!(denominator, 0);
        if denominator < 0 {
            numerator = -numerator;
            denominator = -denominator;
        }
        let mut a = numerator.abs();
        let mut b = denominator;
        while b != 0 {
            (a, b) = (b, a % b);
        }
        Self {
            numerator: numerator / a,
            denominator: denominator / a,
        }
    }
    fn parts(self) -> (i64, i64) {
        (self.numerator, self.denominator)
    }
}
impl Add for Rational {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(
            self.numerator * rhs.denominator + rhs.numerator * self.denominator,
            self.denominator * rhs.denominator,
        )
    }
}
impl Mul for Rational {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.numerator * rhs.numerator,
            self.denominator * rhs.denominator,
        )
    }
}
impl Div for Rational {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        Self::new(
            self.numerator * rhs.denominator,
            self.denominator * rhs.numerator,
        )
    }
}
impl Neg for Rational {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            numerator: -self.numerator,
            denominator: self.denominator,
        }
    }
}
impl Zero for Rational {
    fn zero() -> Self {
        Self::new(0, 1)
    }
    fn is_zero(&self) -> bool {
        self.numerator == 0
    }
}
impl One for Rational {
    fn one() -> Self {
        Self::new(1, 1)
    }
}

type Fraction = (i64, i64);
type RecordedTerm<P> = (Fraction, Fraction, P);

fn terms<P>(sum: Lc<P, Rational>) -> Vec<RecordedTerm<P>> {
    sum.into_terms()
        .map(|t| {
            (
                t.coefficient.re.parts(),
                t.coefficient.im.parts(),
                t.primitive,
            )
        })
        .collect()
}

#[test]
fn products_and_actions_use_exact_fractions_without_clone_copy_or_num() {
    // Compile and exercise every branch with a non-Clone numeric type.
    for a in OPS {
        for b in OPS {
            let product: Lc<Op, Rational> = a.product(&b);
            assert!(product.len() <= 2);
            for (re, im, _) in terms(product) {
                assert!([1, 2, 4].contains(&re.1));
                assert!([1, 2, 4].contains(&im.1));
            }
        }
        for state in [Up, Down] {
            let action: Lc<State, Rational> = state.act(&a);
            assert!(action.len() <= 1);
            for (re, im, _) in terms(action) {
                assert!([1, 2].contains(&re.1));
                assert!([1, 2].contains(&im.1));
            }
        }
    }
    assert_eq!(
        terms(<Op as Product<Rational>>::product(&Sy, &SMinus)),
        [((0, 1), (-1, 4), I), ((0, 1), (-1, 2), Sz),]
    );
    assert_eq!(
        terms(<Op as Product<Rational>>::product(&SPlus, &SMinus)),
        [((1, 2), (0, 1), I), ((1, 1), (0, 1), Sz),]
    );
    assert_eq!(
        terms(<Op as Product<Rational>>::product(&SMinus, &SPlus)),
        [((1, 2), (0, 1), I), ((-1, 1), (0, 1), Sz),]
    );
    assert_eq!(
        terms(<State as Act<Op, Rational>>::act(&Down, &Sy)),
        [((0, 1), (-1, 2), Up)]
    );
    assert!(<Op as Product<Rational>>::product(&SPlus, &SPlus).is_empty());
    assert!(<State as Act<Op, Rational>>::act(&Up, &SPlus).is_empty());
}

#[test]
fn local_structure_and_all_basis_markers_accept_minimal_fractional_scalars() {
    fn state_basis<B: StateBasis<SpinHalf, Rational>>(_: B) {}
    fn operator_basis<B: OperatorBasis<SpinHalf, Rational>>(_: B) {}
    state_basis(SzBasis);
    operator_basis(CartesianBasis);
    operator_basis(LadderBasis);
    assert_eq!(<SpinHalf as LocalStructure<Rational>>::DIMENSION.get(), 2);
    assert_eq!(
        terms(<SpinHalf as LocalStructure<Rational>>::identity()),
        [((1, 1), (0, 1), I)]
    );
    let trace = <Op as Trace<Rational>>::trace(&I);
    assert_eq!(trace.re.parts(), (2, 1));
    assert_eq!(trace.im.parts(), (0, 1));

    for op in OPS {
        let cartesian = <CartesianBasis as Basis<Op, Rational>>::expand(&CartesianBasis, &op);
        let ladder = <LadderBasis as Basis<Op, Rational>>::expand(&LadderBasis, &op);
        assert!(!cartesian.is_empty());
        assert!(!ladder.is_empty());
    }
    assert_eq!(
        terms(<LadderBasis as Basis<Op, Rational>>::expand(
            &LadderBasis,
            &Sx
        )),
        [((1, 2), (0, 1), SPlus), ((1, 2), (0, 1), SMinus),]
    );
    assert_eq!(
        terms(<LadderBasis as Basis<Op, Rational>>::expand(
            &LadderBasis,
            &Sy
        )),
        [((0, 1), (-1, 2), SPlus), ((0, 1), (1, 2), SMinus),]
    );
    let elements: Vec<_> = <LadderBasis as Basis<Op, Rational>>::elements(&LadderBasis).collect();
    assert_eq!(elements, [I, SPlus, SMinus, Sz]);
}
