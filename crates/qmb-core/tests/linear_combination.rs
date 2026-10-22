use std::{cell::Cell, iter::FusedIterator, ops, rc::Rc};

use qmb_core::{
    Complex, LinearCombination, Term,
    num_traits::{Num, One, Zero},
};

fn term(primitive: u8, re: f64, im: f64) -> Term<u8, f64> {
    Term {
        coefficient: Complex::new(re, im),
        operand: primitive,
    }
}

#[test]
fn construction_keeps_order_duplicates_and_zero_coefficients() {
    for n in 0..=4 {
        let terms: Vec<_> = (0..n).map(|i| term(i % 2, f64::from(i), -1.0)).collect();
        let sum = LinearCombination::from_terms(terms.clone());
        assert_eq!(sum.len(), usize::from(n));
        assert_eq!(sum.is_empty(), n == 0);
        assert_eq!(sum.terms().cloned().collect::<Vec<_>>(), terms);
        assert_eq!(sum.clone().into_terms().collect::<Vec<_>>(), terms);
        assert_eq!(sum, terms.into_iter().collect());
    }
    let primitive = LinearCombination::<_, f64>::primitive(7);
    assert_eq!(
        primitive,
        LinearCombination::term(Complex::new(1.0, 0.0), 7)
    );
    assert_eq!(
        LinearCombination::<u8, f64>::default(),
        LinearCombination::zero()
    );
    let zero_term = LinearCombination::term(Complex::new(0.0, 0.0), 0);
    assert!(!zero_term.is_empty());
    assert_ne!(zero_term, LinearCombination::zero());
}

#[test]
fn map_moves_opaque_coefficients_and_operands_once_in_order() {
    struct Scalar(u8);
    struct Operand(u8);
    let original = LinearCombination::from_terms((0..4).map(|i| Term {
        coefficient: Complex::new(Scalar(i), Scalar(7 - i)),
        operand: Operand(i % 2),
    }));
    let mut seen = Vec::new();
    let mapped = original.map(|operand| {
        seen.push(operand.0);
        operand.0 + 10
    });
    assert_eq!(seen, vec![0, 1, 0, 1]);
    for (i, term) in mapped.into_iter().enumerate() {
        assert_eq!(term.coefficient.re.0, i as u8);
        assert_eq!(term.coefficient.im.0, 7 - i as u8);
        assert_eq!(term.operand, 10 + (i as u8) % 2);
    }
    let zero = LinearCombination::<Operand, Scalar>::zero().map(|_| panic!("empty map"));
    assert!(zero.is_empty());
}

#[test]
fn sum_parts_support_borrowed_and_owned_mixed_end_iteration() {
    use qmb_core::expr::parts::SumExpr;
    for n in 0..5 {
        let sum: SumExpr<_> = (0..n).collect();
        assert_eq!(sum.len(), n);
        assert_eq!(
            sum.iter().copied().collect::<Vec<_>>(),
            (0..n).collect::<Vec<_>>()
        );
        let mut iter = sum.into_iter();
        if n > 0 {
            assert_eq!(iter.next_back(), Some(n - 1));
        }
        if n > 1 {
            assert_eq!(iter.next(), Some(0));
        }
        assert_eq!(
            iter.collect::<Vec<_>>(),
            (1..n.saturating_sub(1)).collect::<Vec<_>>()
        );
    }
}

fn check_iterator<I>(mut iter: I, expected: &[u8])
where
    I: DoubleEndedIterator<Item = u8> + ExactSizeIterator + FusedIterator,
{
    let mut front = 0;
    let mut back = expected.len();
    while front < back {
        assert_eq!(iter.len(), back - front);
        assert_eq!(iter.size_hint(), (back - front, Some(back - front)));
        assert_eq!(iter.next(), Some(expected[front]));
        front += 1;
        if front < back {
            back -= 1;
            assert_eq!(iter.next_back(), Some(expected[back]));
        }
    }
    assert_eq!(iter.len(), 0);
    assert_eq!(iter.size_hint(), (0, Some(0)));
    for _ in 0..3 {
        assert_eq!(iter.next(), None);
        assert_eq!(iter.next_back(), None);
    }
}

#[test]
fn borrowed_and_owned_iterators_support_mixed_ends_and_exhaustion() {
    for n in 0..=6 {
        let expected: Vec<_> = (0..n).collect();
        let sum =
            LinearCombination::<_, f64>::from_terms(expected.iter().map(|&p| term(p, 0.0, 0.0)));
        check_iterator(sum.terms().map(|t| t.operand), &expected);
        check_iterator((&sum).into_iter().map(|t| t.operand), &expected);
        check_iterator(sum.clone().into_terms().map(|t| t.operand), &expected);
        check_iterator(sum.into_iter().map(|t| t.operand), &expected);
    }
}

#[test]
fn concatenation_and_owned_operations_require_no_scalar_or_primitive_traits() {
    struct Primitive(u8);
    struct Scalar;
    fn combination(start: u8, len: u8) -> LinearCombination<Primitive, Scalar> {
        LinearCombination::from_terms((start..start + len).map(|p| Term {
            coefficient: Complex::new(Scalar, Scalar),
            operand: Primitive(p),
        }))
    }
    let _: LinearCombination<Primitive, Scalar> = LinearCombination::default();
    let one = LinearCombination::term(Complex::new(Scalar, Scalar), Primitive(9));
    assert_eq!(one.into_terms().next().unwrap().operand.0, 9);
    // Covers all nine Zero / Term / Sum pairs with non-Clone, non-Debug,
    // non-PartialEq primitives and scalars that have no numeric capabilities.
    for left in 0..=3 {
        for right in 0..=3 {
            let sum = combination(0, left) + combination(10, right);
            let expected: Vec<_> = (0..left).chain(10..10 + right).collect();
            assert_eq!(sum.len(), expected.len());
            assert_eq!(
                sum.into_iter().map(|t| t.operand.0).collect::<Vec<_>>(),
                expected
            );
        }
    }
}

#[test]
fn equality_is_structural_and_follows_scalar_partial_eq() {
    let a = term(0, 1.0, 0.0);
    let b = term(1, 2.0, 0.0);
    let sum = LinearCombination::from_terms([a.clone(), b.clone()]);
    assert_eq!(sum, sum.clone());
    assert_ne!(sum, LinearCombination::from_terms([b, a.clone()]));
    assert_ne!(
        LinearCombination::from_terms([a.clone(), a]),
        LinearCombination::term(Complex::new(2.0, 0.0), 0)
    );
    let nan = LinearCombination::term(Complex::new(f64::NAN, 0.0), 0);
    assert_ne!(nan, nan.clone());
    assert!(!format!("{sum:?}").is_empty());
}

#[test]
fn coefficient_one_only_needs_zero_and_one_without_clone_or_num() {
    struct MinimalScalar(u8);
    impl ops::Add for MinimalScalar {
        type Output = Self;
        fn add(self, rhs: Self) -> Self {
            Self(self.0 + rhs.0)
        }
    }
    impl ops::Mul for MinimalScalar {
        type Output = Self;
        fn mul(self, rhs: Self) -> Self {
            Self(self.0 * rhs.0)
        }
    }
    impl Zero for MinimalScalar {
        fn zero() -> Self {
            Self(0)
        }
        fn is_zero(&self) -> bool {
            self.0 == 0
        }
    }
    impl One for MinimalScalar {
        fn one() -> Self {
            Self(1)
        }
    }
    struct Primitive;
    let sum = LinearCombination::<_, MinimalScalar>::primitive(Primitive);
    let term = sum.into_terms().next().unwrap();
    assert_eq!(term.coefficient.re.0, 1);
    assert_eq!(term.coefficient.im.0, 0);
}

#[derive(Clone, Debug, PartialEq)]
struct HeapReal(Box<f64>);

impl HeapReal {
    fn new(value: f64) -> Self {
        Self(Box::new(value))
    }
}

macro_rules! arithmetic {
    ($trait:ident, $method:ident, $op:tt) => {
        impl ops::$trait for HeapReal {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self { Self::new(*self.0 $op *rhs.0) }
        }
    };
}
arithmetic!(Add, add, +);
arithmetic!(Sub, sub, -);
arithmetic!(Mul, mul, *);
arithmetic!(Div, div, /);
arithmetic!(Rem, rem, %);

impl Zero for HeapReal {
    fn zero() -> Self {
        Self::new(0.0)
    }
    fn is_zero(&self) -> bool {
        *self.0 == 0.0
    }
}
impl One for HeapReal {
    fn one() -> Self {
        Self::new(1.0)
    }
}
impl Num for HeapReal {
    type FromStrRadixErr = <f64 as Num>::FromStrRadixErr;
    fn from_str_radix(s: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        f64::from_str_radix(s, radix).map(Self::new)
    }
}

#[test]
fn scaling_needs_neither_copy_scalars_nor_clone_primitives() {
    struct Primitive(u8);
    for n in 0..=3 {
        let sum = LinearCombination::from_terms((0..n).map(|p| Term {
            coefficient: Complex::new(HeapReal::new(f64::from(p)), HeapReal::new(2.0)),
            operand: Primitive(p),
        }));
        let scaled = sum.scaled(Complex::new(HeapReal::new(3.0), HeapReal::new(-1.0)));
        assert_eq!(scaled.len(), usize::from(n));
        for (p, t) in scaled.into_iter().enumerate() {
            assert_eq!(usize::from(t.operand.0), p);
            assert_eq!(*t.coefficient.re.0, f64::from(t.operand.0) * 3.0 + 2.0);
            assert_eq!(*t.coefficient.im.0, 6.0 - f64::from(t.operand.0));
        }
    }
    let primitive = LinearCombination::<_, HeapReal>::primitive(Primitive(8));
    assert_eq!(*primitive.terms().next().unwrap().coefficient.re.0, 1.0);
    let scaled = primitive.scaled(Complex::new(HeapReal::zero(), HeapReal::zero()));
    assert_eq!(scaled.len(), 1);
    assert_eq!(
        scaled.terms().next().unwrap().coefficient,
        Complex::new(HeapReal::zero(), HeapReal::zero())
    );
}

#[test]
fn borrowed_non_send_primitives_and_partial_iterator_drop_are_supported() {
    struct Primitive<'a>(&'a Rc<Cell<usize>>);
    impl Drop for Primitive<'_> {
        fn drop(&mut self) {
            self.0.set(self.0.get() + 1);
        }
    }
    let drops = Rc::new(Cell::new(0));
    for n in 0..=4 {
        drops.set(0);
        let sum = LinearCombination::<_, f64>::from_terms((0..n).map(|_| Term {
            coefficient: Complex::new(0.0, 0.0),
            operand: Primitive(&drops),
        }));
        let mut iter = sum.into_terms();
        drop(iter.next());
        drop(iter.next_back());
        drop(iter);
        assert_eq!(drops.get(), n);
    }
}
