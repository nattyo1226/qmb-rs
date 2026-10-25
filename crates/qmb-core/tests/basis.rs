//! The unified basis trait does not prescribe numeric capabilities.

use qmb_core::{Basis, Complex, LinearCombination};

enum Primitive {
    First,
    Second,
    Sum,
}
struct Scalar(u8);
struct CustomBasis;

impl Basis<Primitive, Scalar> for CustomBasis {
    fn expand(&self, primitive: &Primitive) -> LinearCombination<Primitive, Scalar> {
        let term = |operand| LinearCombination::term(Complex::new(Scalar(1), Scalar(0)), operand);
        match primitive {
            Primitive::First => term(Primitive::First),
            Primitive::Second => term(Primitive::Second),
            Primitive::Sum => term(Primitive::First) + term(Primitive::Second),
        }
    }

    fn index(&self, primitive: &Primitive) -> Option<usize> {
        match primitive {
            Primitive::First => Some(0),
            Primitive::Second => Some(1),
            Primitive::Sum => None,
        }
    }

    fn elements(&self) -> impl Iterator<Item = Primitive> + '_ {
        [Primitive::First, Primitive::Second].into_iter()
    }
}

#[test]
fn one_basis_implementation_supports_all_methods_without_scalar_traits() {
    let basis = CustomBasis;
    let elements: Vec<_> = basis.elements().collect();
    assert_eq!(elements.len(), 2);
    assert_eq!(basis.index(&elements[0]), Some(0));
    assert_eq!(basis.index(&elements[1]), Some(1));
    assert_eq!(basis.index(&Primitive::Sum), None);
    let expansion = basis.expand(&Primitive::Sum);
    assert_eq!(expansion.len(), 2);
    for (i, term) in expansion.into_iter().enumerate() {
        assert_eq!(basis.index(&term.operand), Some(i));
        assert_eq!((term.coefficient.re.0, term.coefficient.im.0), (1, 0));
    }
}
