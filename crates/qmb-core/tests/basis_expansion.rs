//! Basis coordinates and reconstruction against independent three-level algebra.

use num_complex::Complex;
use qmb_core::{Basis, BasisError, Coordinates, HilbertSpace, LocalOperator, LocalState};

#[path = "support/three_level.rs"]
mod fixture;
use fixture::*;

type Ls<T> = LocalState<Ket, Operator, Complex<T>>;
type Lo<T> = LocalOperator<Operator, Complex<T>>;

fn state_vector<T: Real, B: Basis<Ket, Complex<T>>, H>(
    coordinates: &Coordinates<'_, B, Complex<T>, H>,
) -> Vec<Complex<T>> {
    let basis: Vec<_> = coordinates.basis().elements().collect();
    let mut result = vec![c(0.0, 0.0); 3usize.pow(coordinates.factor_count() as u32)];
    for (digits, coefficient) in coordinates.entries() {
        let mut vector = vec![*coefficient];
        for digit in digits {
            vector = vector
                .into_iter()
                .flat_map(|a| basis[*digit].coordinates::<T>().map(|b| a * b))
                .collect();
        }
        for (a, b) in result.iter_mut().zip(vector) {
            *a = *a + b;
        }
    }
    result
}

fn assert_values<T: Real>(
    actual: impl IntoIterator<Item = Complex<T>>,
    expected: impl IntoIterator<Item = Complex<T>>,
) {
    let actual: Vec<_> = actual.into_iter().collect();
    let expected: Vec<_> = expected.into_iter().collect();
    assert_eq!(actual.len(), expected.len());
    for (a, b) in actual.into_iter().zip(expected) {
        assert_complex(a, b);
    }
}

fn check_three_level<T: Real>(kind: qmb_core::TensorKind) {
    let p = Ls::<T>::primitive(Ket::Phase).scaled(c(0.3, -0.7))
        + Ls::primitive(Ket::Sum).scaled(c(-0.2, 0.6));
    let reference = std::array::from_fn::<_, 3, _>(|i| {
        Ket::Phase.coordinates::<T>()[i] * c(0.3, -0.7)
            + Ket::Sum.coordinates::<T>()[i] * c(-0.2, 0.6)
    });
    let expanded = p.expand(&NonOrthogonal).unwrap();
    assert_values(
        state_vector(&expanded.coordinates(&Canonical).unwrap()),
        reference,
    );
    assert_values(
        state_vector(&p.coordinates(&NonOrthogonal).unwrap()),
        reference,
    );
    assert_complex(expanded.inner_product(&expanded), dot(reference, reference));
    assert_complex(p.inner_product(&p), expanded.inner_product(&expanded));
    let b = Ls::<T>::primitive(Ket::Unit(Axis::B));
    let coords = b.coordinates(&NonOrthogonal).unwrap();
    let entries: Vec<_> = coords
        .flat_entries()
        .map(Result::unwrap)
        .map(|(i, c)| (i, *c))
        .collect();
    assert_eq!(entries.len(), 2);
    assert_eq!(entries[0].0, 0);
    assert_complex(entries[0].1, c(-1.0, 0.0));
    assert_complex(entries[1].1, c(1.0, 0.0));
    // Overlap with b0 is zero, whereas its expansion coordinate is -1.
    assert_complex(
        Ls::<T>::primitive(Ket::Unit(Axis::A)).inner_product(&b),
        c(0.0, 0.0),
    );
    assert!(std::ptr::eq(coords.basis(), &NonOrthogonal));
    assert_eq!(coords.dim().unwrap(), 3);

    for a in [Operator::Identity, Operator::Unit(Axis::A, Axis::B)] {
        for b in [Operator::Identity, Operator::Unit(Axis::B, Axis::C)] {
            let op = (Lo::<T>::primitive(a).clone() * &Lo::primitive(b))
                .adjoint()
                .scaled(c(0.4, -0.2));
            let expected = matmul(matrix::<T>(a), matrix(b));
            let expected = std::array::from_fn::<_, 3, _>(|i| {
                std::array::from_fn::<_, 3, _>(|j| expected[j][i].conj() * c(0.4, -0.2))
            });
            let expanded = op.expand(&MatrixUnits).unwrap();
            let coords = expanded.coordinates(&MatrixUnits).unwrap();
            let mut actual = [[c(0.0, 0.0); 3]; 3];
            for (index, coefficient) in coords.flat_entries().map(Result::unwrap) {
                actual[index / 3][index % 3] = *coefficient;
            }
            assert_matrix(actual, expected);
            let acted = p.act(&op).expand(&NonOrthogonal).unwrap();
            assert_values(
                state_vector(&acted.coordinates(&Canonical).unwrap()),
                matvec(expected, reference),
            );
        }
    }

    let f = HilbertSpace::factor(3u8);
    let g = HilbertSpace::factor(8u8);
    let space = HilbertSpace::tensor_all_with(
        kind,
        [
            HilbertSpace::tensor_all_with(kind, [f.clone()]).unwrap(),
            g.clone(),
        ],
    )
    .unwrap();
    let a = Ls::<T>::primitive(Ket::Phase).at(&f).unwrap();
    let b = Ls::<T>::primitive(Ket::Sum).at(&g).unwrap();
    let other = qmb_core::State::tensor_in(
        &space,
        [
            Ls::<T>::primitive(Ket::Unit(Axis::B)).at(&f).unwrap(),
            Ls::primitive(Ket::Unit(Axis::C)).at(&g).unwrap(),
        ],
    )
    .unwrap();
    let state = (qmb_core::State::tensor_in(&space, [a, b])
        .unwrap()
        .scaled(c(0.3, -0.7))
        .clone()
        + &other.scaled(c(-0.2, 0.6)))
        .unwrap();
    let placed = Lo::<T>::primitive(Operator::Unit(Axis::C, Axis::A))
        .at(&f)
        .unwrap();
    let state = state.act(&placed).unwrap();
    let reference: Vec<_> = (0..3)
        .flat_map(|i| {
            (0..3).map(move |j| {
                if i == 2 {
                    Ket::Sum.coordinates::<T>()[j] * c(0.3, -0.7)
                } else {
                    c(0.0, 0.0)
                }
            })
        })
        .collect();
    let expanded = state.expand(&NonOrthogonal).unwrap();
    assert!(expanded.space().shares_node_with(state.space()));
    let coords = expanded.coordinates(&Canonical).unwrap();
    assert_values(state_vector(&coords), reference);
    assert_values(
        state_vector(&state.coordinates(&NonOrthogonal).unwrap()),
        state_vector(&coords),
    );
    assert_complex(
        expanded.inner_product(&expanded).unwrap(),
        state.inner_product(&state).unwrap(),
    );

    // Identity is not a MatrixUnits basis primitive. Both spectator identities
    // must expand into sums of three diagonal units, rather than a fake index.
    let identity = qmb_core::Operator::<_, Operator, Complex<T>>::identity();
    let coords = identity.coordinates_in(&MatrixUnits, &space).unwrap();
    assert_eq!(coords.dim().unwrap(), 81);
    assert_eq!(coords.len(), 9);
    for (digits, value) in coords.entries() {
        assert!(digits.iter().all(|i| [0, 4, 8].contains(i)));
        assert_complex(*value, c(1.0, 0.0));
    }
    let expanded = identity.expand_in(&MatrixUnits, &space).unwrap();
    assert!(expanded.space().unwrap().shares_node_with(&space));
    assert_complex(expanded.trace().unwrap(), c(9.0, 0.0));
    let placed = Lo::<T>::primitive(Operator::Unit(Axis::A, Axis::B))
        .at(&f)
        .unwrap();
    let x = placed.expand_in(&MatrixUnits, &space).unwrap();
    assert_eq!(x.coordinates(&MatrixUnits).unwrap().len(), 3);
    assert_complex(x.trace().unwrap(), c(0.0, 0.0));

    let basis = space
        .basis_states::<Ket, Operator, Complex<T>, _>(&Canonical)
        .unwrap();
    for (index, state) in basis.enumerate() {
        let state = state.unwrap();
        let coords = state.coordinates(&Canonical).unwrap();
        assert_eq!(coords.len(), 1);
        let (actual, coefficient) = coords.flat_entries().next().unwrap().unwrap();
        assert_eq!(actual, index);
        assert_complex(*coefficient, c(1.0, 0.0));
    }
    assert_eq!(
        space
            .basis_operators::<Operator, Complex<T>, _>(&MatrixUnits)
            .unwrap()
            .count(),
        81
    );
    let regrouped = f.tensor_with(kind, &g).unwrap();
    let left: Vec<_> = space
        .basis_states::<Ket, Operator, Complex<T>, _>(&NonOrthogonal)
        .unwrap()
        .map(|s| state_vector(&s.unwrap().coordinates(&Canonical).unwrap()))
        .collect();
    let right: Vec<_> = regrouped
        .basis_states::<Ket, Operator, Complex<T>, _>(&NonOrthogonal)
        .unwrap()
        .map(|s| state_vector(&s.unwrap().coordinates(&Canonical).unwrap()))
        .collect();
    for (a, b) in left.into_iter().zip(right) {
        assert_values(a, b);
    }
}

#[test]
fn independent_three_level_expansion_at_both_precisions() {
    for kind in [qmb_core::TensorKind::Ordinary, qmb_core::TensorKind::Graded] {
        check_three_level::<f32>(kind);
        check_three_level::<f64>(kind);
    }
}

struct BrokenBasis(u8);
impl Basis<Ket, Complex<f64>> for BrokenBasis {
    fn expand(&self, _: &Ket) -> qmb_core::LinearCombination<Ket, Complex<f64>> {
        qmb_core::LinearCombination::primitive(Ket::Phase)
    }
    fn index(&self, _: &Ket) -> Option<usize> {
        match self.0 {
            0 => None,
            1 => Some(3),
            _ => Some(0),
        }
    }
    fn elements(&self) -> impl Iterator<Item = Ket> {
        [Ket::Phase].into_iter()
    }
}

#[test]
fn malformed_expansion_and_enumeration_are_errors_even_under_zero_scale() {
    let state = Ls::<f64>::primitive(Ket::Phase).scaled(c(0.0, 0.0));
    assert!(matches!(
        state.coordinates(&BrokenBasis(0)),
        Err(BasisError::NonBasisElement)
    ));
    assert!(matches!(
        state.expand(&BrokenBasis(1)),
        Err(BasisError::IndexOutOfBounds {
            index: 3,
            dimension: 3
        })
    ));
    let factor = HilbertSpace::factor(0u8);
    let mut elements = factor
        .basis_states::<Ket, Operator, Complex<f64>, _>(&BrokenBasis(2))
        .unwrap();
    assert!(elements.next().unwrap().is_ok());
    assert!(matches!(
        elements.next().unwrap(),
        Err(BasisError::InvalidEnumeration)
    ));
}

#[test]
fn zero_state_factor_does_not_hide_invalid_basis_expansion_in_another_factor() {
    let a = HilbertSpace::factor(0u8);
    let b = HilbertSpace::factor(1u8);
    let state = qmb_core::State::tensor_in(
        &a.tensor(&b).unwrap(),
        [
            qmb_core::State::<_, Ket, Operator, Complex<f64>>::zero(&a),
            Ls::<f64>::primitive(Ket::Phase).at(&b).unwrap(),
        ],
    )
    .unwrap();
    assert!(matches!(
        state.coordinates(&BrokenBasis(0)),
        Err(BasisError::NonBasisElement)
    ));
    assert!(matches!(
        state.expand(&BrokenBasis(1)),
        Err(BasisError::IndexOutOfBounds { .. })
    ));
}

struct BrokenOperatorBasis;
impl Basis<Operator, Complex<f64>> for BrokenOperatorBasis {
    fn expand(&self, p: &Operator) -> qmb_core::LinearCombination<Operator, Complex<f64>> {
        qmb_core::LinearCombination::primitive(*p)
    }
    fn index(&self, p: &Operator) -> Option<usize> {
        matches!(p, Operator::Identity).then_some(0)
    }
    fn elements(&self) -> impl Iterator<Item = Operator> {
        std::iter::once(Operator::Identity)
    }
}

#[test]
fn zero_operator_does_not_hide_invalid_basis_expansion_on_a_later_factor() {
    let a = HilbertSpace::factor(0u8);
    let b = HilbertSpace::factor(1u8);
    let space = a.tensor(&b).unwrap();
    let operator = Lo::<f64>::primitive(Operator::Unit(Axis::A, Axis::B))
        .at(&b)
        .unwrap()
        .scaled(c(0.0, 0.0));
    assert!(matches!(
        operator.coordinates_in(&BrokenOperatorBasis, &space),
        Err(BasisError::NonBasisElement)
    ));
    assert!(matches!(
        operator.expand_in(&BrokenOperatorBasis, &space),
        Err(BasisError::NonBasisElement)
    ));
}
