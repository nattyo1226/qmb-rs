//! Basis expansion checked against literal spin matrices and sparse large spaces.

use num_complex::Complex;
use num_traits::{FromPrimitive, Num, ToPrimitive};
use qmb_core::{
    BasisError, DimensionOverflow, ExpressionError, HilbertSpace, LocalOperator, LocalState,
    Operator, State,
};
use qmb_spin_half::{
    CartesianBasis, LadderBasis, SpinHalfOperator as Op, SpinHalfState as Ket, SzBasis,
};
use std::{fmt::Debug, ops::Neg};

trait Real: Num + Copy + Neg<Output = Self> + Debug + FromPrimitive + ToPrimitive {
    const TOLERANCE: f64;
}
impl Real for f32 {
    const TOLERANCE: f64 = 2e-5;
}
impl Real for f64 {
    const TOLERANCE: f64 = 1e-12;
}
fn c<T: Real>(r: f64, i: f64) -> Complex<T> {
    Complex::new(T::from_f64(r).unwrap(), T::from_f64(i).unwrap())
}
fn close<T: Real>(a: Complex<T>, b: Complex<T>) {
    for (a, b) in [(a.re, b.re), (a.im, b.im)] {
        let a = a.to_f64().unwrap();
        let b = b.to_f64().unwrap();
        assert!(
            (a - b).abs() <= T::TOLERANCE * (1.0 + b.abs()),
            "{a} != {b}"
        );
    }
}
type Lo<T> = LocalOperator<Op, Complex<T>>;
type Ls<T> = LocalState<Ket, Op, Complex<T>>;
type Matrix<T> = [[Complex<T>; 2]; 2];
fn matrix<T: Real>(p: Op) -> Matrix<T> {
    use Op::*;
    match p {
        I => [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]],
        Sx => [[c(0.0, 0.0), c(0.5, 0.0)], [c(0.5, 0.0), c(0.0, 0.0)]],
        Sy => [[c(0.0, 0.0), c(0.0, -0.5)], [c(0.0, 0.5), c(0.0, 0.0)]],
        Sz => [[c(0.5, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(-0.5, 0.0)]],
        SPlus => [[c(0.0, 0.0), c(1.0, 0.0)], [c(0.0, 0.0), c(0.0, 0.0)]],
        SMinus => [[c(0.0, 0.0), c(0.0, 0.0)], [c(1.0, 0.0), c(0.0, 0.0)]],
    }
}
fn check_local<T: Real>() {
    use Op::*;
    let basis = CartesianBasis;
    for a in [I, Sx, Sy, Sz, SPlus, SMinus] {
        for b in [I, Sx, Sy, Sz, SPlus, SMinus] {
            let x = Lo::<T>::primitive(a);
            let y = Lo::<T>::primitive(b);
            let op = ((x.clone() * &y) + x.scaled(c(0.3, -0.7)) + y.scaled(c(-0.2, 0.6)))
                .adjoint()
                .scaled(c(0.4, 0.2));
            let am = matrix::<T>(a);
            let bm = matrix::<T>(b);
            let m: Matrix<T> = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    ((am[j][0] * bm[0][i] + am[j][1] * bm[1][i])
                        + am[j][i] * c(0.3, -0.7)
                        + bm[j][i] * c(-0.2, 0.6))
                    .conj()
                        * c(0.4, 0.2)
                })
            });
            let expected = [
                (m[0][0] + m[1][1]) / c(2.0, 0.0),
                m[0][1] + m[1][0],
                c::<T>(0.0, 1.0) * (m[0][1] - m[1][0]),
                m[0][0] - m[1][1],
            ];
            let expanded = op.expand(&basis).unwrap();
            let via_ladder = op.expand(&LadderBasis).unwrap().expand(&basis).unwrap();
            for expr in [&op, &expanded, &via_ladder] {
                let coords = expr.coordinates(&basis).unwrap();
                assert_eq!(coords.dim().unwrap(), 4);
                let mut actual = [c(0.0, 0.0); 4];
                for (i, value) in coords.flat_entries().map(Result::unwrap) {
                    actual[i] = *value;
                }
                for (a, b) in actual.into_iter().zip(expected) {
                    close(a, b);
                }
            }
            let psi = Ls::<T>::primitive(Ket::Up).scaled(c(0.2, 0.7))
                + Ls::primitive(Ket::Down).scaled(c(-0.3, 0.4));
            let acted = psi.act(&op).expand(&SzBasis).unwrap();
            let coords = acted.coordinates(&SzBasis).unwrap();
            let mut actual = [c(0.0, 0.0); 2];
            for (i, value) in coords.flat_entries().map(Result::unwrap) {
                actual[i] = *value;
            }
            for i in 0..2 {
                close(actual[i], m[i][0] * c(0.2, 0.7) + m[i][1] * c(-0.3, 0.4));
            }
            close(psi.inner_product(&psi), c(0.78, 0.0)); // The input is unchanged.
        }
    }
    let zero =
        Lo::<T>::primitive(SPlus) + Lo::primitive(SMinus) - Lo::primitive(Sx).scaled(c(2.0, 0.0));
    assert!(zero.coordinates(&basis).unwrap().is_empty());
    assert!(
        zero.expand(&basis)
            .unwrap()
            .coordinates(&LadderBasis)
            .unwrap()
            .is_empty()
    );
    let tiny = Lo::<T>::primitive(Sx).scaled(c(1e-35, 0.0));
    assert_eq!(tiny.coordinates(&basis).unwrap().len(), 1);
    let nan = Lo::<T>::primitive(Sx).scaled(c(f64::NAN, 0.0));
    assert_eq!(nan.coordinates(&basis).unwrap().len(), 1);
}

#[test]
fn spin_rules_reexpansion_and_actions_match_matrices_at_both_precisions() {
    check_local::<f32>();
    check_local::<f64>();
}

fn check_composite<T: Real>() {
    let a = HilbertSpace::factor(0u8);
    let b = HilbertSpace::factor(1u8);
    let h = HilbertSpace::tensor_all([a.clone(), HilbertSpace::tensor_all([b.clone()]).unwrap()])
        .unwrap();
    let x = Lo::<T>::primitive(Op::Sx).at(&a).unwrap();
    let y = Lo::<T>::primitive(Op::Sy)
        .at(&b)
        .unwrap()
        .scaled(c(0.3, -0.7));
    let z = Lo::<T>::primitive(Op::Sz).at(&a).unwrap();
    let l = Lo::<T>::primitive(Op::SMinus)
        .at(&b)
        .unwrap()
        .scaled(c(-0.2, 0.6));
    let op = ((x.clone() + &y).unwrap().clone() * &(z.clone() + &l).unwrap())
        .unwrap()
        .adjoint();
    let am: [[Complex<T>; 4]; 4] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            matrix::<T>(Op::Sx)[i / 2][j / 2] * matrix::<T>(Op::I)[i % 2][j % 2]
                + c::<T>(0.3, -0.7)
                    * matrix::<T>(Op::I)[i / 2][j / 2]
                    * matrix::<T>(Op::Sy)[i % 2][j % 2]
        })
    });
    let bm: [[Complex<T>; 4]; 4] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            matrix::<T>(Op::Sz)[i / 2][j / 2] * matrix::<T>(Op::I)[i % 2][j % 2]
                + c::<T>(-0.2, 0.6)
                    * matrix::<T>(Op::I)[i / 2][j / 2]
                    * matrix::<T>(Op::SMinus)[i % 2][j % 2]
        })
    });
    let expected: [[Complex<T>; 4]; 4] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (0..4)
                .fold(c::<T>(0.0, 0.0), |sum, k| sum + am[j][k] * bm[k][i])
                .conj()
        })
    });
    let expanded = op.expand_in(&CartesianBasis, &h).unwrap();
    for expr in [
        expanded.clone(),
        op.expand_in(&LadderBasis, &h)
            .unwrap()
            .expand(&CartesianBasis)
            .unwrap(),
    ] {
        let coords = expr.coordinates(&CartesianBasis).unwrap();
        let mut actual = [[c(0.0, 0.0); 4]; 4];
        for (digits, coefficient) in coords.entries() {
            let ops = [Op::I, Op::Sx, Op::Sy, Op::Sz];
            for (i, row) in actual.iter_mut().enumerate() {
                for (j, v) in row.iter_mut().enumerate() {
                    *v = *v
                        + *coefficient
                            * matrix::<T>(ops[digits[0]])[i / 2][j / 2]
                            * matrix::<T>(ops[digits[1]])[i % 2][j % 2];
                }
            }
        }
        for (a, b) in actual
            .into_iter()
            .flatten()
            .zip(expected.into_iter().flatten())
        {
            close(a, b);
        }
        close(
            expr.trace().unwrap(),
            (0..4).fold(c(0.0, 0.0), |sum, i| sum + expected[i][i]),
        );
    }
    let ud = State::tensor_in(
        &h,
        [
            Ls::<T>::primitive(Ket::Up).at(&a).unwrap(),
            Ls::primitive(Ket::Down).at(&b).unwrap(),
        ],
    )
    .unwrap();
    let du = State::tensor_in(
        &h,
        [
            Ls::<T>::primitive(Ket::Down).at(&a).unwrap(),
            Ls::primitive(Ket::Up).at(&b).unwrap(),
        ],
    )
    .unwrap();
    let psi = (ud.scaled(c(0.3, 0.7)).clone() + &du.scaled(c(-0.2, 0.6))).unwrap();
    let acted = psi.act(&expanded).unwrap().expand(&SzBasis).unwrap();
    let coords = acted.coordinates(&SzBasis).unwrap();
    let mut actual = [c(0.0, 0.0); 4];
    for (i, c) in coords.flat_entries().map(Result::unwrap) {
        actual[i] = *c;
    }
    for i in 0..4 {
        close(
            actual[i],
            expected[i][1] * c(0.3, 0.7) + expected[i][2] * c(-0.2, 0.6),
        );
    }
    for (index, expr) in h
        .basis_operators::<Op, Complex<T>, _>(&CartesianBasis)
        .unwrap()
        .enumerate()
    {
        let expr = expr.unwrap();
        let coords = expr.coordinates(&CartesianBasis).unwrap();
        assert_eq!(coords.len(), 1);
        let (actual, value) = coords.flat_entries().next().unwrap().unwrap();
        assert_eq!(actual, index);
        close(*value, c(1.0, 0.0));
    }
}

#[test]
fn composite_products_adjoints_spectators_and_actions_match_four_by_four_matrices() {
    check_composite::<f32>();
    check_composite::<f64>();
}

#[test]
fn spaces_binding_spectators_zero_and_unit_are_preserved() {
    let a = HilbertSpace::factor("a");
    let b = HilbertSpace::factor("b");
    let h = a.tensor(&b).unwrap();
    let identity = Operator::<_, Op, Complex<f64>>::identity();
    assert!(matches!(
        identity.expand(&CartesianBasis),
        Err(BasisError::Expression(ExpressionError::UnboundOperator))
    ));
    assert!(matches!(
        identity.coordinates(&CartesianBasis),
        Err(BasisError::Expression(ExpressionError::UnboundOperator))
    ));
    let bound = identity.bind(&a).unwrap();
    let expanded = bound.expand_in(&CartesianBasis, &h).unwrap();
    assert!(expanded.space().unwrap().shares_node_with(&h));
    assert!(bound.space().unwrap().shares_node_with(&a));
    let coords = bound.coordinates_in(&CartesianBasis, &h).unwrap();
    assert_eq!(coords.dim().unwrap(), 16);
    assert_eq!(coords.len(), 1);
    assert_eq!(
        coords.entries().next().unwrap(),
        (&[0, 0][..], &c(1.0, 0.0))
    );
    assert!(coords.space().shares_node_with(&h));
    let zero = Operator::<_, Op, Complex<f64>>::zero()
        .expand_in(&CartesianBasis, &h)
        .unwrap();
    assert!(zero.space().unwrap().shares_node_with(&h));
    assert!(zero.coordinates(&CartesianBasis).unwrap().is_empty());
    let state = State::<_, Ket, Op, Complex<f64>>::zero(&h)
        .expand(&SzBasis)
        .unwrap();
    assert!(state.space().shares_node_with(&h));
    assert!(state.coordinates(&SzBasis).unwrap().is_empty());
    let invalid = Lo::<f64>::primitive(Op::Sx)
        .at(&b)
        .unwrap()
        .scaled(c(0.0, 0.0));
    assert!(matches!(
        invalid.coordinates_in(&CartesianBasis, &a),
        Err(BasisError::Expression(ExpressionError::MissingIndex))
    ));

    let unit = State::<u8, Ket, Op, Complex<f64>>::unit().scaled(c(0.3, 0.7));
    let unit = unit.expand(&SzBasis).unwrap();
    let coords = unit.coordinates(&SzBasis).unwrap();
    assert_eq!(coords.dim_local(), None);
    assert_eq!(coords.dim().unwrap(), 1);
    assert_eq!(coords.entries().next().unwrap(), (&[][..], &c(0.3, 0.7)));
    assert_eq!(
        coords.flat_entries().next().unwrap().unwrap(),
        (0, &c(0.3, 0.7))
    );
    let unit_space = HilbertSpace::<u8>::unit();
    let op = Operator::<_, Op, Complex<f64>>::scalar(c(0.4, -0.2))
        .expand_in(&CartesianBasis, &unit_space)
        .unwrap();
    assert_eq!(op.trace().unwrap(), c(0.4, -0.2));
    assert_eq!(
        unit_space
            .basis_states::<Ket, Op, Complex<f64>, _>(&SzBasis)
            .unwrap()
            .count(),
        1
    );
    assert_eq!(
        unit_space
            .basis_operators::<Op, Complex<f64>, _>(&CartesianBasis)
            .unwrap()
            .count(),
        1
    );

    // Coordinates borrow only their basis; they survive destruction of the input.
    let basis = CartesianBasis;
    let coords = Lo::<f64>::primitive(Op::SPlus).coordinates(&basis).unwrap();
    assert_eq!(
        coords.into_entries().collect::<Vec<_>>(),
        vec![(vec![1], c(1.0, 0.0)), (vec![2], c(0.0, 1.0))]
    );
}

#[test]
fn sparse_expansion_survives_global_dimension_and_flat_index_overflow() {
    let factors: Vec<_> = (0..usize::BITS as usize + 1)
        .map(HilbertSpace::factor)
        .collect();
    let space = HilbertSpace::tensor_all(factors.iter().cloned()).unwrap();
    let state = State::tensor_in(
        &space,
        factors
            .iter()
            .map(|f| Ls::<f64>::primitive(Ket::Down).at(f).unwrap()),
    )
    .unwrap();
    let expanded = state.expand(&SzBasis).unwrap();
    let coords = expanded.coordinates(&SzBasis).unwrap();
    assert_eq!(coords.len(), 1);
    assert_eq!(coords.dim(), Err(DimensionOverflow));
    assert_eq!(
        coords.flat_entries().next().unwrap(),
        Err(DimensionOverflow)
    );
    assert!(coords.entries().next().unwrap().0.iter().all(|i| *i == 1));
    let identity = Operator::<_, Op, Complex<f64>>::identity();
    let coords = identity
        .coordinates_in(&CartesianBasis, state.space())
        .unwrap();
    assert_eq!(coords.len(), 1);
    assert_eq!(coords.dim(), Err(DimensionOverflow));
    assert_eq!(
        coords.flat_entries().next().unwrap().unwrap(),
        (0, &c(1.0, 0.0))
    );
    // Producing a prefix does not enumerate or materialize the exponential basis.
    let prefix: Vec<_> = state
        .space()
        .basis_states::<Ket, Op, Complex<f64>, _>(&SzBasis)
        .unwrap()
        .take(2)
        .collect();
    assert_eq!(prefix.len(), 2);
    assert!(prefix.iter().all(Result::is_ok));
}

#[test]
fn deep_expansion_and_destruction_fit_a_small_stack() {
    std::thread::Builder::new()
        .stack_size(128 * 1024)
        .spawn(|| {
            let mut state = Ls::<f64>::primitive(Ket::Down);
            let mut op = Lo::<f64>::primitive(Op::SPlus);
            let mut space = HilbertSpace::factor(0u8);
            for _ in 0..20_000 {
                state = state.scaled(c(1.0, 0.0));
                op = op.adjoint();
                space = HilbertSpace::tensor_all([space]).unwrap();
            }
            assert_eq!(
                state
                    .expand(&SzBasis)
                    .unwrap()
                    .coordinates(&SzBasis)
                    .unwrap()
                    .len(),
                1
            );
            assert_eq!(
                op.expand(&CartesianBasis)
                    .unwrap()
                    .coordinates(&CartesianBasis)
                    .unwrap()
                    .len(),
                2
            );
            let f = space.children().next().unwrap();
            // A local placement must use the leaf, not the surrounding singleton groups.
            let mut leaf = f;
            while leaf.factor_index().is_none() {
                leaf = leaf.children().next().unwrap();
            }
            let placed = state.at(leaf).unwrap();
            let mut placed = placed;
            for _ in 0..20_000 {
                placed = State::tensor_in(&space, [State::unit(), placed]).unwrap();
            }
            let expanded = placed.expand(&SzBasis).unwrap();
            assert_eq!(expanded.coordinates(&SzBasis).unwrap().len(), 1);
            assert!(expanded.space().shares_node_with(placed.space()));
            let op = op
                .at(leaf)
                .unwrap()
                .expand_in(&CartesianBasis, &space)
                .unwrap();
            assert_eq!(op.coordinates(&CartesianBasis).unwrap().len(), 2);
        })
        .unwrap()
        .join()
        .unwrap();
}
