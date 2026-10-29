#[path = "support/graded.rs"]
mod graded;
use graded::OperatorOf;
use graded::*;
use qmb_core::{HilbertSpace as H, *};

#[test]
fn normalized_kinds_grouping_subsystems_and_errors() {
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        let [a, b, d] = [H::factor(8), H::factor(2), H::factor(5)];
        let ab = H::tensor_all_with(kind, [a.clone(), b.clone()]).unwrap();
        let abc = ab.tensor_with(kind, &d).unwrap();
        let regrouped = a
            .tensor_with(kind, &b.tensor_with(kind, &d).unwrap())
            .unwrap();
        assert!(abc.is_compatible_with(&regrouped));
        assert!(!abc.is_compatible_with(
            &H::tensor_all_with(kind, [d.clone(), b.clone(), a.clone()]).unwrap()
        ));
        assert!(
            H::tensor_all_with(kind, [a.clone(), d.clone()])
                .unwrap()
                .is_subsystem_of(&abc)
        );
        assert!(
            !H::tensor_all_with(kind, [d.clone(), a.clone()])
                .unwrap()
                .is_subsystem_of(&abc)
        );
        assert!(a.is_subsystem_of(&abc));
        assert!(H::unit().is_subsystem_of(&abc));
        assert!(matches!(
            ab.tensor_with(kind, &a),
            Err(SpaceError::DuplicateIndex)
        ));
    }
    let a = H::factor(0);
    let b = H::factor(1);
    let d = H::factor(2);
    let g = H::tensor_all_g([a.clone(), b.clone()]).unwrap();
    let o = H::tensor_all([a.clone(), b.clone()]).unwrap();
    assert!(!g.is_compatible_with(&o));
    assert!(!o.is_subsystem_of(&g));
    assert!(H::<u8>::unit().is_compatible_with(&H::tensor_all_g([]).unwrap()));
    assert!(
        H::tensor_all([H::tensor_all_g([]).unwrap(), a.clone()])
            .unwrap()
            .is_compatible_with(&H::tensor_all_g([a.clone()]).unwrap())
    );
    let wrapper = H::tensor_all([H::unit(), g.clone()]).unwrap();
    assert_eq!(wrapper.kind(), Some(TensorKind::Ordinary));
    assert_eq!(wrapper.effective_kind(), Some(TensorKind::Graded));
    assert_eq!(wrapper.children().len(), 2);
    assert!(wrapper.is_compatible_with(&g));
    assert!(matches!(
        H::tensor_all([d.clone(), H::tensor_all([g.clone()]).unwrap()]),
        Err(SpaceError::MixedTensorKinds)
    ));
    assert!(matches!(
        H::tensor_all_g([d, o]),
        Err(SpaceError::MixedTensorKinds)
    ));
    assert!(matches!(
        State::tensor_in(
            &H::tensor_all([a.clone(), H::factor(1), H::factor(3)]).unwrap(),
            [
                ket::<f64>(0, &[a.clone(), b], TensorKind::Graded),
                ket(0, &[H::factor(3)], TensorKind::Ordinary)
            ]
        ),
        Err(ExpressionError::IncompatibleSpaces)
    ));
    let bound = placed::<f64>(A, &a).bind(&g).unwrap();
    assert!(matches!(
        bound.clone()
            * &Operator::identity()
                .bind(&H::tensor_all([a.clone(), H::factor(1)]).unwrap())
                .unwrap(),
        Err(ExpressionError::IncompatibleSpaces)
    ));
    let bad = placed::<f64>(A, &H::factor(99)).scaled(c(0., 0.));
    assert!(matches!(
        bad.trace_in(&g),
        Err(ExpressionError::MissingIndex)
    ));
    assert!(bad.expand_in(&Units, &g).is_err());
}

fn matrix_of<T: Real>(op: &OperatorOf<T>, space: &H<u8>) -> Vec<num_complex::Complex<T>> {
    let n = space.dim::<Ket>().unwrap();
    let mut result = vec![c(0., 0.); n * n];
    for (j, input) in space.basis_states(&States).unwrap().enumerate() {
        let col = vector(&input.unwrap().act(op).unwrap());
        for i in 0..n {
            result[i * n + j] = col[i];
        }
    }
    result
}
fn check<T: Real>() {
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        for n in [2, 3] {
            // Deliberately not sorted by Index value.
            let factors: Vec<_> = [9, 1, 4].into_iter().take(n).map(H::factor).collect();
            let h = H::tensor_all_with(kind, factors.clone()).unwrap();
            let inputs: Vec<_> = h
                .basis_states::<Ket, Op, num_complex::Complex<T>, _>(&States)
                .unwrap()
                .map(Result::unwrap)
                .collect();
            for i in 0..n {
                for j in 0..n {
                    for a in OPS {
                        for b in OPS {
                            let lhs = placed::<T>(a, &factors[i]);
                            let rhs = placed(b, &factors[j]);
                            let op = (lhs.clone() * &rhs).unwrap().bind(&h).unwrap();
                            let expected =
                                product(&embedding(a, i, n, kind), &embedding(b, j, n, kind));
                            values(&matrix_of(&op, &h), &expected);
                            values(&matrix_of(&op.adjoint(), &h), &adjoint(&expected));
                            values(&matrix_of(&op.adjoint().adjoint(), &h), &expected);
                            close(op.trace().unwrap(), trace(&expected));
                            close(
                                (rhs.clone() * &lhs).unwrap().trace_in(&h).unwrap(),
                                trace(&expected),
                            );
                            close(
                                lhs.hs_inner_product_in(&rhs, &h).unwrap(),
                                trace(&product(
                                    &adjoint(&embedding(a, i, n, kind)),
                                    &embedding(b, j, n, kind),
                                )),
                            );
                            for basis in [false, true] {
                                let expanded = if basis {
                                    op.expand(&Parities).unwrap()
                                } else {
                                    op.expand(&Units).unwrap()
                                };
                                values(&matrix_of(&expanded, &h), &expected);
                                close(expanded.trace().unwrap(), trace(&expected));
                                let first = if basis {
                                    op.coordinates(&Parities)
                                        .unwrap()
                                        .into_entries()
                                        .collect::<Vec<_>>()
                                } else {
                                    op.coordinates(&Units).unwrap().into_entries().collect()
                                };
                                let again = if basis {
                                    expanded
                                        .expand(&Parities)
                                        .unwrap()
                                        .coordinates(&Parities)
                                        .unwrap()
                                        .into_entries()
                                        .collect::<Vec<_>>()
                                } else {
                                    expanded
                                        .expand(&Units)
                                        .unwrap()
                                        .coordinates(&Units)
                                        .unwrap()
                                        .into_entries()
                                        .collect()
                                };
                                assert_eq!(first.len(), again.len());
                                for ((ix, x), (iy, y)) in first.iter().zip(again) {
                                    assert_eq!(ix, &iy);
                                    close(*x, y);
                                }
                            }
                            let psi = &inputs[(1 << n) - 1];
                            values(
                                &vector(&psi.act(&op).unwrap()),
                                &vector(&psi.act(&rhs).unwrap().act(&lhs).unwrap()),
                            );
                            close(
                                inputs[0].inner_product(&psi.act(&op).unwrap()).unwrap(),
                                inputs[0]
                                    .act(&op.adjoint())
                                    .unwrap()
                                    .inner_product(psi)
                                    .unwrap(),
                            );
                        }
                    }
                }
            }
            // Mixed parity local operators, complex scales, and an entangled input.
            let mixed = (LocalOperator::<Op, num_complex::Complex<T>>::primitive(A)
                .scaled(c(0.2, 0.7))
                + LocalOperator::primitive(Op::P))
            .at(&factors[n - 1])
            .unwrap();
            let psi = (inputs[0].scaled(c(0.3, 0.2)).clone()
                + &inputs[(1 << n) - 1].scaled(c(-0.7, 0.4)))
                .unwrap();
            let m: Vec<_> = embedding::<T>(A, n - 1, n, kind)
                .into_iter()
                .zip(embedding(Op::P, n - 1, n, kind))
                .map(|(a, b)| a * c::<T>(0.2, 0.7) + b)
                .collect();
            let acted = psi.act(&mixed).unwrap();
            values(&vector(&acted), &act(&m, &vector(&psi)));
            values(&matrix_of(&mixed.expand_in(&Units, &h).unwrap(), &h), &m);
            values(
                &vector(&acted.expand(&States).unwrap().expand(&States).unwrap()),
                &vector(&acted),
            );
            close(psi.inner_product(&psi).unwrap(), c(0.78, 0.));
            for i in 0..n {
                for j in 0..n {
                    let a = placed::<T>(A, &factors[i]);
                    let b = placed(CREATE, &factors[j]);
                    if kind == TensorKind::Graded || i == j {
                        let car = a.anticommutator(&b).unwrap();
                        let expected = if i == j {
                            embedding(Op::I, 0, n, kind)
                        } else {
                            vec![c(0., 0.); 1 << (2 * n)]
                        };
                        values(&matrix_of(&car, &h), &expected);
                        values(
                            &matrix_of(&a.anticommutator(&placed(A, &factors[j])).unwrap(), &h),
                            &vec![c(0., 0.); 1 << (2 * n)],
                        );
                    }
                }
            }
            close(
                OperatorOf::<T>::identity().trace_in(&h).unwrap(),
                c((1 << n) as f64, 0.),
            );
        }
    }
}
#[test]
fn independent_two_and_three_factor_reference_f32() {
    check::<f32>();
}
#[test]
fn independent_two_and_three_factor_reference_f64() {
    check::<f64>();
}

#[test]
fn scope_and_intrinsic_coordinates_are_preserved() {
    let factors = [H::factor(0), H::factor(1)];
    let [a, b] = factors.clone();
    let h = H::tensor_all_g(factors.clone()).unwrap();
    let op = placed::<f64>(A, &b);
    let left = ket(1, &[a], TensorKind::Graded);
    let right = ket(1, &[b], TensorKind::Graded);
    let inner = State::tensor_in(&h, [left.clone(), right.act(&op).unwrap()]).unwrap();
    let inner_ready = State::tensor_in(
        &h,
        [left.clone(), right.act(&op).unwrap().evaluate().unwrap()],
    )
    .unwrap();
    let outer = State::tensor_in(&h, [left, right])
        .unwrap()
        .act(&op)
        .unwrap();
    values(&vector(&inner), &vector(&inner_ready));
    values(&vector(&outer), &vector(&inner.scaled(c(-1., 0.))));
    values(&vector(&inner.expand(&States).unwrap()), &vector(&inner));
    let coords = op.coordinates_in(&Parities, &h).unwrap();
    let entries: Vec<_> = coords.entries().collect();
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].0, &[0, 2]);
    close(*entries[0].1, c(1., 0.));
    values(&matrix_of(&op, &h), &kron(&local(Op::P), &local(A)));
    let old = op.bind(&H::tensor_all(factors.clone()).unwrap()).unwrap();
    values(
        &vector(
            &ket(3, &factors, TensorKind::Graded)
                .act_in(&old, &h)
                .unwrap(),
        ),
        &[c(0., 0.), c(0., 0.), c(-1., 0.), c(0., 0.)],
    );
    close(
        (op.adjoint().clone() * &op).unwrap().trace_in(&h).unwrap(),
        c(2., 0.),
    );
}

#[test]
fn constructors_fallibility_and_deep_graded_trees() {
    type Q = qmb_core::ExpressionFactory<Ket, Op, num_complex::Complex<f64>>;
    let a = H::factor(0);
    let b = H::factor(1);
    let g = Q::product_state(
        &H::tensor_all_g([a.clone(), b.clone()]).unwrap(),
        [Ket(1), Ket(1)],
    )
    .unwrap();
    let explicit = Q::product_state(g.space(), [Ket(1), Ket(1)]).unwrap();
    assert!(g.space().is_compatible_with(explicit.space()));
    values(&vector(&g), &vector(&explicit));
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        let x = Q::state(Ket(1)).at(&a).unwrap();
        let y = Q::state(Ket(1)).at(&b).unwrap();
        let target = H::tensor_all_with(kind, [a.clone(), b.clone()]).unwrap();
        let candidates = [
            Q::tensor_in(&target, [x.clone(), y.clone()]).unwrap(),
            Q::try_tensor_in(&target, [Ok(x.clone()), Ok(y.clone())]).unwrap(),
            State::tensor_in(&target, [x.clone(), y.clone()]).unwrap(),
            State::try_tensor_in(&target, [Ok(x), Ok(y)]).unwrap(),
            Q::product_state(&target, [Ket(1), Ket(1)]).unwrap(),
        ];
        for state in candidates {
            assert!(state.space().is_compatible_with(&target));
            values(
                &vector(&state),
                &[c(0., 0.), c(0., 0.), c(0., 0.), c(1., 0.)],
            );
        }
    }
    let single = Q::state(Ket(0)).at(&a).unwrap();
    let mut count = 0;
    let inputs = [
        Ok(single.clone()),
        Ok(single),
        Err(ExpressionError::MissingIndex),
        Err(ExpressionError::NotAFactor),
    ]
    .into_iter()
    .inspect(|_| count += 1);
    assert!(matches!(
        Q::try_tensor_in(g.space(), inputs),
        Err(ExpressionError::MissingIndex)
    ));
    assert_eq!(count, 3);
    std::thread::Builder::new()
        .stack_size(128 * 1024)
        .spawn(|| {
            let factors = [H::factor(0), H::factor(1)];
            let mut psi = ket::<f64>(3, &factors, TensorKind::Graded);
            for _ in 0..20_000 {
                psi = State::tensor_in(
                    &H::tensor_all_g(factors.clone()).unwrap(),
                    [State::unit(), psi],
                )
                .unwrap();
            }
            assert_eq!(psi.space().kind(), Some(TensorKind::Graded));
            assert!(
                psi.space()
                    .is_compatible_with(&H::tensor_all_g(factors.clone()).unwrap())
            );
            let out = psi
                .act(&placed(A, &factors[1]))
                .unwrap()
                .evaluate()
                .unwrap();
            close(
                out.norm_squared().map(|x| c::<f64>(x, 0.)).unwrap(),
                c(1., 0.),
            );
            drop(out);
            drop(psi);
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn graded_action_needs_neither_product_nor_numeric_coefficient_operations() {
    use std::cell::Cell;
    #[derive(Clone)]
    struct Scalar;
    impl Conjugate for Scalar {
        type Output = Self;
        fn conjugate(&self) -> Self {
            self.clone()
        }
    }
    enum Action {
        Flip,
        P,
    }
    impl HomogeneousParity for Action {
        fn parity(&self) -> Parity {
            match self {
                Self::Flip => Parity::Odd,
                Self::P => Parity::Even,
            }
        }
    }
    impl ParityOperator for Action {
        fn parity_operator() -> Self {
            Self::P
        }
    }
    impl Adjoint<Scalar> for Action {
        type Term = Self;
        fn adjoint(&self) -> LinearCombination<Self, Scalar> {
            LinearCombination::term(
                Scalar,
                match self {
                    Self::Flip => Self::Flip,
                    Self::P => Self::P,
                },
            )
        }
    }
    struct Borrowed<'a>(&'a Cell<usize>);
    impl<'a> Act<Scalar> for Borrowed<'a> {
        type Operator = Action;
        type Term = Self;
        fn act(&self, op: &Action) -> LinearCombination<Self, Scalar> {
            if matches!(op, Action::P) {
                self.0.set(self.0.get() + 1);
            }
            LinearCombination::term(Scalar, Self(self.0))
        }
    }
    let count = Cell::new(0);
    let factors = [H::factor(0), H::factor(1), H::factor(2)];
    let space = H::tensor_all_g(factors.iter().cloned()).unwrap();
    let psi = State::tensor_in(
        &space,
        factors
            .iter()
            .map(|f| LocalState::primitive(Borrowed(&count)).at(f).unwrap()),
    )
    .unwrap();
    let op = LocalOperator::primitive(Action::Flip)
        .at(&factors[2])
        .unwrap();
    let result = psi.act(&op.adjoint()).unwrap().evaluate().unwrap();
    assert_eq!(count.get(), 2);
    assert!(result.space().is_compatible_with(psi.space()));
}

#[test]
fn large_separable_states_and_traces_do_not_need_global_dimensions_or_bases() {
    let factors: Vec<_> = (0..80).map(H::factor).collect();
    let psi = State::<_, Ket, Op, num_complex::Complex<f64>>::tensor_in(
        &H::tensor_all_g(factors.iter().cloned()).unwrap(),
        factors.iter().map(|f| {
            (LocalState::primitive(Ket(0)) + LocalState::primitive(Ket(1)))
                .at(f)
                .unwrap()
        }),
    )
    .unwrap();
    assert!(psi.space().dim::<Ket>().is_err());
    let a = placed(A, &factors[79]);
    let result = psi.act(&a).unwrap().evaluate().unwrap();
    close(result.inner_product(&result).unwrap(), c(2f64.powi(79), 0.));
    close(
        (a.adjoint().clone() * &a)
            .unwrap()
            .trace_in(psi.space())
            .unwrap(),
        c(2f64.powi(79), 0.),
    );
}

#[test]
fn local_expression_homomorphism_three_odd_factors_and_empty_bases() {
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        let factors = [H::factor(0), H::factor(1), H::factor(2)];
        let h = H::tensor_all_with(kind, factors.clone()).unwrap();
        let local = qmb_core::LocalOperator::<Op, num_complex::Complex<f64>>::primitive(A)
            .scaled(c(0.3, 0.4))
            + LocalOperator::primitive(Op::P);
        let local_product = (local.clone() * &local.adjoint()).adjoint();
        let placed_local = local_product.at(&factors[2]).unwrap();
        let lifted = local.at(&factors[2]).unwrap();
        let placed_product = (lifted.clone() * &lifted.adjoint()).unwrap().adjoint();
        values(
            &matrix_of(&placed_local, &h),
            &matrix_of(&placed_product, &h),
        );
        values(
            &matrix_of(&placed_local.expand_in(&Units, &h).unwrap(), &h),
            &matrix_of(&placed_product, &h),
        );
        for order in [[0, 1, 2], [2, 1, 0], [1, 2, 0]] {
            let mut op = OperatorOf::<f64>::identity();
            let mut m = embedding(Op::I, 0, 3, kind);
            for i in order {
                op = (op.clone() * &placed(A, &factors[i])).unwrap();
                m = product(&m, &embedding(A, i, 3, kind));
            }
            values(&matrix_of(&op, &h), &m);
            values(
                &matrix_of(&op.adjoint().expand_in(&Units, &h).unwrap(), &h),
                &adjoint(&m),
            );
        }
        for (i, op) in h
            .basis_operators::<Op, num_complex::Complex<f64>, _>(&Parities)
            .unwrap()
            .enumerate()
        {
            let coords = op.unwrap().coordinates(&Parities).unwrap();
            assert_eq!(coords.len(), 1);
            let (index, value) = coords.flat_entries().next().unwrap().unwrap();
            assert_eq!(index, i);
            close(*value, c(1., 0.));
        }
        let unit = H::<u8>::tensor_all_with(kind, []).unwrap();
        let one = State::<u8, Ket, Op, num_complex::Complex<f64>>::tensor_in(&unit, []).unwrap();
        assert_eq!(one.space().kind(), Some(kind));
        close(one.inner_product(&one).unwrap(), c(1., 0.));
        let state = unit
            .basis_states::<Ket, Op, num_complex::Complex<f64>, _>(&States)
            .unwrap()
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(state.coordinates(&States).unwrap().factor_count(), 0);
        assert_eq!(
            unit.basis_operators::<Op, num_complex::Complex<f64>, _>(&Units)
                .unwrap()
                .count(),
            1
        );
        close(
            OperatorOf::<f64>::identity().trace_in(&unit).unwrap(),
            c(1., 0.),
        );
    }
}
