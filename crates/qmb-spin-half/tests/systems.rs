use num_complex::Complex;
use qmb_core::{ExpressionError, HilbertSpace as H, SpaceError, State, StateOf};
use qmb_geometry::{Boundary, Chain, ChainError};
use qmb_spin_half::{SpinHalf, SpinHalfOperator as Op, SpinHalfState as Ket, UtilityError};

type Spin = SpinHalf<Complex<f64>>;

#[test]
fn local_and_composite_expressions_connect_to_existing_types() {
    let a = H::factor("a");
    let b = H::factor("b");
    let up: qmb_core::LocalStateOf<Spin> = Spin::state(Ket::Up);
    let sz: qmb_core::LocalOperatorOf<Spin> = Spin::operator(Op::Sz);
    let psi: StateOf<Spin, _> = State::tensor_in(
        &a.tensor(&b).unwrap(),
        [up.at(&a).unwrap(), Spin::state(Ket::Down).at(&b).unwrap()],
    )
    .unwrap();
    let zz: qmb_core::OperatorOf<Spin, _> =
        (sz.at(&a).unwrap().clone() * &sz.at(&b).unwrap()).unwrap();
    assert_eq!(psi.expectation(&zz).unwrap(), Complex::new(-0.25, 0.));
    // Construction fixes f32 even without downstream evaluation or annotation.
    let single = SpinHalf::<Complex<f32>>::state(Ket::Up);
    let operator = SpinHalf::<Complex<f32>>::operator(Op::Sx);
    assert_eq!(
        single
            .act(&operator)
            .inner_product(&SpinHalf::<Complex<f32>>::state(Ket::Down)),
        Complex::new(0.5_f32, 0.)
    );
    let chain = Chain::new(2, Boundary::Open);
    let sites = [a, b];
    let h = qmb_spin_half::models::SpinHalfModels::<SpinHalf<Complex<f32>>>::heisenberg(
        &chain, &sites, 1.0,
    )
    .unwrap();
    let psi = qmb_spin_half::models::SpinHalfModels::<SpinHalf<Complex<f32>>>::neel(
        &chain,
        &sites,
        Ket::Up,
    )
    .unwrap();
    assert_eq!(psi.expectation(&h).unwrap(), Complex::new(0.25_f32, 0.));
}

#[test]
fn models_accept_multiple_leaves_and_alternate_by_site() {
    let a = H::factor(0);
    let b = H::factor(1);
    let c = H::factor(2);
    let nested =
        H::tensor_all([H::unit(), H::tensor_all([a.clone(), b.clone()]).unwrap()]).unwrap();
    let sites = [nested, c.clone()];
    let chain = Chain::new(2, Boundary::Open);
    let h = qmb_spin_half::models::SpinHalfModels::<Spin>::heisenberg(&chain, &sites, 1.0).unwrap();
    let xyz =
        qmb_spin_half::models::SpinHalfModels::<Spin>::xyz(&chain, &sites, 1.0, 1.0, 1.0).unwrap();
    let psi = qmb_spin_half::models::SpinHalfModels::<Spin>::neel(&chain, &sites, Ket::Up).unwrap();
    let expected = State::tensor_in(
        psi.space(),
        [
            Spin::state(Ket::Up).at(&a).unwrap(),
            Spin::state(Ket::Up).at(&b).unwrap(),
            Spin::state(Ket::Down).at(&c).unwrap(),
        ],
    )
    .unwrap();
    assert_eq!(psi.inner_product(&expected).unwrap(), Complex::new(1., 0.));
    assert_eq!(psi.expectation(&h).unwrap(), Complex::new(0.5, 0.));
    // Compare every matrix element against explicitly placed cross-site terms.
    let mut terms = Vec::new();
    for left in [&a, &b] {
        for axis in [Op::Sx, Op::Sy, Op::Sz] {
            let local = Spin::operator(axis);
            terms.push(
                (local.at(left).unwrap().clone() * &local.at(&c).unwrap())
                    .unwrap()
                    .scaled(Complex::new(-1., 0.)),
            );
        }
    }
    let reference = qmb_core::OperatorOf::<Spin, _>::sum(terms).unwrap();
    let basis: Vec<_> = (0..8)
        .map(|bits| {
            State::tensor_in(
                psi.space(),
                [&a, &b, &c].into_iter().enumerate().map(|(i, leaf)| {
                    Spin::state(if bits >> (2 - i) & 1 == 0 {
                        Ket::Up
                    } else {
                        Ket::Down
                    })
                    .at(leaf)
                    .unwrap()
                }),
            )
            .unwrap()
        })
        .collect();
    for left in &basis {
        for right in &basis {
            let value = left.inner_product(&right.act(&h).unwrap()).unwrap();
            assert_eq!(
                value,
                left.inner_product(&right.act(&reference).unwrap()).unwrap()
            );
            assert_eq!(
                value,
                left.inner_product(&right.act(&xyz).unwrap()).unwrap()
            );
        }
    }
    let full = h.space().unwrap();
    assert_eq!(full.indices().copied().collect::<Vec<_>>(), [0, 1, 2]);
    assert!(full.children().next().unwrap().shares_node_with(&sites[0]));
    assert_eq!(full.children().next().unwrap().children().count(), 2);
}

#[test]
fn empty_sites_keep_edges_and_parity_and_all_empty_inputs_give_units() {
    let sites = [H::factor(0), H::unit(), H::factor(1)];
    let chain = Chain::new(3, Boundary::Open);
    let psi = qmb_spin_half::models::SpinHalfModels::<Spin>::neel(&chain, &sites, Ket::Up).unwrap();
    let expected = State::tensor_in(
        psi.space(),
        [
            Spin::state(Ket::Up).at(&sites[0]).unwrap(),
            Spin::state(Ket::Up).at(&sites[2]).unwrap(),
        ],
    )
    .unwrap();
    assert_eq!(psi.inner_product(&expected).unwrap(), Complex::new(1., 0.));
    assert_eq!(
        psi.expectation(
            &qmb_spin_half::models::SpinHalfModels::<Spin>::heisenberg(&chain, &sites, 1.0)
                .unwrap()
        )
        .unwrap(),
        Complex::new(0., 0.)
    );
    // The periodic closing edge still couples site 0 to site 2.
    assert_eq!(
        psi.expectation(
            &qmb_spin_half::models::SpinHalfModels::<Spin>::heisenberg(
                &Chain::new(3, Boundary::Periodic),
                &sites,
                1.0
            )
            .unwrap()
        )
        .unwrap(),
        Complex::new(-0.25, 0.)
    );
    for count in 0..=3 {
        let sites: Vec<_> = (0..count).map(|_| H::<u8>::unit()).collect();
        let chain = Chain::new(count, Boundary::Periodic);
        let psi =
            qmb_spin_half::models::SpinHalfModels::<Spin>::neel(&chain, &sites, Ket::Down).unwrap();
        let h =
            qmb_spin_half::models::SpinHalfModels::<Spin>::heisenberg(&chain, &sites, 1.0).unwrap();
        assert_eq!(
            psi.inner_product(&State::unit()).unwrap(),
            Complex::new(1., 0.)
        );
        assert_eq!(h.trace().unwrap(), Complex::new(0., 0.));
        assert_eq!(h.space().unwrap().children().count(), count);
    }
    // A single nonempty site is retained in a Hamiltonian with no terms.
    let h = qmb_spin_half::models::SpinHalfModels::<Spin>::heisenberg(
        &Chain::new(1, Boundary::Open),
        &sites[..1],
        1.0,
    )
    .unwrap();
    assert_eq!(h.space().unwrap().indices().count(), 1);
}

#[test]
fn every_model_validates_length_and_duplicate_indices_without_edges() {
    for n in 0..=3 {
        for count in 0..=4 {
            let sites: Vec<_> = (0..count).map(H::factor).collect();
            let chain = Chain::new(n, Boundary::Open);
            let expected =
                (n != count).then_some(UtilityError::Chain(ChainError::SiteCountMismatch));
            assert_eq!(
                qmb_spin_half::models::SpinHalfModels::<Spin>::neel(&chain, &sites, Ket::Up).err(),
                expected
            );
            assert_eq!(
                qmb_spin_half::models::SpinHalfModels::<Spin>::heisenberg(&chain, &sites, 1.0)
                    .err(),
                expected
            );
            assert_eq!(
                qmb_spin_half::models::SpinHalfModels::<Spin>::xyz(&chain, &sites, 1.0, 2.0, 3.0)
                    .err(),
                expected
            );
        }
    }
    for boundary in [Boundary::Open, Boundary::Periodic] {
        // A middle empty site prevents all edges; duplicates must still fail.
        let sites = [
            H::factor(0),
            H::unit(),
            H::tensor_all([H::factor(0)]).unwrap(),
        ];
        let chain = Chain::new(3, boundary);
        let expected = Some(UtilityError::Space(SpaceError::DuplicateIndex));
        assert_eq!(
            qmb_spin_half::models::SpinHalfModels::<Spin>::neel(&chain, &sites, Ket::Up).err(),
            expected
        );
        assert_eq!(
            qmb_spin_half::models::SpinHalfModels::<Spin>::heisenberg(&chain, &sites, 1.0).err(),
            expected
        );
        assert_eq!(
            qmb_spin_half::models::SpinHalfModels::<Spin>::xyz(&chain, &sites, 1.0, 2.0, 3.0).err(),
            expected
        );
    }
}

#[test]
fn model_construction_needs_only_operation_specific_coefficient_capabilities() {
    struct NoCapabilities;
    #[derive(Clone)]
    struct Coupling(i32);
    impl From<i32> for Coupling {
        fn from(value: i32) -> Self {
            Self(value)
        }
    }
    impl core::ops::Neg for Coupling {
        type Output = Self;
        fn neg(self) -> Self {
            Self(-self.0)
        }
    }
    // No numeric traits, Clone or Act are required for the Neel coefficient.
    let sites = [H::factor(0), H::factor(1)];
    let chain = Chain::new(2, Boundary::Open);
    let _ = qmb_spin_half::models::SpinHalfModels::<SpinHalf<NoCapabilities>>::neel(
        &chain,
        &sites,
        Ket::Up,
    )
    .unwrap();
    // Hamiltonian construction needs no addition, multiplication, Zero or One.
    // XYZ's real coupling source need not implement Clone either.
    struct Real(i32);
    impl From<Real> for Coupling {
        fn from(value: Real) -> Self {
            Self(value.0)
        }
    }
    let _ = qmb_spin_half::models::SpinHalfModels::<SpinHalf<Coupling>>::xyz(
        &chain,
        &sites,
        Real(1),
        Real(2),
        Real(3),
    )
    .unwrap();
    let _ =
        qmb_spin_half::models::SpinHalfModels::<SpinHalf<Coupling>>::heisenberg(&chain, &sites, 1)
            .unwrap();
}

#[test]
fn non_clone_indices_and_deep_site_grouping_are_supported_iteratively() {
    #[derive(PartialEq, Eq)]
    struct Index(u8);
    std::thread::Builder::new()
        .stack_size(128 * 1024)
        .spawn(|| {
            let mut first = H::factor(Index(1));
            for _ in 0..20_000 {
                first = H::tensor_all([first]).unwrap();
            }
            let sites = [first, H::factor(Index(0))];
            let chain = Chain::new(2, Boundary::Periodic);
            let state =
                qmb_spin_half::models::SpinHalfModels::<Spin>::neel(&chain, &sites, Ket::Up)
                    .unwrap();
            let h = qmb_spin_half::models::SpinHalfModels::<Spin>::heisenberg(&chain, &sites, 1.0)
                .unwrap();
            assert!(
                h.space()
                    .unwrap()
                    .children()
                    .next()
                    .unwrap()
                    .shares_node_with(&sites[0])
            );
            assert_eq!(
                state.space().indices().map(|i| i.0).collect::<Vec<_>>(),
                [1, 0]
            );
            assert_eq!(state.expectation(&h).unwrap(), Complex::new(0.25, 0.));
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn model_errors_preserve_sources_and_messages() {
    for (error, source) in [
        (
            UtilityError::from(ExpressionError::NotAFactor),
            ExpressionError::NotAFactor.to_string(),
        ),
        (
            UtilityError::from(SpaceError::DuplicateIndex),
            SpaceError::DuplicateIndex.to_string(),
        ),
        (
            UtilityError::from(ChainError::SiteCountMismatch),
            ChainError::SiteCountMismatch.to_string(),
        ),
    ] {
        assert_eq!(error.to_string(), source);
        assert_eq!(
            std::error::Error::source(&error).unwrap().to_string(),
            source
        );
    }
}

macro_rules! construction_precision {
    ($name:ident, $real:ty) => {
        #[test]
        fn $name() {
            type Q = SpinHalf<Complex<$real>>;
            let a = H::factor(0u8);
            let b = H::factor(1);
            let space = a.tensor(&b).unwrap();
            let state = Q::product_state(
                &qmb_core::HilbertSpace::tensor_all([(a).clone(), (b).clone()]).unwrap(),
                [Ket::Up, Ket::Down],
            )
            .unwrap();
            let explicit = Q::tensor_in(
                &space,
                [
                    Q::state(Ket::Up).at(&a).unwrap(),
                    Q::state(Ket::Down).at(&b).unwrap(),
                ],
            )
            .unwrap();
            let fallible = Q::try_tensor_in(
                &space,
                [Q::state(Ket::Up).at(&a), Q::state(Ket::Down).at(&b)],
            )
            .unwrap();
            assert_eq!(
                state.inner_product(&explicit).unwrap(),
                Complex::new(1., 0.)
            );
            assert_eq!(
                state.inner_product(&fallible).unwrap(),
                Complex::new(1., 0.)
            );
            let h = Q::try_operator_sum(
                [Op::Sx, Op::Sy, Op::Sz]
                    .into_iter()
                    .map(|axis| Q::operator_product_at([(&a, axis), (&b, axis)])),
            )
            .unwrap()
            .bind(&space)
            .unwrap();
            assert_eq!(state.expectation(&h).unwrap(), Complex::new(-0.25, 0.));
            let sum = Q::state_sum(&space, [state.scaled(2.0), state.scaled(-1.0)]).unwrap();
            let tried = Q::try_state_sum(&space, [Ok(sum)]).unwrap();
            assert_eq!(state.inner_product(&tried).unwrap(), Complex::new(1., 0.));
            // Both scalar input forms retain the precision selected by Q.
            let half = Q::state(Ket::Up).scaled(0.5);
            let complex = Q::state(Ket::Up).scaled(Complex::new(0.5, 0.0));
            assert_eq!(half.inner_product(&complex), Complex::new(0.25, 0.));
            let placed = state.scaled(0.5);
            assert_eq!(placed.inner_product(&state).unwrap(), Complex::new(0.5, 0.));
            let sx = Q::operator(Op::Sx).scaled(0.5);
            assert_eq!(sx.trace(), Complex::new(0., 0.));
            assert_eq!(
                state.expectation(&h.scaled(0.5)).unwrap(),
                Complex::new(-0.125, 0.)
            );
            let lc = qmb_core::LinearCombination::<_, Complex<$real>>::primitive(0u8).scaled(0.5);
            assert_eq!(
                lc.terms().next().unwrap().coefficient,
                Complex::new(0.5, 0.)
            );
            // Same-index products preserve order: Sx Sy = i Sz / 2.
            let up = Q::state(Ket::Up).at(&a).unwrap();
            let xy = Q::operator_product_at([(&a, Op::Sx), (&a, Op::Sy)]).unwrap();
            let yx = Q::operator_product_at([(&a, Op::Sy), (&a, Op::Sx)]).unwrap();
            assert_eq!(up.expectation(&xy).unwrap(), Complex::new(0., 0.25));
            assert_eq!(up.expectation(&yx).unwrap(), Complex::new(0., -0.25));
            let explicit = Q::operator_product([
                Q::operator(Op::Sx).at(&a).unwrap(),
                Q::operator(Op::Sy).at(&a).unwrap(),
            ])
            .unwrap();
            let fallible =
                Q::try_operator_product([Q::operator(Op::Sx).at(&a), Q::operator(Op::Sy).at(&a)])
                    .unwrap();
            assert_eq!(
                up.expectation(&explicit).unwrap(),
                up.expectation(&xy).unwrap()
            );
            assert_eq!(
                up.expectation(&fallible).unwrap(),
                up.expectation(&xy).unwrap()
            );
            for unit in [
                Q::tensor_in(&H::<u8>::unit(), []),
                Q::try_tensor_in(&H::<u8>::unit(), []),
                Q::product_state(&H::<u8>::unit(), []),
            ] {
                let unit = unit.unwrap();
                assert_eq!(unit.inner_product(&unit).unwrap(), Complex::new(1., 0.));
            }
            for zero in [Q::state_sum(&space, []), Q::try_state_sum(&space, [])] {
                assert_eq!(zero.unwrap().norm_squared().unwrap(), 0.);
            }
            for zero in [Q::operator_sum::<u8>([]), Q::try_operator_sum::<u8>([])] {
                assert_eq!(
                    zero.unwrap().bind(&space).unwrap().trace().unwrap(),
                    Complex::new(0., 0.)
                );
            }
            for identity in [
                Q::operator_product::<u8>([]),
                Q::try_operator_product::<u8>([]),
                Q::operator_product_at::<u8>([]),
            ] {
                assert_eq!(
                    identity.unwrap().bind(&space).unwrap().trace().unwrap(),
                    Complex::new(4., 0.)
                );
            }
        }
    };
}
construction_precision!(construction_f32, f32);
construction_precision!(construction_f64, f64);
