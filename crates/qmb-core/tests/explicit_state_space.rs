use qmb_core::{ExpressionError as E, HilbertSpace as H, State, StateOf, TensorKind};

// Construction requires none of Clone, numeric operations or evaluation rules.
struct Primitive;
struct Coeff;
type System = qmb_core::ExpressionFactory<Primitive, (), Coeff>;

#[derive(PartialEq, Eq)]
struct Index(u8);

fn leaf<I>(space: &H<I>) -> StateOf<System, I> {
    System::state(Primitive).at(space).unwrap()
}

#[test]
fn all_five_entries_share_the_exact_target_and_accept_different_grouping() {
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        let a = H::factor(Index(8));
        let b = H::factor(Index(2));
        let target = H::tensor_all_with(
            kind,
            [
                H::unit(),
                H::tensor_all_with(kind, [a.clone()]).unwrap(),
                b.clone(),
            ],
        )
        .unwrap();
        for state in [
            State::tensor_in(&target, [leaf(&a), leaf(&b)]),
            State::try_tensor_in(&target, [Ok(leaf(&a)), Ok(leaf(&b))]),
            System::tensor_in(&target, [leaf(&a), leaf(&b)]),
            System::try_tensor_in(&target, [Ok(leaf(&a)), Ok(leaf(&b))]),
            System::product_state(&target, [Primitive, Primitive]),
        ] {
            let state = state.unwrap();
            assert!(state.space().shares_node_with(&target));
            assert_eq!(state.space().children().len(), 3);
            assert_eq!(
                state.space().indices().map(|i| i.0).collect::<Vec<_>>(),
                [8, 2]
            );
        }
        // Equality of indices, not pointer equality, determines coverage.
        assert!(
            System::tensor_in(
                &target,
                [leaf(&H::factor(Index(8))), leaf(&H::factor(Index(2)))]
            )
            .is_ok()
        );
    }
}

#[test]
fn every_coverage_failure_uses_incompatible_spaces() {
    let a = H::factor(9);
    let b = H::factor(1);
    let c = H::factor(4);
    let target = H::tensor_all([a.clone(), b.clone()]).unwrap();
    for states in [
        vec![],
        vec![leaf(&a)],
        vec![leaf(&b), leaf(&a)],
        vec![leaf(&a), leaf(&a)],
        vec![leaf(&a), leaf(&c)],
        vec![leaf(&a), leaf(&b), leaf(&c)],
    ] {
        assert_eq!(
            System::tensor_in(&target, states).err(),
            Some(E::IncompatibleSpaces)
        );
    }
    for primitives in [
        vec![],
        vec![Primitive],
        vec![Primitive, Primitive, Primitive],
    ] {
        assert_eq!(
            System::product_state(&target, primitives).err(),
            Some(E::IncompatibleSpaces)
        );
    }
}

#[test]
fn unit_and_singleton_keep_the_given_handle_without_reinterpreting_kinds() {
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        let unit = H::<u8>::tensor_all_with(kind, [H::unit()]).unwrap();
        for state in [
            StateOf::<System, _>::tensor_in(&unit, []),
            StateOf::<System, _>::try_tensor_in(&unit, []),
            System::tensor_in(&unit, []),
            System::try_tensor_in(&unit, []),
            System::product_state(&unit, []),
        ] {
            assert!(state.unwrap().space().shares_node_with(&unit));
        }
        assert_eq!(
            System::product_state(&unit, [Primitive]).err(),
            Some(E::IncompatibleSpaces)
        );
        let a = H::factor(0);
        let wrapper = H::tensor_all_with(kind, [H::unit(), a.clone()]).unwrap();
        let state = System::tensor_in(&wrapper, [leaf(&a)]).unwrap();
        assert!(state.space().shares_node_with(&wrapper));
        assert_eq!(
            System::tensor_in(&wrapper, []).err(),
            Some(E::IncompatibleSpaces)
        );
        let zero = StateOf::<System, _>::zero(&wrapper);
        assert!(
            System::tensor_in(&wrapper, [zero])
                .unwrap()
                .space()
                .shares_node_with(&wrapper)
        );
    }
    assert!(
        StateOf::<System, u8>::unit()
            .space()
            .is_compatible_with(&H::unit())
    );
}

#[test]
fn multi_factor_inputs_must_match_the_target_kind_even_when_zero() {
    let a = H::factor(0);
    let b = H::factor(1);
    let ordinary = H::tensor_all([a.clone(), b.clone()]).unwrap();
    let graded = H::tensor_all_g([a.clone(), b.clone()]).unwrap();
    for (target, other) in [(&ordinary, &graded), (&graded, &ordinary)] {
        for state in [
            System::product_state(other, [Primitive, Primitive]).unwrap(),
            StateOf::<System, _>::zero(other),
        ] {
            assert_eq!(
                System::tensor_in(target, [state]).err(),
                Some(E::IncompatibleSpaces)
            );
        }
        let source = System::product_state(target, [Primitive, Primitive]).unwrap();
        let wrapper = H::tensor_all([H::unit(), target.clone()]).unwrap();
        assert!(
            System::tensor_in(&wrapper, [State::unit(), source])
                .unwrap()
                .space()
                .shares_node_with(&wrapper)
        );
    }
}

#[test]
fn fallible_entries_short_circuit_before_any_space_validation() {
    let a = H::factor(0);
    let target = H::factor(1);
    for via_system in [false, true] {
        let mut pulls = 0;
        let input = (0..).map(|n| {
            pulls += 1;
            match n {
                0 | 1 => Ok(leaf(&a)), // Both wrong index and duplicate.
                2 => Err(E::MissingIndex),
                _ => panic!("consumed after the first input error"),
            }
        });
        let result = if via_system {
            System::try_tensor_in(&target, input)
        } else {
            State::try_tensor_in(&target, input)
        };
        assert_eq!(result.err(), Some(E::MissingIndex));
        assert_eq!(pulls, 3);
    }
}

#[test]
fn composite_input_can_be_reused_with_an_explicit_compatible_target() {
    let factors = [H::factor(0), H::factor(1), H::factor(2)];
    let ab = H::tensor_all_g(factors[..2].iter().cloned()).unwrap();
    let abc = H::tensor_all_g([ab.clone(), factors[2].clone()]).unwrap();
    let block = System::product_state(&ab, [Primitive, Primitive]).unwrap();
    let state = System::tensor_in(&abc, [block, leaf(&factors[2])]).unwrap();
    assert!(state.space().shares_node_with(&abc));
}

#[test]
fn product_state_does_not_borrow_the_target_or_primitive_container() {
    let name = String::from("factor");
    let state = {
        let target = H::factor(name.as_str());
        let primitives = vec![Primitive];
        System::product_state(&target, primitives).unwrap()
    };
    assert_eq!(
        state.space().indices().copied().collect::<Vec<_>>(),
        ["factor"]
    );
}
