use qmb_core::{
    HilbertSpace, LocalOperator, LocalOperatorOf, LocalState, LocalStateOf, OperatorOf,
    QuantumSystem, StateOf,
};

// These types deliberately provide no capabilities (including Clone).
struct Ket;
struct Op;
struct Coeff;
struct System;
struct Equivalent;
macro_rules! system {
    ($q:ty) => {
        impl QuantumSystem for $q {
            type State = Ket;
            type Operator = Op;
            type Coeff = Coeff;
        }
    };
}
system!(System);
system!(Equivalent);

#[test]
fn independent_construction_and_aliases_need_no_evaluation_capabilities() {
    fn state<Q: QuantumSystem>(primitive: Q::State) -> LocalStateOf<Q> {
        Q::state(primitive)
    }
    let local: LocalState<Ket, Op, Coeff> = state::<System>(Ket);
    let operator: LocalOperator<Op, Coeff> = System::operator(Op);
    let _: LocalOperatorOf<System> = operator.clone();
    let _: LocalStateOf<System> = local.clone();
    let factor = HilbertSpace::factor(0);
    let state: StateOf<System, _> = local.at(&factor).unwrap();
    let op: OperatorOf<System, _> = operator.at(&factor).unwrap();
    // The chosen system is not an expression identity tag.
    let other = Equivalent::state(Ket).at(&factor).unwrap();
    let _: StateOf<Equivalent, _> = state.add(&other).unwrap();
    let _: OperatorOf<Equivalent, _> = op
        .product(&Equivalent::operator(Op).at(&factor).unwrap())
        .unwrap();
    let _ = local.act(&operator);
}

#[test]
fn all_system_builders_accept_non_clone_primitives_coefficients_and_indices() {
    #[derive(PartialEq, Eq)]
    struct Index(u8);
    let a = HilbertSpace::factor(Index(3));
    let b = HilbertSpace::factor(Index(1));
    let space = a.tensor_product(&b).unwrap();
    let tensor = System::state_tensor_at([(&a, Ket), (&b, Ket)]).unwrap();
    assert_eq!(
        tensor.space().indices().map(|i| i.0).collect::<Vec<_>>(),
        [3, 1]
    );
    let parts = || {
        [
            System::state(Ket).at(&a).unwrap(),
            System::state(Ket).at(&b).unwrap(),
        ]
    };
    assert!(
        System::state_tensor(parts())
            .unwrap()
            .space()
            .is_compatible_with(&space)
    );
    assert!(
        System::try_state_tensor(parts().map(Ok))
            .unwrap()
            .space()
            .is_compatible_with(&space)
    );
    let sum = System::state_sum(&space, [tensor.clone(), tensor.clone()]).unwrap();
    assert!(sum.space().shares_node_with(&space));
    let sum = System::try_state_sum(&space, [Ok(tensor)]).unwrap();
    assert!(sum.space().shares_node_with(&space));
    let op = System::operator_product_at([(&a, Op), (&a, Op), (&b, Op)]).unwrap();
    assert!(op.space().is_none());
    for constructed in [
        System::operator_sum([op.clone(), op.clone()]),
        System::try_operator_sum([Ok(op.clone()), Ok(op.clone())]),
        System::operator_product([op.clone(), op.clone()]),
        System::try_operator_product([Ok(op.clone()), Ok(op)]),
    ] {
        assert!(
            constructed
                .unwrap()
                .bind(&space)
                .unwrap()
                .space()
                .unwrap()
                .shares_node_with(&space)
        );
    }
}

#[test]
fn empty_system_builders_keep_the_existing_space_contracts() {
    let space = HilbertSpace::factor(0);
    for zero in [
        System::state_sum(&space, []),
        System::try_state_sum(&space, []),
    ] {
        assert!(zero.unwrap().space().shares_node_with(&space));
    }
    for unit in [
        System::state_tensor::<u8>([]),
        System::try_state_tensor::<u8>([]),
        System::state_tensor_at::<u8>([]),
    ] {
        assert_eq!(unit.unwrap().space().indices().count(), 0);
    }
    for op in [
        System::operator_sum::<u8>([]),
        System::try_operator_sum::<u8>([]),
        System::operator_product::<u8>([]),
        System::try_operator_product::<u8>([]),
        System::operator_product_at::<u8>([]),
    ] {
        assert!(op.unwrap().space().is_none());
    }
}

#[test]
fn placed_builders_validate_leaves_and_reject_only_state_duplicates() {
    use qmb_core::{ExpressionError as E, SpaceError};
    let a = HilbertSpace::factor(0);
    assert_eq!(
        System::state_tensor_at([(&a, Ket), (&a, Ket)]).err(),
        Some(E::Space(SpaceError::DuplicateIndex))
    );
    assert!(System::operator_product_at([(&a, Op), (&a, Op)]).is_ok());
    for nonleaf in [
        HilbertSpace::unit(),
        HilbertSpace::tensor([a.clone()]).unwrap(),
    ] {
        assert_eq!(
            System::state_tensor_at([(&nonleaf, Ket)]).err(),
            Some(E::NotAFactor)
        );
        assert_eq!(
            System::operator_product_at([(&nonleaf, Op)]).err(),
            Some(E::NotAFactor)
        );
    }
    let unit = HilbertSpace::unit();
    // Placement errors take precedence over duplicates among earlier entries.
    assert_eq!(
        System::state_tensor_at([(&a, Ket), (&a, Ket), (&unit, Ket)]).err(),
        Some(E::NotAFactor)
    );
    for operators in [false, true] {
        let pulls = core::cell::Cell::new(0);
        let targets = [&a, &unit, &a];
        if operators {
            assert_eq!(
                System::operator_product_at(targets.into_iter().map(|factor| {
                    pulls.set(pulls.get() + 1);
                    (factor, Op)
                }))
                .err(),
                Some(E::NotAFactor)
            );
        } else {
            assert_eq!(
                System::state_tensor_at(targets.into_iter().map(|factor| {
                    pulls.set(pulls.get() + 1);
                    (factor, Ket)
                }))
                .err(),
                Some(E::NotAFactor)
            );
        }
        assert_eq!(pulls.get(), 2);
    }
}

#[test]
fn builders_share_spaces_without_borrowing_assignment_storage() {
    #[derive(PartialEq, Eq)]
    struct BorrowedIndex<'a>(&'a str);
    let names = [String::from("a"), String::from("b")];
    let (state, operator) = {
        let a = HilbertSpace::factor(BorrowedIndex(&names[0]));
        let b = HilbertSpace::factor(BorrowedIndex(&names[1]));
        let space = a.tensor_product(&b).unwrap();
        let state = System::state_tensor_at([(&a, Ket), (&b, Ket)]).unwrap();
        let op = System::operator_product_at([(&a, Op), (&b, Op)])
            .unwrap()
            .bind(&space)
            .unwrap();
        (state, op)
    };
    assert_eq!(
        state.space().indices().map(|i| i.0).collect::<Vec<_>>(),
        ["a", "b"]
    );
    assert!(operator.space().unwrap().is_compatible_with(state.space()));
}
