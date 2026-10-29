use qmb_core::{
    HilbertSpace, LocalOperator, LocalOperatorOf, LocalState, LocalStateOf, OperatorOf,
    QuantumSystem, StateOf,
};

// These types deliberately provide no capabilities (including Clone).
struct Ket;
struct Op;
struct Coeff;
type System = qmb_core::ExpressionFactory<Ket, Op, Coeff>;
struct Equivalent;
impl QuantumSystem for Equivalent {
    type State = Ket;
    type Operator = Op;
    type Coeff = Coeff;
}

#[test]
fn independent_construction_and_aliases_need_no_evaluation_capabilities() {
    fn state<Q: QuantumSystem>(primitive: Q::State) -> LocalStateOf<Q> {
        LocalStateOf::<Q>::primitive(primitive)
    }
    let local: LocalState<Ket, Op, Coeff> = state::<System>(Ket);
    let operator: LocalOperator<Op, Coeff> = System::operator(Op);
    let _: LocalOperatorOf<System> = operator.clone();
    let _: LocalStateOf<System> = local.clone();
    let factor = HilbertSpace::factor(0);
    let state: StateOf<System, _> = local.at(&factor).unwrap();
    let op: OperatorOf<System, _> = operator.at(&factor).unwrap();
    // The chosen system is not an expression identity tag.
    let other = LocalStateOf::<Equivalent>::primitive(Ket)
        .at(&factor)
        .unwrap();
    let _: StateOf<Equivalent, _> = (state.clone() + &other).unwrap();
    let _: OperatorOf<Equivalent, _> = (op.clone()
        * &LocalOperatorOf::<Equivalent>::primitive(Op)
            .at(&factor)
            .unwrap())
        .unwrap();
    let _ = local.act(&operator);
}

#[test]
fn all_system_builders_accept_non_clone_primitives_coefficients_and_indices() {
    #[derive(PartialEq, Eq)]
    struct Index(u8);
    let a = HilbertSpace::factor(Index(3));
    let b = HilbertSpace::factor(Index(1));
    let space = a.tensor(&b).unwrap();
    let tensor = System::product_state(&space, [Ket, Ket]).unwrap();
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
        System::tensor_in(&space, parts())
            .unwrap()
            .space()
            .is_compatible_with(&space)
    );
    assert!(
        System::try_tensor_in(&space, parts().map(Ok))
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
    let unit_space = HilbertSpace::<u8>::tensor_all_g([HilbertSpace::unit()]).unwrap();
    for unit in [
        System::tensor_in(&unit_space, []),
        System::try_tensor_in(&unit_space, []),
        System::product_state(&unit_space, []),
    ] {
        assert!(unit.unwrap().space().shares_node_with(&unit_space));
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
fn explicit_builders_validate_counts_and_operators_still_validate_placement() {
    use qmb_core::ExpressionError as E;
    let a = HilbertSpace::factor(0);
    for states in [vec![], vec![Ket, Ket]] {
        assert_eq!(
            System::product_state(&a, states).err(),
            Some(E::IncompatibleSpaces)
        );
    }
    assert_eq!(
        System::tensor_in(
            &a,
            [
                System::state(Ket).at(&a).unwrap(),
                System::state(Ket).at(&a).unwrap()
            ]
        )
        .err(),
        Some(E::IncompatibleSpaces)
    );
    assert!(System::operator_product_at([(&a, Op), (&a, Op)]).is_ok());
    for nonleaf in [
        HilbertSpace::unit(),
        HilbertSpace::tensor_all([a.clone()]).unwrap(),
    ] {
        assert_eq!(
            System::operator_product_at([(&nonleaf, Op)]).err(),
            Some(E::NotAFactor)
        );
    }
    let unit = HilbertSpace::unit();
    let pulls = core::cell::Cell::new(0);
    assert_eq!(
        System::operator_product_at([&a, &unit, &a].into_iter().map(|factor| {
            pulls.set(pulls.get() + 1);
            (factor, Op)
        }))
        .err(),
        Some(E::NotAFactor)
    );
    assert_eq!(pulls.get(), 2);
}

#[test]
fn builders_share_spaces_without_borrowing_assignment_storage() {
    #[derive(PartialEq, Eq)]
    struct BorrowedIndex<'a>(&'a str);
    let names = [String::from("a"), String::from("b")];
    let (state, operator) = {
        let a = HilbertSpace::factor(BorrowedIndex(&names[0]));
        let b = HilbertSpace::factor(BorrowedIndex(&names[1]));
        let space = a.tensor(&b).unwrap();
        let state = System::product_state(&space, [Ket, Ket]).unwrap();
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
