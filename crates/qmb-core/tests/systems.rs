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
