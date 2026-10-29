//! The same input-error priority applies to every fallible collection entry.
use core::cell::Cell;
use qmb_core::{
    ExpressionError as E, HilbertSpace as H, LocalOperator, LocalState, Operator, State,
};

type System = qmb_core::ExpressionFactory<(), (), ()>;

// Panics on any next() after the first Err, including an end-of-input probe.
fn failing_at<'a, T: 'a>(
    position: usize,
    pulls: &'a Cell<usize>,
    mut make: impl FnMut() -> T + 'a,
) -> impl Iterator<Item = Result<T, E>> + 'a {
    (0..).map(move |i| {
        assert!(i <= position, "input consumed after the first error");
        pulls.set(pulls.get() + 1);
        if i == position {
            Err(E::NotAFactor)
        } else {
            Ok(make())
        }
    })
}

#[test]
fn every_try_collection_stops_at_the_first_input_error() {
    let a = H::factor(0u8);
    for position in 0..=3 {
        let pulls = Cell::new(0);
        macro_rules! check {
            ($call:expr) => {{
                pulls.set(0);
                assert_eq!($call.err(), Some(E::NotAFactor));
                assert_eq!(pulls.get(), position + 1);
            }};
        }
        check!(LocalState::<(), (), ()>::try_sum(failing_at(
            position,
            &pulls,
            || LocalState::primitive(())
        )));
        check!(LocalOperator::<(), ()>::try_sum(failing_at(
            position,
            &pulls,
            || LocalOperator::primitive(())
        )));
        check!(State::<_, (), (), ()>::try_sum(
            &a,
            failing_at(position, &pulls, || System::state(()).at(&a).unwrap())
        ));
        check!(State::<_, (), (), ()>::try_tensor_in(
            &a,
            failing_at(position, &pulls, || System::state(()).at(&a).unwrap())
        ));
        check!(Operator::<_, (), ()>::try_sum(failing_at(
            position,
            &pulls,
            || System::operator(()).at(&a).unwrap()
        )));
        check!(Operator::<_, (), ()>::try_product(failing_at(
            position,
            &pulls,
            || System::operator(()).at(&a).unwrap()
        )));
        check!(System::try_state_sum(
            &a,
            failing_at(position, &pulls, || System::state(()).at(&a).unwrap())
        ));
        check!(System::try_tensor_in(
            &a,
            failing_at(position, &pulls, || System::state(()).at(&a).unwrap())
        ));
        check!(System::try_operator_sum(failing_at(
            position,
            &pulls,
            || System::operator(()).at(&a).unwrap()
        )));
        check!(System::try_operator_product(failing_at(
            position,
            &pulls,
            || System::operator(()).at(&a).unwrap()
        )));
    }
}

#[test]
fn input_errors_precede_space_and_duplicate_errors() {
    let a = H::factor(0u8);
    let b = H::factor(1);
    let state = || System::state(()).at(&a).unwrap();
    assert_eq!(
        System::try_state_sum(&b, [Ok(state()), Err(E::MissingIndex)]).err(),
        Some(E::MissingIndex)
    );
    assert_eq!(
        System::try_tensor_in(&a, [Ok(state()), Ok(state()), Err(E::MissingIndex)]).err(),
        Some(E::MissingIndex)
    );
    assert_eq!(
        System::try_state_sum(&b, [Ok(state())]).err(),
        Some(E::IncompatibleSpaces)
    );
    assert_eq!(
        System::try_tensor_in(&a, [Ok(state()), Ok(state())]).err(),
        Some(E::IncompatibleSpaces)
    );
    let left = || Operator::<_, (), ()>::identity().bind(&a).unwrap();
    let right = || Operator::<_, (), ()>::identity().bind(&b).unwrap();
    assert_eq!(
        System::try_operator_sum([Ok(left()), Ok(right()), Err(E::MissingIndex)]).err(),
        Some(E::MissingIndex)
    );
    assert_eq!(
        System::try_operator_product([Ok(left()), Ok(right()), Err(E::MissingIndex)]).err(),
        Some(E::MissingIndex)
    );
    assert_eq!(
        System::try_operator_sum([Ok(left()), Ok(right())]).err(),
        Some(E::IncompatibleSpaces)
    );
    assert_eq!(
        System::try_operator_product([Ok(left()), Ok(right())]).err(),
        Some(E::IncompatibleSpaces)
    );
    // A later unbound operand changes the whole binding decision.
    assert!(
        System::try_operator_sum([Ok(left()), Ok(right()), Ok(Operator::identity())])
            .unwrap()
            .space()
            .is_none()
    );
    assert!(
        System::try_operator_product([Ok(left()), Ok(right()), Ok(Operator::identity())])
            .unwrap()
            .space()
            .is_none()
    );
    let same = System::try_operator_sum([Ok(left()), Ok(left())]).unwrap();
    assert!(same.space().unwrap().shares_node_with(&a));
    let same = System::try_operator_product([Ok(left()), Ok(left())]).unwrap();
    assert!(same.space().unwrap().shares_node_with(&a));
}

#[test]
fn eager_inputs_are_not_rolled_back_and_successful_inputs_are_fully_consumed() {
    let a = H::factor(0u8);
    let b = H::factor(1);
    let pulls = Cell::new(0);
    let eager: Vec<_> = (0..3)
        .map(|i| {
            pulls.set(pulls.get() + 1);
            if i == 1 {
                Err(E::NotAFactor)
            } else {
                Ok(System::state(()).at(&a).unwrap())
            }
        })
        .collect();
    assert_eq!(System::try_state_sum(&a, eager).err(), Some(E::NotAFactor));
    assert_eq!(pulls.get(), 3);
    pulls.set(0);
    let all_ok = (0..3).map(|_| {
        pulls.set(pulls.get() + 1);
        Ok(System::state(()).at(&a).unwrap())
    });
    assert_eq!(
        System::try_state_sum(&b, all_ok).err(),
        Some(E::IncompatibleSpaces)
    );
    assert_eq!(pulls.get(), 3);
}
