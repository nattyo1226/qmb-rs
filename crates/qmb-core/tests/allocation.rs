//! Check the promised allocation behavior independently of the representation.

use std::{
    alloc::{GlobalAlloc, Layout, System},
    cell::Cell,
};

use qmb_core::{Complex, LinearCombination, Term};

struct CountingAllocator;
thread_local! {
    static ALLOCATIONS: Cell<Option<usize>> = const { Cell::new(None) };
}

fn count_allocation() {
    let _ = ALLOCATIONS.try_with(|count| {
        if let Some(n) = count.get() {
            count.set(Some(n + 1));
        }
    });
}

// The wrapper forwards each allocation unchanged to System; the counter only
// observes this thread while a check is active and never allocates itself.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        count_allocation();
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        count_allocation();
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        count_allocation();
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

fn without_allocation(f: impl FnOnce()) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            ALLOCATIONS.with(|count| count.set(None));
        }
    }
    ALLOCATIONS.with(|count| count.set(Some(0)));
    let _reset = Reset;
    f();
    assert_eq!(ALLOCATIONS.with(Cell::get), Some(0));
}

#[test]
fn empty_and_single_term_storage_and_iteration_do_not_allocate() {
    without_allocation(|| {
        let zero = LinearCombination::<u8, f64>::zero();
        assert_eq!(zero.terms().next(), None);
        assert_eq!(zero.into_terms().next(), None);
        let empty = LinearCombination::<u8, f64>::from_terms([]);
        assert_eq!(empty.into_iter().len(), 0);
        let single = LinearCombination::from_terms([Term {
            coefficient: Complex::new(2.0, 0.0),
            primitive: 1,
        }]);
        let single = single + LinearCombination::zero();
        assert_eq!(single.terms().len(), 1);
        assert_eq!(single.into_terms().next().unwrap().primitive, 1);
        let single = LinearCombination::<_, f64>::primitive(2).scaled(Complex::new(0.0, 0.0));
        assert_eq!(single.into_iter().next_back().unwrap().primitive, 2);
    });
}

#[test]
fn owning_iteration_reuses_the_existing_multiple_term_storage() {
    let sum = LinearCombination::from_terms((0..5).map(|p| Term {
        coefficient: Complex::new(1.0, 0.0),
        primitive: p,
    }));
    without_allocation(|| {
        let mut iter = sum.into_terms();
        assert_eq!(iter.next().unwrap().primitive, 0);
        assert_eq!(iter.next_back().unwrap().primitive, 4);
        assert_eq!(iter.len(), 3);
    });
}
