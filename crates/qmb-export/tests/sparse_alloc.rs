#![cfg(feature = "sprs")]

use num_complex::Complex;
use qmb_core::{HilbertSpace, LocalState, Operator, State};
use qmb_export::{Export, ExportError, Exported, MatrixArgs, MatrixMetadata, VectorMetadata};
use qmb_spin_half::{SpinHalfOperator as O, SpinHalfState as S, SzBasis};
use std::{
    alloc::{GlobalAlloc, Layout, System},
    cell::Cell,
    ptr,
};

struct BoundedAllocator;
thread_local! { static LIMIT: Cell<Option<usize>>=const {Cell::new(None)}; }
fn allowed(size: usize) -> bool {
    LIMIT
        .try_with(|limit| limit.get().is_none_or(|max| size <= max))
        .unwrap_or(true)
}
// Forwards ordinary allocations to System. During one test-thread check,
// oversized requests return null so an accidental dense allocation fails safely.
unsafe impl GlobalAlloc for BoundedAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if allowed(layout.size()) {
            unsafe { System.alloc(layout) }
        } else {
            ptr::null_mut()
        }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        if allowed(layout.size()) {
            unsafe { System.alloc_zeroed(layout) }
        } else {
            ptr::null_mut()
        }
    }
    unsafe fn realloc(&self, p: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        if allowed(size) {
            unsafe { System.realloc(p, layout, size) }
        } else {
            ptr::null_mut()
        }
    }
    unsafe fn dealloc(&self, p: *mut u8, layout: Layout) {
        unsafe { System.dealloc(p, layout) }
    }
}
#[global_allocator]
static ALLOCATOR: BoundedAllocator = BoundedAllocator;
fn bounded<R>(f: impl FnOnce() -> R) -> R {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            LIMIT.with(|v| v.set(None));
        }
    }
    LIMIT.with(|v| v.set(Some(64 * 1024)));
    let _guard = Reset;
    f()
}

#[test]
fn sparse_exports_need_no_dense_intermediate_and_dense_allocation_failure_is_reported() {
    let basis = SzBasis;
    let space = HilbertSpace::tensor((0..40u8).map(HilbertSpace::factor)).unwrap();
    let state = State::tensor(space.indices().map(|&i| {
        LocalState::<S, O, Complex<f64>>::primitive(S::Up)
            .at(&HilbertSpace::factor(i))
            .unwrap()
    }))
    .unwrap();
    let sparse: Exported<qmb_export::sprs::CsVec<_>, VectorMetadata<'_, _, _, S>> =
        bounded(|| (&state).export(&basis)).unwrap();
    assert_eq!(sparse.data().dim(), 1usize << 40);
    assert_eq!(sparse.data().nnz(), 1);
    assert_eq!(sparse.data().indices(), &[0]);
    let dense: Result<Exported<Vec<_>, VectorMetadata<'_, _, _, S>>, _> =
        bounded(|| (&state).export(&basis));
    assert!(matches!(dense, Err(ExportError::Allocation(_))));
    let space = HilbertSpace::tensor((0..9u8).map(HilbertSpace::factor)).unwrap();
    let identity = Operator::<u8, O, Complex<f64>>::identity();
    let matrix: Exported<qmb_export::sprs::CsMat<_>, MatrixMetadata<'_, S, _, _>> =
        bounded(|| (&identity).export(MatrixArgs::new(&basis).in_space(&space))).unwrap();
    assert_eq!(matrix.data().shape(), (512, 512));
    assert_eq!(matrix.data().nnz(), 512);
    assert!(matrix.data().is_csc());
}
