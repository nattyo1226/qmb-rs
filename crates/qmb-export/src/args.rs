use qmb_core::HilbertSpace;
use std::marker::PhantomData;

/// Arguments selecting an operator's state type, state basis and optional space.
pub struct MatrixArgs<'b, S, B, A = ()> {
    pub(crate) basis: &'b B,
    pub(crate) ambient: A,
    state: PhantomData<fn() -> S>,
}

impl<'b, S, B> MatrixArgs<'b, S, B> {
    /// Uses a local basis and, for placed operators, an existing binding.
    pub fn new(basis: &'b B) -> Self {
        Self {
            basis,
            ambient: (),
            state: PhantomData,
        }
    }
}

impl<'b, S, B, A> MatrixArgs<'b, S, B, A> {
    /// Selects an explicit ambient space in preference to an old binding.
    pub fn in_space<I>(self, ambient: &HilbertSpace<I>) -> MatrixArgs<'b, S, B, &HilbertSpace<I>> {
        MatrixArgs {
            basis: self.basis,
            ambient,
            state: PhantomData,
        }
    }
}

pub(crate) trait Ambient<I> {
    fn explicit(&self) -> Option<&HilbertSpace<I>>;
}
impl<I> Ambient<I> for () {
    fn explicit(&self) -> Option<&HilbertSpace<I>> {
        None
    }
}
impl<I> Ambient<I> for &HilbertSpace<I> {
    fn explicit(&self) -> Option<&HilbertSpace<I>> {
        Some(self)
    }
}
