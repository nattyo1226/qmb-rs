use qmb_core::HilbertSpace;
use std::marker::PhantomData;

/// Numerical data paired with the metadata needed to interpret its coordinates.
pub struct Exported<D, M> {
    pub(crate) data: D,
    pub(crate) metadata: M,
}

impl<D, M> Exported<D, M> {
    /// Borrows the numerical container.
    pub fn data(&self) -> &D {
        &self.data
    }
    /// Borrows its basis, shape and factor-space interpretation.
    pub fn metadata(&self) -> &M {
        &self.metadata
    }
    /// Takes the numerical container, discarding the metadata.
    pub fn into_data(self) -> D {
        self.data
    }
    /// Takes both the numerical container and its interpretation.
    pub fn into_parts(self) -> (D, M) {
        (self.data, self.metadata)
    }
}

/// Dense matrix memory layout, separate from tensor-factor basis order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layout {
    /// Element `(row, col)` is at `row * ncols + col`.
    RowMajor,
}

/// Borrowed local basis, shared space and shape of a coordinate vector.
///
/// `H` is `()` for local expressions, and `HilbertSpace<I>` for placed ones.
/// `P` records a source state's primitive type; direct Coordinates exports use
/// `P = ()` because Coordinates does not retain a primitive type.
pub struct VectorMetadata<'b, B, H = (), P = ()> {
    pub(crate) basis: &'b B,
    pub(crate) space: H,
    pub(crate) dimension: usize,
    pub(crate) factor_count: usize,
    pub(crate) local_dimension: Option<usize>,
    pub(crate) primitive: PhantomData<fn() -> P>,
}

impl<'b, B, H, P> VectorMetadata<'b, B, H, P> {
    /// Borrows the uniform local basis.
    pub fn basis(&self) -> &'b B {
        self.basis
    }
    /// Borrows the selected codomain space, or `()` for a local expression.
    pub fn space(&self) -> &H {
        &self.space
    }
    /// Borrows the vector's codomain interpretation.
    pub fn codomain(&self) -> &H {
        &self.space
    }
    /// Returns the one-dimensional container shape.
    pub fn shape(&self) -> [usize; 1] {
        [self.dimension]
    }
    /// Returns the ket's scalar-domain dimension.
    pub fn domain_dimension(&self) -> usize {
        1
    }
    /// Returns the number of basic factors in tensor order.
    pub fn factor_count(&self) -> usize {
        self.factor_count
    }
    /// Returns the local basis size, or None for the zero-factor unit space.
    pub fn local_dimension(&self) -> Option<usize> {
        self.local_dimension
    }
    pub(crate) fn with_primitive<Q>(self) -> VectorMetadata<'b, B, H, Q> {
        VectorMetadata {
            basis: self.basis,
            space: self.space,
            dimension: self.dimension,
            factor_count: self.factor_count,
            local_dimension: self.local_dimension,
            primitive: PhantomData,
        }
    }
}

impl<B, I, P> VectorMetadata<'_, B, HilbertSpace<I>, P> {
    /// Returns the zero-factor scalar domain of a ket.
    pub fn domain(&self) -> HilbertSpace<I> {
        HilbertSpace::unit()
    }
}

/// Shape, borrowed state basis and selected action space of a matrix.
///
/// Input and output use the same state basis and ordered factor space.
pub struct MatrixMetadata<'b, S, B, H = ()> {
    pub(crate) basis: &'b B,
    pub(crate) space: H,
    pub(crate) dimension: usize,
    pub(crate) factor_count: usize,
    pub(crate) local_dimension: Option<usize>,
    pub(crate) state: PhantomData<fn() -> S>,
}

impl<'b, S, B, H> MatrixMetadata<'b, S, B, H> {
    /// Borrows the common local input and output state basis.
    pub fn basis(&self) -> &'b B {
        self.basis
    }
    /// Borrows the selected action space, or `()` for a local operator.
    pub fn space(&self) -> &H {
        &self.space
    }
    /// Borrows the selected input space.
    pub fn domain(&self) -> &H {
        &self.space
    }
    /// Borrows the selected output space.
    pub fn codomain(&self) -> &H {
        &self.space
    }
    /// Returns the number of rows and columns.
    pub fn shape(&self) -> [usize; 2] {
        [self.dimension; 2]
    }
    /// Returns the dense export's memory layout.
    pub fn layout(&self) -> Layout {
        Layout::RowMajor
    }
    /// Returns the number of basic factors in tensor order.
    pub fn factor_count(&self) -> usize {
        self.factor_count
    }
    /// Returns the local state-basis size, or None on the unit space.
    pub fn local_dimension(&self) -> Option<usize> {
        self.local_dimension
    }
}

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
