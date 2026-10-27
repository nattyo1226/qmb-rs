use crate::{Export, ExportError, Exported, VectorMetadata, error::reserve};
use num_traits::{One, Zero};
use qmb_core::{
    Act, Adjoint, Basis, Conjugate, Coordinates, DimensionOverflow, HilbertSpace, LocalState,
    SpanDimension, State,
};
use std::{marker::PhantomData, ops::Mul};

// Storage adapters share coordinate extraction, metadata and precision handling.
pub(crate) trait VectorStorage<C>: Sized {
    fn build(
        dimension: usize,
        entries: impl ExactSizeIterator<Item = Result<(usize, C), ExportError>>,
    ) -> Result<Self, ExportError>;
}

impl<C: Zero> VectorStorage<C> for Vec<C> {
    fn build(
        dimension: usize,
        entries: impl ExactSizeIterator<Item = Result<(usize, C), ExportError>>,
    ) -> Result<Self, ExportError> {
        let mut values = Vec::new();
        reserve(&mut values, dimension)?;
        values.resize_with(dimension, C::zero);
        for entry in entries {
            let (index, value) = entry?;
            *values.get_mut(index).ok_or(ExportError::InvalidShape)? = value;
        }
        Ok(values)
    }
}

pub(crate) fn flat_index(digits: &[usize], local: Option<usize>) -> Result<usize, ExportError> {
    let mut index = 0usize;
    for &digit in digits {
        let radix = local.ok_or(ExportError::InvalidShape)?;
        if digit >= radix {
            return Err(ExportError::InvalidShape);
        }
        index = index
            .checked_mul(radix)
            .and_then(|n| n.checked_add(digit))
            .ok_or(DimensionOverflow)?;
    }
    Ok(index)
}

pub(crate) fn owned_entries<B, C, H>(
    coordinates: Coordinates<'_, B, C, H>,
) -> impl ExactSizeIterator<Item = Result<(usize, C), ExportError>> {
    let local = coordinates.local_dimension();
    coordinates
        .into_entries()
        .map(move |(digits, c)| Ok((flat_index(&digits, local)?, c)))
}

impl<'b, B, C, H: Clone, D: VectorStorage<C>> Export<Exported<D, VectorMetadata<'b, B, H>>>
    for Coordinates<'b, B, C, H>
{
    type Error = ExportError;
    fn export(self, (): ()) -> Result<Exported<D, VectorMetadata<'b, B, H>>, ExportError> {
        let metadata = VectorMetadata {
            basis: self.basis(),
            space: self.space().clone(),
            dimension: self.dimension()?,
            factor_count: self.factor_count(),
            local_dimension: self.local_dimension(),
            primitive: PhantomData,
        };
        let data = D::build(metadata.dimension, owned_entries(self))?;
        Ok(Exported { data, metadata })
    }
}

impl<'b, B, C: Clone, H: Clone, D: VectorStorage<C>> Export<Exported<D, VectorMetadata<'b, B, H>>>
    for &Coordinates<'b, B, C, H>
{
    type Error = ExportError;
    fn export(self, (): ()) -> Result<Exported<D, VectorMetadata<'b, B, H>>, ExportError> {
        let metadata = VectorMetadata {
            basis: self.basis(),
            space: self.space().clone(),
            dimension: self.dimension()?,
            factor_count: self.factor_count(),
            local_dimension: self.local_dimension(),
            primitive: PhantomData,
        };
        let data = D::build(
            metadata.dimension,
            self.flat_entries().map(|entry| {
                let (index, c) = entry?;
                Ok((index, c.clone()))
            }),
        )?;
        Ok(Exported { data, metadata })
    }
}

impl<'b, S, O, C, B, D> Export<Exported<D, VectorMetadata<'b, B, (), S>>, &'b B>
    for &LocalState<S, O, C>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
    S: SpanDimension + Act<C, Operator = O, Term = S>,
    O: Adjoint<C, Term = O>,
    B: Basis<S, C>,
    D: VectorStorage<C>,
{
    type Error = ExportError;
    fn export(
        self,
        basis: &'b B,
    ) -> Result<Exported<D, VectorMetadata<'b, B, (), S>>, ExportError> {
        let result: Exported<D, _> = self.coordinates_in(basis)?.export(())?;
        Ok(Exported {
            data: result.data,
            metadata: result.metadata.with_primitive(),
        })
    }
}

impl<'b, I: Eq, S, O, C, B, D> Export<Exported<D, VectorMetadata<'b, B, HilbertSpace<I>, S>>, &'b B>
    for &State<I, S, O, C>
where
    C: Clone + Zero + One + Conjugate<Output = C>,
    S: SpanDimension + Act<C, Operator = O, Term = S>,
    O: Adjoint<C, Term = O>,
    B: Basis<S, C>,
    D: VectorStorage<C>,
{
    type Error = ExportError;
    fn export(
        self,
        basis: &'b B,
    ) -> Result<Exported<D, VectorMetadata<'b, B, HilbertSpace<I>, S>>, ExportError> {
        let result: Exported<D, _> = self.coordinates_in(basis)?.export(())?;
        Ok(Exported {
            data: result.data,
            metadata: result.metadata.with_primitive(),
        })
    }
}
