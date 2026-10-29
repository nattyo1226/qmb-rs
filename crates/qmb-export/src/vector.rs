use crate::{
    Export, ExportError, Exported, VectorMetadata, backend::VectorStorage,
    coordinates::into_flat_entries,
};
use num_traits::{One, Zero};
use qmb_core::{
    Act, Adjoint, Basis, Conjugate, Coordinates, Dimension, HilbertSpace, LocalState, State,
};
use std::{marker::PhantomData, ops::Mul};

impl<'b, B, C, H: Clone, D: VectorStorage<C>> Export<Exported<D, VectorMetadata<'b, B, H>>>
    for Coordinates<'b, B, C, H>
{
    type Error = ExportError;
    fn export(self, (): ()) -> Result<Exported<D, VectorMetadata<'b, B, H>>, ExportError> {
        let metadata = VectorMetadata {
            basis: self.basis(),
            space: self.space().clone(),
            dimension: self.dim()?,
            factor_count: self.factor_count(),
            dim_local: self.dim_local(),
            primitive: PhantomData,
        };
        let data = D::build(metadata.dimension, into_flat_entries(self))?;
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
            dimension: self.dim()?,
            factor_count: self.factor_count(),
            dim_local: self.dim_local(),
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
    S: Dimension + Act<C, Operator = O, Term = S>,
    O: Adjoint<C, Term = O>,
    B: Basis<S, C>,
    D: VectorStorage<C>,
{
    type Error = ExportError;
    fn export(
        self,
        basis: &'b B,
    ) -> Result<Exported<D, VectorMetadata<'b, B, (), S>>, ExportError> {
        let result: Exported<D, _> = self.coordinates(basis)?.export(())?;
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
    S: Dimension + Act<C, Operator = O, Term = S>,
    O: Adjoint<C, Term = O> + qmb_core::HomogeneousParity + qmb_core::ParityOperator,
    B: Basis<S, C>,
    D: VectorStorage<C>,
{
    type Error = ExportError;
    fn export(
        self,
        basis: &'b B,
    ) -> Result<Exported<D, VectorMetadata<'b, B, HilbertSpace<I>, S>>, ExportError> {
        let result: Exported<D, _> = self.coordinates(basis)?.export(())?;
        Ok(Exported {
            data: result.data,
            metadata: result.metadata.with_primitive(),
        })
    }
}
