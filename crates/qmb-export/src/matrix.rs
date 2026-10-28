use crate::{
    Export, ExportError, Exported, MatrixArgs, MatrixMetadata,
    args::Ambient,
    backend::{MatrixBuilder, MatrixStorage},
    coordinates::owned_entries,
};
use num_traits::{One, Zero};
use qmb_core::{
    Act, Adjoint, Basis, BasisError, Conjugate, Dimension, ExpressionError, HilbertSpace,
    LocalOperator, LocalState, Operator,
};
use std::{marker::PhantomData, ops::Mul};

impl<'b, S, O, C, B, D> Export<Exported<D, MatrixMetadata<'b, S, B>>, MatrixArgs<'b, S, B>>
    for &LocalOperator<O, C>
where
    C: Clone + Zero + Mul<Output = C> + Conjugate<Output = C>,
    S: Dimension + Act<C, Operator = O, Term = S>,
    O: Adjoint<C, Term = O>,
    B: Basis<S, C>,
    D: MatrixStorage<C>,
{
    type Error = ExportError;
    fn export(
        self,
        args: MatrixArgs<'b, S, B>,
    ) -> Result<Exported<D, MatrixMetadata<'b, S, B>>, ExportError> {
        let dimension = S::dim()?.get();
        let mut builder = D::builder(dimension)?;
        let mut count = 0;
        for primitive in args.basis.elements() {
            if count >= dimension || args.basis.index(&primitive) != Some(count) {
                return Err(BasisError::InvalidEnumeration.into());
            }
            let column = LocalState::primitive(primitive)
                .act(self)
                .coordinates_in(args.basis)?;
            builder.column(count, owned_entries(column))?;
            count += 1;
        }
        if count != dimension {
            return Err(BasisError::InvalidEnumeration.into());
        }
        Ok(Exported {
            data: builder.finish()?,
            metadata: MatrixMetadata {
                basis: args.basis,
                space: (),
                dimension,
                factor_count: 1,
                local_dimension: Some(dimension),
                state: PhantomData,
            },
        })
    }
}

impl<'b, I: Eq, S, O, C, B, A, D>
    Export<Exported<D, MatrixMetadata<'b, S, B, HilbertSpace<I>>>, MatrixArgs<'b, S, B, A>>
    for &Operator<I, O, C>
where
    C: Clone + Zero + One + Conjugate<Output = C>,
    S: Dimension + Act<C, Operator = O, Term = S>,
    O: Adjoint<C, Term = O>,
    B: Basis<S, C>,
    A: Ambient<I>,
    D: MatrixStorage<C>,
{
    type Error = ExportError;
    fn export(
        self,
        args: MatrixArgs<'b, S, B, A>,
    ) -> Result<Exported<D, MatrixMetadata<'b, S, B, HilbertSpace<I>>>, ExportError> {
        let space = args
            .ambient
            .explicit()
            .or_else(|| self.space())
            .ok_or(ExpressionError::UnboundOperator)?;
        // Validate every reference even for a zero operator, before allocating output.
        let bound = self.bind(space)?;
        let dimension = space.dimension::<S>()?;
        let factor_count = space.indices().count();
        let local_dimension = if factor_count == 0 {
            None
        } else {
            Some(S::dim()?.get())
        };
        let mut builder = D::builder(dimension)?;
        for (col, input) in space.basis_states::<S, O, C, _>(args.basis)?.enumerate() {
            let column = input?.act(&bound)?.coordinates_in(args.basis)?;
            builder.column(col, owned_entries(column))?;
        }
        Ok(Exported {
            data: builder.finish()?,
            metadata: MatrixMetadata {
                basis: args.basis,
                space: space.clone(),
                dimension,
                factor_count,
                local_dimension,
                state: PhantomData,
            },
        })
    }
}
