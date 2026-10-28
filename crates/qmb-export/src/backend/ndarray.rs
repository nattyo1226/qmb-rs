use crate::{
    Export, ExportError, Exported, MatrixMetadata, VectorMetadata,
    backend::{MatrixBuilder, MatrixStorage, VectorStorage, vec::DenseBuilder},
};
use ndarray::{Array1, Array2};
use num_traits::Zero;

fn vector<C>(values: Vec<C>, dimension: usize) -> Result<Array1<C>, ExportError> {
    if values.len() != dimension {
        return Err(ExportError::InvalidShape);
    }
    Ok(Array1::from_vec(values))
}

fn matrix<C>(values: Vec<C>, dimension: usize) -> Result<Array2<C>, ExportError> {
    Ok(Array2::from_shape_vec((dimension, dimension), values)?)
}

impl<'b, B, C, H, P> Export<Exported<Array1<C>, VectorMetadata<'b, B, H, P>>>
    for Exported<Vec<C>, VectorMetadata<'b, B, H, P>>
{
    type Error = ExportError;
    fn export(
        self,
        (): (),
    ) -> Result<Exported<Array1<C>, VectorMetadata<'b, B, H, P>>, ExportError> {
        Ok(Exported {
            data: vector(self.data, self.metadata.dimension)?,
            metadata: self.metadata,
        })
    }
}

impl<'b, S, B, C, H> Export<Exported<Array2<C>, MatrixMetadata<'b, S, B, H>>>
    for Exported<Vec<C>, MatrixMetadata<'b, S, B, H>>
{
    type Error = ExportError;
    fn export(
        self,
        (): (),
    ) -> Result<Exported<Array2<C>, MatrixMetadata<'b, S, B, H>>, ExportError> {
        let d = self.metadata.dimension;
        let data = matrix(self.data, d)?;
        Ok(Exported {
            data,
            metadata: self.metadata,
        })
    }
}

impl<C: Zero> VectorStorage<C> for Array1<C> {
    fn build(
        dimension: usize,
        entries: impl ExactSizeIterator<Item = Result<(usize, C), ExportError>>,
    ) -> Result<Self, ExportError> {
        vector(Vec::build(dimension, entries)?, dimension)
    }
}

pub(crate) struct ArrayBuilder<C>(DenseBuilder<C>);
impl<C: Zero> MatrixStorage<C> for Array2<C> {
    type Builder = ArrayBuilder<C>;
    fn builder(dimension: usize) -> Result<Self::Builder, ExportError> {
        Ok(ArrayBuilder(Vec::builder(dimension)?))
    }
}
impl<C> MatrixBuilder<C> for ArrayBuilder<C> {
    type Output = Array2<C>;
    fn column(
        &mut self,
        col: usize,
        entries: impl ExactSizeIterator<Item = Result<(usize, C), ExportError>>,
    ) -> Result<(), ExportError> {
        self.0.column(col, entries)
    }
    fn finish(self) -> Result<Array2<C>, ExportError> {
        let d = self.0.dimension;
        matrix(self.0.finish()?, d)
    }
}
