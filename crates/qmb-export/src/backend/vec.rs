//! Dense Vec storage.

use super::{MatrixBuilder, MatrixStorage, VectorStorage, reserve};
use crate::ExportError;
use num_traits::Zero;

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

pub(crate) struct DenseBuilder<C> {
    pub(crate) values: Vec<C>,
    pub(crate) dimension: usize,
}

impl<C: Zero> MatrixStorage<C> for Vec<C> {
    type Builder = DenseBuilder<C>;
    fn builder(dimension: usize) -> Result<Self::Builder, ExportError> {
        let count = dimension
            .checked_mul(dimension)
            .ok_or(ExportError::SizeOverflow)?;
        let mut values = Vec::new();
        reserve(&mut values, count)?;
        values.resize_with(count, C::zero);
        Ok(DenseBuilder { values, dimension })
    }
}

impl<C> MatrixBuilder<C> for DenseBuilder<C> {
    type Output = Vec<C>;
    fn column(
        &mut self,
        col: usize,
        entries: impl ExactSizeIterator<Item = Result<(usize, C), ExportError>>,
    ) -> Result<(), ExportError> {
        if col >= self.dimension {
            return Err(ExportError::InvalidShape);
        }
        for entry in entries {
            let (row, c) = entry?;
            if row >= self.dimension {
                return Err(ExportError::InvalidShape);
            }
            self.values[row * self.dimension + col] = c;
        }
        Ok(())
    }
    fn finish(self) -> Result<Vec<C>, ExportError> {
        Ok(self.values)
    }
}
