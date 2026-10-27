use crate::{
    ExportError,
    error::reserve,
    matrix::{MatrixBuilder, MatrixStorage},
    vector::VectorStorage,
};
use sprs::{CsMat, CsVec};

impl<C> VectorStorage<C> for CsVec<C> {
    fn build(
        dimension: usize,
        entries: impl ExactSizeIterator<Item = Result<(usize, C), ExportError>>,
    ) -> Result<Self, ExportError> {
        let mut indices = Vec::new();
        let mut values = Vec::new();
        reserve(&mut indices, entries.len())?;
        reserve(&mut values, entries.len())?;
        for entry in entries {
            let (index, value) = entry?;
            indices.push(index);
            values.push(value);
        }
        Self::try_new(dimension, indices, values).map_err(|(_, _, error)| error.into())
    }
}

pub(crate) struct SparseBuilder<C> {
    dimension: usize,
    indptr: Vec<usize>,
    indices: Vec<usize>,
    values: Vec<C>,
}

impl<C> MatrixStorage<C> for CsMat<C> {
    type Builder = SparseBuilder<C>;
    fn builder(dimension: usize) -> Result<Self::Builder, ExportError> {
        let pointers = dimension.checked_add(1).ok_or(ExportError::SizeOverflow)?;
        let mut indptr = Vec::new();
        reserve(&mut indptr, pointers)?;
        indptr.push(0);
        Ok(SparseBuilder {
            dimension,
            indptr,
            indices: Vec::new(),
            values: Vec::new(),
        })
    }
}

impl<C> MatrixBuilder<C> for SparseBuilder<C> {
    type Output = CsMat<C>;
    fn column(
        &mut self,
        col: usize,
        entries: impl ExactSizeIterator<Item = Result<(usize, C), ExportError>>,
    ) -> Result<(), ExportError> {
        if col >= self.dimension || col != self.indptr.len() - 1 {
            return Err(ExportError::InvalidShape);
        }
        reserve(&mut self.indices, entries.len())?;
        reserve(&mut self.values, entries.len())?;
        for entry in entries {
            let (row, value) = entry?;
            self.indices.push(row);
            self.values.push(value);
        }
        self.indptr.push(self.values.len());
        Ok(())
    }
    fn finish(self) -> Result<CsMat<C>, ExportError> {
        CsMat::try_new_csc(
            (self.dimension, self.dimension),
            self.indptr,
            self.indices,
            self.values,
        )
        .map_err(|(_, _, _, error)| error.into())
    }
}
