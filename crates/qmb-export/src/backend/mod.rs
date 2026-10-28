//! Internal storage construction contracts and allocation checks.

use crate::ExportError;
use std::alloc::Layout;

#[cfg(feature = "ndarray")]
mod ndarray;
#[cfg(feature = "sprs")]
mod sprs;
pub(crate) mod vec;

// Storage adapters share coordinate extraction, metadata and precision handling.
pub(crate) trait VectorStorage<C>: Sized {
    fn build(
        dimension: usize,
        entries: impl ExactSizeIterator<Item = Result<(usize, C), ExportError>>,
    ) -> Result<Self, ExportError>;
}

pub(crate) trait MatrixStorage<C>: Sized {
    type Builder: MatrixBuilder<C, Output = Self>;
    fn builder(dimension: usize) -> Result<Self::Builder, ExportError>;
}

pub(crate) trait MatrixBuilder<C> {
    type Output;
    fn column(
        &mut self,
        col: usize,
        entries: impl ExactSizeIterator<Item = Result<(usize, C), ExportError>>,
    ) -> Result<(), ExportError>;
    fn finish(self) -> Result<Self::Output, ExportError>;
}

pub(crate) fn reserve<T>(values: &mut Vec<T>, additional: usize) -> Result<(), ExportError> {
    let length = values
        .len()
        .checked_add(additional)
        .ok_or(ExportError::SizeOverflow)?;
    Layout::array::<T>(length).map_err(|_| ExportError::SizeOverflow)?;
    values.try_reserve(additional)?;
    Ok(())
}
