//! Shared conversion of sparse coordinates to flat indices.

use crate::ExportError;
use qmb_core::{Coordinates, DimensionOverflow};

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
