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
