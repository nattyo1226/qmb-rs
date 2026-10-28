/// Boundary condition shared by geometries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// Connect consecutive sites only.
    Open,
    /// Also connect the last site to the first, without loops or duplicate edges.
    Periodic,
}
