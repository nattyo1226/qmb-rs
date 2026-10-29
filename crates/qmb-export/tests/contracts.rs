#[cfg(feature = "sprs")]
use qmb_core::Coordinates;
use qmb_core::{
    Act, Adjoint, Basis, BasisError, Dimension, DimensionOverflow, ExpressionError, HilbertSpace,
    LinearCombination, LocalOperator, LocalState, Operator, State,
};
use qmb_export::{Export, ExportError, Exported, MatrixArgs, MatrixMetadata, VectorMetadata};
use std::num::NonZeroUsize;

// These primitives deliberately lack Clone/Eq, and the operator lacks Product/Identity.
struct Ket(u8);
enum Op {
    Flip,
    Parity,
}
#[allow(non_upper_case_globals)]
const Op: Op = Op::Flip;
impl qmb_core::HomogeneousParity for Op {
    fn parity(&self) -> qmb_core::Parity {
        qmb_core::Parity::Even
    }
}
impl qmb_core::ParityOperator for Op {
    fn parity_operator() -> Self {
        Self::Parity
    }
}
struct Ordered;
impl Dimension for Ket {
    fn dim() -> Result<NonZeroUsize, DimensionOverflow> {
        Ok(NonZeroUsize::new(2).unwrap())
    }
}
impl Act<f64> for Ket {
    type Operator = Op;
    type Term = Self;
    fn act(&self, op: &Op) -> LinearCombination<Self, f64> {
        LinearCombination::primitive(Ket(if matches!(op, Op::Parity) {
            self.0
        } else {
            1 - self.0
        }))
    }
}
impl Adjoint<f64> for Op {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, f64> {
        LinearCombination::primitive(match self {
            Self::Flip => Self::Flip,
            Self::Parity => Self::Parity,
        })
    }
}
impl Basis<Ket, f64> for Ordered {
    fn expand(&self, p: &Ket) -> LinearCombination<Ket, f64> {
        LinearCombination::primitive(Ket(p.0))
    }
    fn index(&self, p: &Ket) -> Option<usize> {
        Some(usize::from(p.0))
    }
    fn elements(&self) -> impl Iterator<Item = Ket> {
        [Ket(0), Ket(1)].into_iter()
    }
}

#[derive(PartialEq, Eq)]
struct Index(u8); // No Clone, Hash or Ord.

#[test]
fn exports_need_neither_primitive_cloning_nor_operator_products_or_index_cloning() {
    let op = LocalOperator::primitive(Op);
    let local: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, ()>> =
        (&op).export(MatrixArgs::new(&Ordered)).unwrap();
    assert_eq!(local.data(), &[0., 1., 1., 0.]);
    let a = HilbertSpace::factor(Index(0));
    let b = HilbertSpace::factor(Index(1));
    let space = a.tensor_product(&b).unwrap();
    let op = op.at(&a).unwrap();
    let matrix: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> = (&op)
        .export(MatrixArgs::new(&Ordered).in_space(&space))
        .unwrap();
    assert_eq!(matrix.metadata().shape(), [4, 4]);
    assert_eq!(matrix.data()[2], 1.);
    assert_eq!(matrix.data()[7], 1.);
    assert_eq!(matrix.data().iter().sum::<f64>(), 4.);
    let state = qmb_core::State::tensor_in(
        &space,
        [
            LocalState::<Ket, Op, f64>::primitive(Ket(1))
                .at(&a)
                .unwrap(),
            LocalState::primitive(Ket(0)).at(&b).unwrap(),
        ],
    )
    .unwrap();
    let v: Exported<Vec<_>, VectorMetadata<'_, _, _, Ket>> = (&state).export(&Ordered).unwrap();
    assert_eq!(v.data(), &[0., 0., 1., 0.]);
    #[cfg(feature = "sprs")]
    {
        let sparse: Exported<qmb_export::sprs::CsMat<_>, MatrixMetadata<'_, Ket, _, _>> = (&op)
            .export(MatrixArgs::new(&Ordered).in_space(&space))
            .unwrap();
        assert_eq!(sparse.data().nnz(), 4);
    }
}

struct Broken(u8);
impl Basis<Ket, f64> for Broken {
    fn expand(&self, _: &Ket) -> LinearCombination<Ket, f64> {
        LinearCombination::primitive(Ket(if self.0 == 3 { 3 } else { 0 }))
    }
    fn index(&self, p: &Ket) -> Option<usize> {
        if self.0 == 4 {
            None
        } else {
            Some(usize::from(p.0))
        }
    }
    fn elements(&self) -> impl Iterator<Item = Ket> {
        match self.0 {
            0 => vec![Ket(1), Ket(0)],
            1 => vec![Ket(0)],
            2 => vec![Ket(0), Ket(1), Ket(0)],
            _ => vec![Ket(0), Ket(1)],
        }
        .into_iter()
    }
}

#[test]
fn malformed_bases_and_unbound_or_invalid_spaces_are_errors_even_for_zero() {
    let zero = LocalOperator::<Op, f64>::zero();
    for kind in 0..3 {
        let basis = Broken(kind);
        let result: Result<Exported<Vec<_>, MatrixMetadata<'_, Ket, _, ()>>, _> =
            (&zero).export(MatrixArgs::new(&basis));
        assert!(matches!(
            result,
            Err(ExportError::Basis(BasisError::InvalidEnumeration))
        ));
    }
    let state = LocalState::<Ket, Op, f64>::primitive(Ket(0)).scaled(0.);
    let basis = Broken(3);
    let result: Result<Exported<Vec<_>, VectorMetadata<'_, _, (), Ket>>, _> =
        (&state).export(&basis);
    assert!(matches!(
        result,
        Err(ExportError::Basis(BasisError::IndexOutOfBounds { .. }))
    ));
    let basis = Broken(4);
    let result: Result<Exported<Vec<_>, VectorMetadata<'_, _, (), Ket>>, _> =
        (&state).export(&basis);
    assert!(matches!(
        result,
        Err(ExportError::Basis(BasisError::NonBasisElement))
    ));

    let a = HilbertSpace::factor(0u8);
    let b = HilbertSpace::factor(1u8);
    let op = Operator::<u8, Op, f64>::zero();
    let result: Result<Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>>, _> =
        (&op).export(MatrixArgs::new(&Ordered));
    assert!(matches!(
        result,
        Err(ExportError::Expression(ExpressionError::UnboundOperator))
    ));
    let op = LocalOperator::<Op, f64>::primitive(Op)
        .at(&a)
        .unwrap()
        .scaled(0.);
    let result: Result<Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>>, _> =
        (&op).export(MatrixArgs::new(&Ordered).in_space(&b));
    assert!(matches!(
        result,
        Err(ExportError::Expression(ExpressionError::MissingIndex))
    ));
    let bound = op.bind(&a).unwrap();
    let larger = a.tensor_product(&b).unwrap();
    let matrix: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> = (&bound)
        .export(MatrixArgs::new(&Ordered).in_space(&larger))
        .unwrap();
    assert_eq!(matrix.metadata().shape(), [4, 4]);
    assert!(matrix.metadata().space().shares_node_with(&larger));
    assert_eq!(bound.space().unwrap().indices().count(), 1);
}

#[test]
fn unit_zero_and_identity_preserve_spaces_and_shapes() {
    use num_complex::Complex;
    use qmb_spin_half::{SpinHalfOperator as O, SpinHalfState as S, SzBasis};
    let space = HilbertSpace::<u8>::unit();
    let state = State::<u8, S, O, Complex<f64>>::unit();
    let vector: Exported<Vec<_>, VectorMetadata<'_, _, _, S>> = (&state).export(&SzBasis).unwrap();
    assert_eq!(vector.data(), &[Complex::new(1., 0.)]);
    assert_eq!(vector.metadata().factor_count(), 0);
    assert_eq!(vector.metadata().local_dimension(), None);
    let identity = Operator::<u8, O, Complex<f64>>::identity();
    let matrix: Exported<Vec<_>, MatrixMetadata<'_, S, _, _>> = (&identity)
        .export(MatrixArgs::new(&SzBasis).in_space(&space))
        .unwrap();
    assert_eq!(matrix.data(), vector.data());
    assert_eq!(matrix.metadata().shape(), [1, 1]);
    let zero = State::<u8, S, O, Complex<f64>>::zero(&space);
    let vector: Exported<Vec<_>, VectorMetadata<'_, _, _, S>> = (&zero).export(&SzBasis).unwrap();
    assert_eq!(vector.data(), &[Complex::new(0., 0.)]);
    #[cfg(feature = "sprs")]
    {
        let vector: Exported<qmb_export::sprs::CsVec<_>, VectorMetadata<'_, _, _, S>> =
            (&zero).export(&SzBasis).unwrap();
        assert_eq!(vector.data().dim(), 1);
        assert_eq!(vector.data().nnz(), 0);
        let zero = Operator::<u8, O, Complex<f64>>::zero();
        let matrix: Exported<qmb_export::sprs::CsMat<_>, MatrixMetadata<'_, S, _, _>> = (&zero)
            .export(MatrixArgs::new(&SzBasis).in_space(&space))
            .unwrap();
        assert_eq!(matrix.data().shape(), (1, 1));
        assert_eq!(matrix.data().nnz(), 0);
    }
    let space =
        HilbertSpace::tensor([HilbertSpace::factor(0u8), HilbertSpace::factor(1u8)]).unwrap();
    let identity: Exported<Vec<_>, MatrixMetadata<'_, S, _, _>> = (&identity)
        .export(MatrixArgs::new(&SzBasis).in_space(&space))
        .unwrap();
    for (i, &v) in identity.data().iter().enumerate() {
        assert_eq!(v, Complex::new(if i / 4 == i % 4 { 1. } else { 0. }, 0.));
    }
}

#[test]
fn export_detects_axis_and_dense_storage_overflow_before_allocating() {
    let space = HilbertSpace::tensor((0..usize::BITS).map(HilbertSpace::factor)).unwrap();
    let op = Operator::<u32, Op, f64>::identity();
    let result: Result<Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>>, _> =
        (&op).export(MatrixArgs::new(&Ordered).in_space(&space));
    assert!(matches!(result, Err(ExportError::Dimension(_))));
    let state = State::<u32, Ket, Op, f64>::zero(&space);
    let result: Result<Exported<Vec<_>, VectorMetadata<'_, _, _, Ket>>, _> =
        (&state).export(&Ordered);
    assert!(matches!(result, Err(ExportError::Dimension(_))));

    // The axis fits, but its square does not; no dense allocation is attempted.
    let space = HilbertSpace::tensor((0..usize::BITS / 2 + 1).map(HilbertSpace::factor)).unwrap();
    let result: Result<Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>>, _> =
        (&op).export(MatrixArgs::new(&Ordered).in_space(&space));
    assert!(matches!(result, Err(ExportError::SizeOverflow)));
}

#[test]
fn sparse_coordinates_preserve_small_values_and_non_finite_values() {
    let tiny = f64::MIN_POSITIVE;
    let state = LocalState::<Ket, Op, f64>::primitive(Ket(0)).scaled(tiny)
        + LocalState::primitive(Ket(1)).scaled(f64::NAN);
    let output: Exported<Vec<_>, VectorMetadata<'_, _, (), Ket>> =
        (&state).export(&Ordered).unwrap();
    assert_eq!(output.data()[0], tiny);
    assert!(output.data()[1].is_nan());
    #[cfg(feature = "sprs")]
    {
        let output: Exported<qmb_export::sprs::CsVec<_>, VectorMetadata<'_, _, (), Ket>> =
            (&state).export(&Ordered).unwrap();
        assert_eq!(output.data().nnz(), 2);
        assert_eq!(output.data().data()[0], tiny);
        assert!(output.data().data()[1].is_nan());
    }
}

#[cfg(feature = "sprs")]
#[test]
fn owned_sparse_coordinate_exports_do_not_add_clone_or_numeric_bounds() {
    // This generic function compiles without Clone/Zero or any arithmetic on C.
    fn convert<'b, B, C, H: Clone>(
        coords: Coordinates<'b, B, C, H>,
    ) -> Exported<qmb_export::sprs::CsVec<C>, VectorMetadata<'b, B, H>> {
        coords.export(()).unwrap()
    }
    let basis = Ordered;
    let state = LocalState::<Ket, Op, f64>::primitive(Ket(1));
    let coords = state.coordinates_in(&basis).unwrap();
    let sparse = convert(coords);
    assert_eq!(sparse.data().indices(), &[1]);
}

#[test]
fn deep_spaces_are_walked_without_recursive_export_calls() {
    let mut space = HilbertSpace::factor(1u8);
    for _ in 0..20_000 {
        space = HilbertSpace::tensor([space]).unwrap();
    }
    std::thread::Builder::new()
        .stack_size(128 * 1024)
        .spawn(move || {
            let op = Operator::<u8, Op, f64>::identity();
            let output: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, _>> = (&op)
                .export(MatrixArgs::new(&Ordered).in_space(&space))
                .unwrap();
            assert_eq!(output.data(), &[1., 0., 0., 1.]);
            drop(output);
            drop(space);
        })
        .unwrap()
        .join()
        .unwrap();
}

// Compile-time lifetime checks: exported data borrows the basis, not the input.
#[test]
fn a_result_outlives_its_source_without_borrowing_it() {
    let basis = Ordered;
    let result: Exported<Vec<_>, VectorMetadata<'_, _, (), Ket>> = {
        let state = LocalState::<Ket, Op, f64>::primitive(Ket(1));
        (&state).export(&basis).unwrap()
    };
    assert_eq!(result.data(), &[0., 1.]);
}

#[cfg(feature = "ndarray")]
#[test]
fn owned_ndarray_conversions_add_no_coefficient_or_primitive_bounds() {
    fn vector<'b, B, C, H, P>(
        source: Exported<Vec<C>, VectorMetadata<'b, B, H, P>>,
    ) -> Exported<qmb_export::ndarray::Array1<C>, VectorMetadata<'b, B, H, P>> {
        source.export(()).unwrap()
    }
    fn matrix<'b, S, B, C, H>(
        source: Exported<Vec<C>, MatrixMetadata<'b, S, B, H>>,
    ) -> Exported<qmb_export::ndarray::Array2<C>, MatrixMetadata<'b, S, B, H>> {
        source.export(()).unwrap()
    }
    let basis = Ordered;
    let state = LocalState::<Ket, Op, f64>::primitive(Ket(1));
    let source: Exported<Vec<_>, VectorMetadata<'_, _, (), Ket>> = (&state).export(&basis).unwrap();
    let pointer = source.data().as_ptr();
    assert_eq!(vector(source).data().as_ptr(), pointer);
    let op = LocalOperator::<Op, f64>::primitive(Op);
    let source: Exported<Vec<_>, MatrixMetadata<'_, Ket, _, ()>> =
        (&op).export(MatrixArgs::new(&basis)).unwrap();
    let pointer = source.data().as_ptr();
    assert_eq!(matrix(source).data().as_ptr(), pointer);
}

#[test]
fn export_errors_retain_the_original_error_as_their_source() {
    use std::error::Error;
    let e = ExportError::from(BasisError::InvalidEnumeration);
    assert!(e.source().unwrap().is::<BasisError>());
    assert_eq!(e.to_string(), BasisError::InvalidEnumeration.to_string());
    let e = ExportError::from(ExpressionError::MissingIndex);
    assert!(e.source().unwrap().is::<ExpressionError>());
    let e = ExportError::from(qmb_core::SpaceError::DuplicateIndex);
    assert!(e.source().unwrap().is::<qmb_core::SpaceError>());
    let e = ExportError::from(DimensionOverflow);
    assert!(e.source().unwrap().is::<DimensionOverflow>());
    let reserve = Vec::<u8>::new().try_reserve(usize::MAX).unwrap_err();
    let e = ExportError::from(reserve);
    assert!(
        e.source()
            .unwrap()
            .is::<std::collections::TryReserveError>()
    );
    #[cfg(feature = "ndarray")]
    {
        let shape = qmb_export::ndarray::Array2::from_shape_vec((2, 2), vec![1.]).unwrap_err();
        let e = ExportError::from(shape);
        assert!(e.source().unwrap().is::<qmb_export::ndarray::ShapeError>());
    }
    #[cfg(feature = "sprs")]
    {
        let (_, _, structure) =
            qmb_export::sprs::CsVec::try_new(2, vec![2usize], vec![1.]).unwrap_err();
        let e = ExportError::from(structure);
        assert!(
            e.source()
                .unwrap()
                .is::<qmb_export::sprs::errors::StructureError>()
        );
    }
}
