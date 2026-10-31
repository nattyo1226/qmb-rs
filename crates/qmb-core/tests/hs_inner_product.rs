//! Normalized HS contractions preserve ambient-space and grading semantics.
#[path = "support/graded.rs"]
mod support;

use num_complex::Complex;
use qmb_core::{
    Adjoint, DimensionOverflow, HilbertSpace, HomogeneousParity, Identity, LinearCombination,
    LocalOperator, Operator, Parity, ParityOperator, Product, TensorKind, Trace,
};
use support::{Op, Real, c, close};

// No Dimension, state type, or basis implementation: HS gets d from Tr(I).
#[derive(Clone, Copy)]
struct Primitive(Op);
impl Identity for Primitive {
    fn identity() -> Self {
        Self(Op::I)
    }
}
impl HomogeneousParity for Primitive {
    fn parity(&self) -> Parity {
        self.0.parity()
    }
}
impl ParityOperator for Primitive {
    fn parity_operator() -> Self {
        Self(Op::P)
    }
}
impl<T: Real> Product<Complex<T>> for Primitive {
    type Term = Self;
    fn product(&self, rhs: &Self) -> LinearCombination<Self, Complex<T>> {
        self.0.product(&rhs.0).map(Self)
    }
}
impl<T: Real> Adjoint<Complex<T>> for Primitive {
    type Term = Self;
    fn adjoint(&self) -> LinearCombination<Self, Complex<T>> {
        self.0.adjoint().map(Self)
    }
}
impl<T: Real> Trace<Complex<T>> for Primitive {
    fn trace(&self) -> Complex<T> {
        self.0.trace()
    }
}

fn contractions<T: Real>() {
    type Local<T> = LocalOperator<Primitive, Complex<T>>;
    for a in support::OPS {
        for b in support::OPS {
            let lhs = Local::<T>::primitive(Primitive(a)).scaled(c(0.3, -0.7));
            let rhs = Local::<T>::primitive(Primitive(b)).scaled(c(-0.2, 0.4));
            let reference = support::trace(&support::product(
                &support::adjoint(&support::local::<T>(a)),
                &support::local(b),
            )) * c::<T>(0.3, 0.7)
                * c::<T>(-0.2, 0.4)
                / c(2.0, 0.0);
            close(lhs.hs_inner_product(&rhs), reference);
            for kind in [TensorKind::Ordinary, TensorKind::Graded] {
                let factors = [HilbertSpace::factor(0u8), HilbertSpace::factor(1)];
                let space = HilbertSpace::tensor_all_with(kind, factors.iter().cloned()).unwrap();
                let left = lhs.at(&factors[1]).unwrap().bind(&space).unwrap();
                let right = rhs.at(&factors[1]).unwrap().bind(&space).unwrap();
                close(left.hs_inner_product(&right).unwrap(), reference);
                // Different placements exercise nontrivial parity strings.
                let other = rhs.at(&factors[0]).unwrap();
                let dense = support::trace(&support::product(
                    &support::adjoint(&support::embedding::<T>(a, 1, 2, kind)),
                    &support::embedding(b, 0, 2, kind),
                )) * c::<T>(0.3, 0.7)
                    * c::<T>(-0.2, 0.4)
                    / c(4.0, 0.0);
                close(left.hs_inner_product_in(&other, &space).unwrap(), dense);
            }
        }
    }
}

#[test]
fn local_and_graded_contractions_match_independent_matrices() {
    contractions::<f32>();
    contractions::<f64>();
}

#[test]
fn explicit_space_overrides_bindings_and_keeps_spectators_normalized() {
    type Local = LocalOperator<Primitive, Complex<f64>>;
    let factors = [
        HilbertSpace::factor(0u8),
        HilbertSpace::factor(1),
        HilbertSpace::factor(2),
    ];
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        let old = HilbertSpace::tensor_all_with(kind, factors[..2].iter().cloned()).unwrap();
        let full = HilbertSpace::tensor_all_with(kind, factors.iter().cloned()).unwrap();
        let op = Local::primitive(Primitive(support::A))
            .at(&factors[1])
            .unwrap()
            .bind(&old)
            .unwrap();
        close(op.hs_inner_product(&op).unwrap(), c(0.5, 0.0));
        close(op.hs_inner_product_in(&op, &full).unwrap(), c(0.5, 0.0));
        let number = Local::primitive(Primitive(Op::U(1, 1)))
            .at(&factors[1])
            .unwrap();
        close(number.trace_in(&old).unwrap(), c(2.0, 0.0));
        close(number.trace_in(&full).unwrap(), c(4.0, 0.0));
    }
}

#[test]
fn scalar_unit_spaces_and_dimensions_beyond_usize_and_float_ranges() {
    type Op = Operator<usize, Primitive, Complex<f64>>;
    for kind in [TensorKind::Ordinary, TensorKind::Graded] {
        let empty = HilbertSpace::<usize>::tensor_all_with(kind, []).unwrap();
        let scalar = Op::scalar(c(3.0, 4.0)).bind(&empty).unwrap();
        close(scalar.hs_inner_product(&scalar).unwrap(), c(25.0, 0.0));
        close(scalar.trace().unwrap(), c(3.0, 4.0));
        let large =
            HilbertSpace::tensor_all_with(kind, (0..2048).map(HilbertSpace::factor)).unwrap();
        assert_eq!(large.dim::<support::Ket>(), Err(DimensionOverflow));
        let identity = Op::identity().bind(&large).unwrap();
        close(identity.hs_inner_product(&identity).unwrap(), c(1.0, 0.0));
        let zero = Op::zero().bind(&large).unwrap();
        close(identity.hs_inner_product(&zero).unwrap(), c(0.0, 0.0));
    }
}
