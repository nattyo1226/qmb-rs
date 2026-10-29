mod support;
use num_complex::Complex;
use qmb_core::{
    Act, Adjoint, Basis, Dimension, HomogeneousParity, Identity, InnerProduct, LinearCombination,
    Parity, ParityOperator, Product, Trace,
};
use qmb_fermion::{
    Fermion, FermionOperator as O, FermionState as S, LadderBasis, MajoranaBasis, OccupationBasis,
};
use support::*;

fn primitive_rules<T: Real>() {
    for (i, s) in STATES.into_iter().enumerate() {
        for (j, t) in STATES.into_iter().enumerate() {
            close::<T>(
                to64(<S as InnerProduct<Complex<T>>>::inner_product(&s, &t)),
                (if i == j { 1. } else { 0. }).into(),
            );
        }
    }
    for op in OPS {
        let m = matrix(op);
        let adj = <O as Adjoint<Complex<T>>>::adjoint(&op);
        values::<T>(&combination(adj), &adjoint(&m));
        close::<T>(to64(<O as Trace<Complex<T>>>::trace(&op)), trace(&m));
        for (j, s) in STATES.into_iter().enumerate() {
            let action: LinearCombination<S, Complex<T>> = s.act(&op);
            let mut v = [Complex::new(0., 0.); 2];
            for term in action {
                assert_eq!(
                    term.operand.parity() == Parity::Odd,
                    (s.parity() == Parity::Odd) != (op.parity() == Parity::Odd)
                );
                v[usize::from(term.operand == S::Occupied)] += to64(term.coefficient);
            }
            values::<T>(&v, &[m[j], m[2 + j]]);
        }
        for rhs in OPS {
            let product: LinearCombination<O, Complex<T>> = op.product(&rhs);
            for term in product.terms() {
                assert_eq!(
                    term.operand.parity() == Parity::Odd,
                    (op.parity() == Parity::Odd) != (rhs.parity() == Parity::Odd)
                );
            }
            values::<T>(&combination(product), &support::product(&m, &matrix(rhs)));
            close::<T>(
                trace(&support::product(&m, &matrix(rhs))),
                trace(&support::product(&matrix(rhs), &m)),
            );
        }
    }
}

fn laws<T: Real>() {
    let op = Fermion::<Complex<T>>::operator;
    for a in OPS {
        for b in OPS {
            let ab = op(a).clone() * &op(b);
            close::<T>(to64(ab.trace()), to64((op(b).clone() * &op(a)).trace()));
            close::<T>(
                to64(ab.trace()),
                trace(&support::product(&matrix(a), &matrix(b))),
            );
            values::<T>(
                &local_value(&ab.adjoint()),
                &local_value(&(op(b).adjoint().clone() * &op(a).adjoint())),
            );
            values::<T>(&local_value(&op(a).adjoint().adjoint()), &matrix(a));
            for s in STATES {
                let ket = Fermion::<Complex<T>>::state(s);
                let combined = expression::<T>(a.product(&b));
                let l = ket.act(&combined);
                let r = ket.act(&op(b)).act(&op(a));
                for target in STATES {
                    let bra = Fermion::<Complex<T>>::state(target);
                    close::<T>(to64(bra.inner_product(&l)), to64(bra.inner_product(&r)));
                    close::<T>(
                        to64(bra.inner_product(&ket.act(&op(a)))),
                        to64(bra.act(&op(a).adjoint()).inner_product(&ket)),
                    );
                }
            }
            for d in OPS {
                values::<T>(
                    &local_value(&(ab.clone() * &op(d))),
                    &local_value(&(op(a).clone() * &(op(b).clone() * &op(d)))),
                );
            }
        }
        values::<T>(&local_value(&(op(O::I).clone() * &op(a))), &matrix(a));
        values::<T>(&local_value(&(op(a).clone() * &op(O::I))), &matrix(a));
    }
    let [a, d, n, p, x, y, z] = [
        O::Annihilate,
        O::Create,
        O::Number,
        O::Parity,
        O::MajoranaX,
        O::MajoranaY,
        O::MajoranaZ,
    ]
    .map(op);
    let zero = vec![0.0.into(); 4];
    values::<T>(&local_value(&(a.clone() * &a)), &zero);
    values::<T>(&local_value(&(d.clone() * &d)), &zero);
    values::<T>(&local_value(&a.anticommutator(&d)), &matrix(O::I));
    values::<T>(&local_value(&(d.clone() * &a)), &matrix(O::Number));
    values::<T>(&local_value(&(n.clone() * &n)), &matrix(O::Number));
    values::<T>(&local_value(&(p.clone() * &p)), &matrix(O::I));
    values::<T>(&local_value(&p.anticommutator(&a)), &zero);
    values::<T>(&local_value(&p.anticommutator(&d)), &zero);
    for (l, r, out) in [(&x, &y, &z), (&y, &z, &x), (&z, &x, &y)] {
        let comm = (l.clone() * r) + (r.clone() * l).scaled(c::<T>(-1., 0.));
        values::<T>(
            &local_value(&comm),
            &local_value(&out.scaled(c::<T>(0., 1.))),
        );
    }
    values::<T>(
        &local_value(&((x.clone() * &x) + (y.clone() * &y) + (z.clone() * &z))),
        &scale(&matrix(O::I), 0.75.into()),
    );
    values::<T>(
        &local_value(&(x.clone() * &y).scaled(c::<T>(0., -2.))),
        &matrix(O::MajoranaZ),
    );
}

fn bases<T: Real>() {
    fn check<T: Real, B: Basis<O, Complex<T>>>(basis: B, expected: &[O]) {
        let elements: Vec<_> = basis.elements().collect();
        assert_eq!(elements, expected);
        for primitive in OPS {
            assert_eq!(
                basis.index(&primitive),
                elements.iter().position(|p| *p == primitive)
            );
            let expansion = basis.expand(&primitive);
            for term in expansion.terms() {
                assert!(elements.contains(&term.operand));
            }
            values::<T>(&combination(expansion), &matrix(primitive));
            let op = Fermion::<Complex<T>>::operator(primitive);
            let expanded = op.expand(&basis).unwrap();
            values::<T>(&local_value(&expanded), &matrix(primitive));
            let round = expanded
                .expand(&LadderBasis)
                .unwrap()
                .expand(&MajoranaBasis)
                .unwrap()
                .expand(&basis)
                .unwrap();
            values::<T>(&local_value(&round), &matrix(primitive));
            close::<T>(to64(expanded.trace()), trace(&matrix(primitive)));
        }
    }
    check::<T, _>(LadderBasis, &[O::I, O::Annihilate, O::Create, O::Number]);
    check::<T, _>(
        MajoranaBasis,
        &[O::I, O::MajoranaX, O::MajoranaY, O::MajoranaZ],
    );
    let states: Vec<_> =
        <OccupationBasis as Basis<S, Complex<T>>>::elements(&OccupationBasis).collect();
    assert_eq!(states, STATES);
    for (j, state) in states.iter().enumerate() {
        assert_eq!(
            <OccupationBasis as Basis<S, Complex<T>>>::index(&OccupationBasis, state),
            Some(j)
        );
        let expanded = <OccupationBasis as Basis<S, Complex<T>>>::expand(&OccupationBasis, state);
        assert_eq!(expanded.into_iter().next().unwrap().operand, *state);
    }
    let majoranas =
        [O::I, O::MajoranaX, O::MajoranaY, O::MajoranaZ].map(Fermion::<Complex<T>>::operator);
    for (i, a) in majoranas.iter().enumerate() {
        for (j, b) in majoranas.iter().enumerate() {
            close::<T>(
                to64(a.hs_inner_product(b)),
                (if i != j {
                    0.
                } else if i == 0 {
                    2.
                } else {
                    0.5
                })
                .into(),
            );
        }
    }
}

#[test]
fn full_rules_laws_and_bases_f32() {
    primitive_rules::<f32>();
    laws::<f32>();
    bases::<f32>();
}
#[test]
fn full_rules_laws_and_bases_f64() {
    primitive_rules::<f64>();
    laws::<f64>();
    bases::<f64>();
}
#[test]
fn public_metadata() {
    fn value<T: Clone + Copy + core::fmt::Debug + Eq + core::hash::Hash>() {}
    value::<O>();
    value::<S>();
    value::<LadderBasis>();
    value::<MajoranaBasis>();
    value::<OccupationBasis>();
    assert_eq!(S::dim().unwrap().get(), 2);
    assert_eq!(O::dim().unwrap().get(), 4);
    assert_eq!(O::identity(), O::I);
    assert_eq!(O::parity_operator(), O::Parity);
    assert_eq!(S::Empty.parity(), Parity::Even);
    assert_eq!(S::Occupied.parity(), Parity::Odd);
    for op in OPS {
        assert_eq!(
            op.parity() == Parity::Odd,
            matches!(op, O::Annihilate | O::Create | O::MajoranaX | O::MajoranaY)
        );
    }
}
