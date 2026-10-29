mod support;

use num_complex::Complex;
use qmb_core::{
    Act, Adjoint, Basis, Dimension, HomogeneousParity, Identity, InnerProduct,
    LinearCombination as Lc, Parity, ParityOperator, Product, Trace,
};
use qmb_spin_one::{
    CartesianBasis, LadderBasis, SpinOne, SpinOneOperator as O, SpinOneState as S, SzBasis,
};
use support::*;

fn lc_product<T: Real>(a: &Lc<O, Complex<T>>, b: &Lc<O, Complex<T>>) -> Lc<O, Complex<T>> {
    Lc::from_terms(a.terms().flat_map(|a| {
        b.terms().flat_map(move |b| {
            <O as Product<Complex<T>>>::product(&a.operand, &b.operand)
                .scaled(a.coefficient * b.coefficient)
        })
    }))
}
fn lc_adjoint<T: Real>(a: &Lc<O, Complex<T>>) -> Lc<O, Complex<T>> {
    Lc::from_terms(
        a.terms().flat_map(|a| {
            <O as Adjoint<Complex<T>>>::adjoint(&a.operand).scaled(a.coefficient.conj())
        }),
    )
}
fn lc_action<T: Real>(op: &Lc<O, Complex<T>>, state: &Lc<S, Complex<T>>) -> Lc<S, Complex<T>> {
    Lc::from_terms(op.terms().flat_map(|a| {
        state.terms().flat_map(move |s| {
            <S as Act<Complex<T>>>::act(&s.operand, &a.operand)
                .scaled(a.coefficient * s.coefficient)
        })
    }))
}
fn primitive_rules<T: Real>() {
    for s in STATES {
        for t in STATES {
            close::<T>(
                to64(<S as InnerProduct<Complex<T>>>::inner_product(&s, &t)),
                f64::from(s == t).into(),
            );
        }
    }
    for a in OPS {
        let m = matrix(a);
        let adj = <O as Adjoint<Complex<T>>>::adjoint(&a);
        assert_eq!(adj.len(), 1);
        values::<T>(&combination(adj.clone()), &adjoint(&m));
        values::<T>(&combination(lc_adjoint(&adj)), &m);
        close::<T>(to64(<O as Trace<Complex<T>>>::trace(&a)), trace(&m));
        for s in STATES {
            let acted = <S as Act<Complex<T>>>::act(&s, &a);
            let expected = apply(&m, &ket(s));
            if expected.iter().all(|v| v.norm() < 1e-14) {
                assert!(acted.is_empty(), "{a:?} {s:?}");
            }
            assert!(acted.terms().all(|t| t.coefficient != c(0., 0.)));
            values::<T>(&state_combination(acted), &expected);
        }
        for b in OPS {
            let actual = <O as Product<Complex<T>>>::product(&a, &b);
            let expected = product(&m, &matrix(b));
            assert!(actual.len() <= 9);
            assert!(actual.terms().all(|t| t.coefficient != c(0., 0.)));
            if a == O::I {
                assert_eq!(actual, Lc::primitive(b));
            }
            if b == O::I {
                assert_eq!(actual, Lc::primitive(a));
            }
            let basis: &[O] = if COMMON.contains(&a) && COMMON.contains(&b) {
                &COMMON
            } else if LADDER.contains(&a) && LADDER.contains(&b) {
                &LADDER
            } else {
                &CART
            };
            let indices: Vec<_> = actual
                .terms()
                .map(|t| {
                    basis
                        .iter()
                        .position(|p| *p == t.operand)
                        .expect("product left its prescribed basis")
                })
                .collect();
            assert!(indices.windows(2).all(|i| i[0] < i[1]));
            if expected.iter().all(|v| v.norm() < 1e-14) {
                assert!(actual.is_empty(), "{a:?} {b:?}");
            }
            values::<T>(&combination(actual), &expected);
        }
    }
}

fn laws<T: Real>() {
    for a in OPS {
        let al = Lc::<_, Complex<T>>::primitive(a);
        let am = matrix(a);
        for s in STATES {
            for t in STATES {
                let psi = SpinOne::<Complex<T>>::state(s);
                let phi = SpinOne::<Complex<T>>::state(t);
                let op = SpinOne::<Complex<T>>::operator(a);
                close::<T>(
                    to64(phi.inner_product(&psi.act(&op))),
                    to64(phi.act(&op.adjoint()).inner_product(&psi)),
                );
            }
        }
        for b in OPS {
            let bl = Lc::primitive(b);
            let ab = lc_product(&al, &bl);
            let bm = matrix(b);
            close::<T>(
                trace(&combination(ab.clone())),
                trace(&combination(lc_product(&bl, &al))),
            );
            values::<T>(
                &combination(lc_adjoint(&ab)),
                &combination(lc_product(&lc_adjoint(&bl), &lc_adjoint(&al))),
            );
            for s in STATES {
                let ket = Lc::primitive(s);
                let lhs = lc_action(&ab, &ket);
                let rhs = lc_action(&al, &lc_action(&bl, &ket));
                let expected = apply(&product(&am, &bm), &support::ket(s));
                values::<T>(&state_combination(lhs), &expected);
                values::<T>(&state_combination(rhs), &expected);
            }
            for d in OPS {
                let dl = Lc::primitive(d);
                let lhs = lc_product(&ab, &dl);
                let rhs = lc_product(&al, &lc_product(&bl, &dl));
                let expected = product(&product(&am, &bm), &matrix(d));
                values::<T>(&combination(lhs), &expected);
                values::<T>(&combination(rhs), &expected);
            }
        }
    }
}

fn spin_relations<T: Real>() {
    let op = SpinOne::<Complex<T>>::operator;
    let [x, y, z, p, m, q0] = [O::Sx, O::Sy, O::Sz, O::SPlus, O::SMinus, O::Q0].map(op);
    for (a, b, result) in [(&x, &y, &z), (&y, &z, &x), (&z, &x, &y)] {
        let comm = (a.clone() * b) - (b.clone() * a);
        values::<T>(
            &local_value(&comm),
            &local_value(&result.scaled(c::<T>(0., 1.))),
        );
    }
    for (ladder, sign) in [(&p, 1.), (&m, -1.)] {
        values::<T>(
            &local_value(&((z.clone() * ladder) - (ladder.clone() * &z))),
            &local_value(&ladder.scaled(c::<T>(sign, 0.))),
        );
        values::<T>(
            &local_value(&((ladder.clone() * ladder).clone() * ladder)),
            &[0.0.into(); 9],
        );
    }
    values::<T>(
        &local_value(&((p.clone() * &m) - (m.clone() * &p))),
        &scale(&matrix(O::Sz), 2.0.into()),
    );
    values::<T>(
        &local_value(&((x.clone() * &x) + (y.clone() * &y) + (z.clone() * &z))),
        &scale(&matrix(O::I), 2.0.into()),
    );
    values::<T>(
        &local_value(&(z.clone() * &z)),
        &local_value(&(op(O::I).scaled(c::<T>(2. / 3., 0.)) + q0.scaled(c::<T>(1. / 3., 0.)))),
    );
    for a in [&x, &y, &z] {
        close::<T>(to64((a.clone() * a).trace()), 2.0.into());
    }
    assert_eq!(
        <O as Product<Complex<T>>>::product(&O::SPlus, &O::SPlus),
        Lc::primitive(O::Q2Plus)
    );
    assert_eq!(
        <O as Product<Complex<T>>>::product(&O::SMinus, &O::SMinus),
        Lc::primitive(O::Q2Minus)
    );
    for (a, sign) in [
        (O::Q1Plus, 1.),
        (O::Q1Minus, -1.),
        (O::Q2Plus, 2.),
        (O::Q2Minus, -2.),
    ] {
        values::<T>(
            &local_value(&((z.clone() * &op(a)) - (op(a).clone() * &z))),
            &scale(&matrix(a), sign.into()),
        );
    }
}

fn bases<T: Real>() {
    fn check<T: Real, B: Basis<O, Complex<T>>>(basis: B, expected: [O; 9], norms: [f64; 9]) {
        let elements: Vec<_> = basis.elements().collect();
        assert_eq!(elements, expected);
        for p in OPS {
            assert_eq!(basis.index(&p), expected.iter().position(|v| *v == p));
            let expansion = basis.expand(&p);
            assert!(expansion.terms().all(|t| expected.contains(&t.operand)));
            values::<T>(&combination(expansion), &matrix(p));
            let expr = SpinOne::<Complex<T>>::operator(p);
            let expanded = expr.expand(&basis).unwrap();
            let round = expanded
                .expand(&LadderBasis)
                .unwrap()
                .expand(&CartesianBasis)
                .unwrap()
                .expand(&basis)
                .unwrap();
            values::<T>(&local_value(&round), &matrix(p));
            close::<T>(to64(expanded.trace()), trace(&matrix(p)));
            let coords = expr.coordinates(&basis).unwrap();
            assert_eq!(coords.dim().unwrap(), 9);
            let restored = coords
                .flat_entries()
                .fold(LocalOp::<T>::zero(), |sum, entry| {
                    let (i, c) = entry.unwrap();
                    sum + SpinOne::<Complex<T>>::operator(expected[i]).scaled(*c)
                });
            values::<T>(&local_value(&restored), &matrix(p));
        }
        for (i, a) in expected.into_iter().enumerate() {
            for (j, b) in expected.into_iter().enumerate() {
                close::<T>(
                    to64(
                        SpinOne::<Complex<T>>::operator(a)
                            .hs_inner_product(&SpinOne::<Complex<T>>::operator(b)),
                    ),
                    if i == j { norms[i] } else { 0. }.into(),
                );
            }
        }
    }
    check::<T, _>(CartesianBasis, CART, [3., 2., 2., 2., 2., 6., 2., 2., 2.]);
    check::<T, _>(LadderBasis, LADDER, [3., 4., 4., 2., 6., 4., 4., 4., 4.]);
    let states: Vec<_> = <SzBasis as Basis<S, Complex<T>>>::elements(&SzBasis).collect();
    assert_eq!(states, STATES);
    for (i, s) in STATES.into_iter().enumerate() {
        assert_eq!(
            <SzBasis as Basis<S, Complex<T>>>::index(&SzBasis, &s),
            Some(i)
        );
        assert_eq!(
            <SzBasis as Basis<S, Complex<T>>>::expand(&SzBasis, &s),
            Lc::primitive(s)
        );
    }
}

#[test]
fn rules_f32() {
    primitive_rules::<f32>();
}
#[test]
fn rules_f64() {
    primitive_rules::<f64>();
}
#[test]
fn algebra_laws_f32() {
    laws::<f32>();
}
#[test]
fn algebra_laws_f64() {
    laws::<f64>();
}
#[test]
fn spin_relations_f32() {
    spin_relations::<f32>();
}
#[test]
fn spin_relations_f64() {
    spin_relations::<f64>();
}
#[test]
fn bases_f32() {
    bases::<f32>();
}
#[test]
fn bases_f64() {
    bases::<f64>();
}
#[test]
fn public_metadata() {
    fn value<T: Clone + Copy + core::fmt::Debug + Eq + core::hash::Hash>() {}
    fn marker<T: Default>() {}
    value::<O>();
    value::<S>();
    value::<CartesianBasis>();
    value::<LadderBasis>();
    value::<SzBasis>();
    marker::<CartesianBasis>();
    marker::<LadderBasis>();
    marker::<SzBasis>();
    assert_eq!(S::dim().unwrap().get(), 3);
    assert_eq!(O::dim().unwrap().get(), 9);
    assert_eq!(O::identity(), O::I);
    assert_eq!(O::parity_operator(), O::I);
    for s in STATES {
        assert_eq!(s.parity(), Parity::Even);
    }
    for op in OPS {
        assert_eq!(op.parity(), Parity::Even);
    }
}
