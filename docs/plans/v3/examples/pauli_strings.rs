//! Executable design prototype; Relation and SupportCombinations are proposed APIs.
//! One site is one spin. All supports up to CAP, including identity, are visited.

use num_complex::Complex64;
use qmb_core::{HilbertSpace, OperatorOf};
use qmb_spin_half::{SpinHalf, SpinHalfOperator as Op};

type Error = Box<dyn std::error::Error>;
type PauliString<I, C> = OperatorOf<SpinHalf<C>, I>;

trait Relation {
    type Group<'a>: AsRef<[usize]>
    where
        Self: 'a;

    fn num_sites(&self) -> usize;
    fn groups(&self) -> impl Iterator<Item = Self::Group<'_>> + '_;
}

#[derive(Clone, Copy)]
struct SiteGroup<const CAP: usize> {
    sites: [usize; CAP],
    len: usize,
}

impl<const CAP: usize> AsRef<[usize]> for SiteGroup<CAP> {
    fn as_ref(&self) -> &[usize] {
        &self.sites[..self.len]
    }
}

struct SupportCombinations<const CAP: usize> {
    num_sites: usize,
}

impl<const CAP: usize> SupportCombinations<CAP> {
    fn new(num_sites: usize) -> Self {
        Self { num_sites }
    }
}

impl<const CAP: usize> Relation for SupportCombinations<CAP> {
    type Group<'a> = SiteGroup<CAP>;

    fn num_sites(&self) -> usize {
        self.num_sites
    }

    fn groups(&self) -> impl Iterator<Item = Self::Group<'_>> + '_ {
        let n = self.num_sites;
        let max_len = CAP.min(n);
        let mut current = Some(SiteGroup {
            sites: [0; CAP],
            len: 0,
        });
        std::iter::from_fn(move || {
            let result = current.take()?;
            let mut next = result;
            for pos in (0..result.len).rev() {
                if next.sites[pos] < n - result.len + pos {
                    next.sites[pos] += 1;
                    for j in pos + 1..result.len {
                        next.sites[j] = next.sites[j - 1] + 1;
                    }
                    current = Some(next);
                    return Some(result);
                }
            }
            if result.len < max_len {
                next.len += 1;
                for j in 0..next.len {
                    next.sites[j] = j;
                }
                current = Some(next);
            }
            Some(result)
        })
    }
}

// All words over Sx/Sy/Sz of the specified length, without computing 3^len.
// Length zero produces one empty word, which becomes the identity operator.
fn axis_assignments(len: usize) -> impl Iterator<Item = Vec<Op>> {
    std::iter::successors(Some(vec![Op::Sx; len]), |axes| {
        let mut next = axes.clone();
        let pos = next.iter().rposition(|axis| *axis != Op::Sz)?;
        next[pos] = if next[pos] == Op::Sx { Op::Sy } else { Op::Sz };
        next[pos + 1..].fill(Op::Sx);
        Some(next)
    })
}

// The caller needs no numeric evaluation bounds, and Index need not be Clone.
// Supports must be distinct sets for the output strings to be globally unique.
fn for_each_pauli_string<C, I, R>(
    relation: &R,
    sites: &[HilbertSpace<I>],
    mut emit: impl FnMut(PauliString<I, C>),
) -> Result<(), Error>
where
    C: From<f64>,
    I: Eq,
    R: Relation,
{
    if relation.num_sites() != sites.len() {
        return Err("site count mismatch".into());
    }
    if sites.iter().any(|site| site.factor_index().is_none()) {
        return Err("each site must be one spin (one leaf)".into());
    }
    // Shares handles; does not clone Index. Also validates index uniqueness.
    let space = HilbertSpace::tensor_all(sites.iter().cloned())?;

    for group in relation.groups() {
        let support = group.as_ref();
        if support.iter().any(|&i| i >= sites.len()) {
            return Err("site number out of range".into());
        }
        if support
            .iter()
            .enumerate()
            .any(|(pos, i)| support[..pos].contains(i))
        {
            return Err("a support must not repeat a site".into());
        }

        for axes in axis_assignments(support.len()) {
            let factors = support.iter().zip(axes).map(|(&i, axis)| {
                // X = 2 Sx, Y = 2 Sy, Z = 2 Sz.
                SpinHalf::<C>::operator(axis).scaled(2.0).at(&sites[i])
            });
            let word = SpinHalf::<C>::try_operator_product(factors)?.bind(&space)?;
            emit(word);
        }
    }
    Ok(())
}

fn main() -> Result<(), Error> {
    let sites: Vec<_> = (0..3).map(HilbertSpace::factor).collect();
    let relation = SupportCombinations::<2>::new(sites.len());
    let mut strings = Vec::new();

    for_each_pauli_string::<Complex64, _, _>(&relation, &sites, |word| {
        strings.push(word);
    })?;

    // I + 3 sites * 3 axes + 3 pairs * 9 axis assignments.
    assert_eq!(strings.len(), 37);
    println!(
        "3 sites, support size <= 2: {} Pauli strings",
        strings.len()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use qmb_export::{Export, Exported, MatrixArgs};
    use qmb_spin_half::{SpinHalfState, SzBasis};

    fn matrix(word: &PauliString<usize, Complex64>) -> Vec<Complex64> {
        let result: Exported<Vec<Complex64>, _> = word
            .export(MatrixArgs::<SpinHalfState, _>::new(&SzBasis))
            .unwrap();
        result.into_data()
    }

    // Independent Pauli matrices. Labels: 0=I, 1=X, 2=Y, 3=Z.
    fn reference(labels: &[usize]) -> Vec<Complex64> {
        let dim = 1usize << labels.len();
        (0..dim * dim)
            .map(|flat| {
                let (row, col) = (flat / dim, flat % dim);
                labels
                    .iter()
                    .enumerate()
                    .fold(Complex64::new(1.0, 0.0), |entry, (pos, &label)| {
                        let shift = labels.len() - 1 - pos;
                        let (r, c) = ((row >> shift) & 1, (col >> shift) & 1);
                        let local = match (label, r, c) {
                            (0, 0, 0) | (0, 1, 1) | (1, 0, 1) | (1, 1, 0) | (3, 0, 0) => {
                                Complex64::new(1.0, 0.0)
                            }
                            (2, 0, 1) => Complex64::new(0.0, -1.0),
                            (2, 1, 0) => Complex64::new(0.0, 1.0),
                            (3, 1, 1) => Complex64::new(-1.0, 0.0),
                            _ => Complex64::new(0.0, 0.0),
                        };
                        entry * local
                    })
            })
            .collect()
    }

    #[test]
    fn all_37_strings_match_independent_pauli_matrices_once() {
        let sites: Vec<_> = (0..3).map(HilbertSpace::factor).collect();
        let mut expected: Vec<_> = (0usize..64)
            .filter_map(|code| {
                let labels = [(code / 16) % 4, (code / 4) % 4, code % 4];
                (labels.iter().filter(|&&x| x != 0).count() <= 2).then(|| reference(&labels))
            })
            .collect();
        assert_eq!(expected.len(), 37);
        let mut count = 0;
        for_each_pauli_string::<Complex64, _, _>(
            &SupportCombinations::<2>::new(3),
            &sites,
            |word| {
                let actual = matrix(&word);
                let matching = expected
                    .iter()
                    .position(|expected| *expected == actual)
                    .expect("missing, duplicate or incorrectly normalized Pauli string");
                expected.remove(matching);
                assert_eq!(word.space().unwrap().dim::<SpinHalfState>().unwrap(), 8);
                count += 1;
            },
        )
        .unwrap();
        assert_eq!(count, 37);
        assert!(expected.is_empty());
    }

    #[test]
    fn zero_capacity_empty_system_and_cap_larger_than_system() {
        for (n, expected_count) in [(0, 1), (2, 16)] {
            let sites: Vec<_> = (0..n).map(HilbertSpace::factor).collect();
            let mut count = 0;
            let relation = SupportCombinations::<4>::new(n);
            for_each_pauli_string::<Complex64, _, _>(&relation, &sites, |word| {
                assert_eq!(
                    word.space().unwrap().dim::<SpinHalfState>().unwrap(),
                    1 << n
                );
                if n == 0 {
                    assert_eq!(matrix(&word), vec![Complex64::new(1.0, 0.0)]);
                }
                count += 1;
            })
            .unwrap();
            assert_eq!(count, expected_count);
        }
        let sites: Vec<_> = (0..2).map(HilbertSpace::factor).collect();
        let mut matrices = Vec::new();
        for_each_pauli_string::<Complex64, _, _>(
            &SupportCombinations::<0>::new(2),
            &sites,
            |word| matrices.push(matrix(&word)),
        )
        .unwrap();
        assert_eq!(matrices, vec![reference(&[0, 0])]);
    }

    #[test]
    fn construction_needs_neither_clone_on_index_or_coefficient_nor_evaluation() {
        #[derive(PartialEq, Eq)]
        struct Index(usize);
        struct Coefficient;
        impl From<f64> for Coefficient {
            fn from(value: f64) -> Self {
                assert_eq!(value, 2.0);
                Self
            }
        }
        let sites = [
            HilbertSpace::factor(Index(10)),
            HilbertSpace::factor(Index(4)),
        ];
        let mut count = 0;
        for_each_pauli_string::<Coefficient, _, _>(
            &SupportCombinations::<2>::new(2),
            &sites,
            |word| {
                assert!(word.space().is_some());
                count += 1;
            },
        )
        .unwrap();
        assert_eq!(count, 16);
    }

    struct Explicit {
        count: usize,
        groups: Vec<Vec<usize>>,
    }
    impl Relation for Explicit {
        type Group<'a> = &'a [usize];
        fn num_sites(&self) -> usize {
            self.count
        }
        fn groups(&self) -> impl Iterator<Item = Self::Group<'_>> + '_ {
            self.groups.iter().map(Vec::as_slice)
        }
    }

    #[test]
    fn invalid_inputs_are_rejected_before_emitting_their_terms() {
        let sites = [HilbertSpace::factor(0usize)];
        for relation in [
            Explicit {
                count: 2,
                groups: vec![],
            },
            Explicit {
                count: 1,
                groups: vec![vec![1]],
            },
            Explicit {
                count: 1,
                groups: vec![vec![0, 0]],
            },
        ] {
            let mut emitted = false;
            assert!(
                for_each_pauli_string::<Complex64, _, _>(&relation, &sites, |_| emitted = true,)
                    .is_err()
            );
            assert!(!emitted);
        }
        let duplicate = [HilbertSpace::factor(0usize), HilbertSpace::factor(0)];
        assert!(
            for_each_pauli_string::<Complex64, _, _>(
                &SupportCombinations::<0>::new(2),
                &duplicate,
                |_| {},
            )
            .is_err()
        );
        let composite = HilbertSpace::tensor_all(sites.iter().cloned()).unwrap();
        assert!(
            for_each_pauli_string::<Complex64, _, _>(
                &SupportCombinations::<0>::new(1),
                &[composite],
                |_| {},
            )
            .is_err()
        );
    }
}
