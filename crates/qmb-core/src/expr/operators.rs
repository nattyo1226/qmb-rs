use core::ops::{Add, Mul, Neg, Sub};
use num_traits::{One, Zero};

use super::{LocalOperator, LocalState};
use crate::{Complex, LocalStructure};

// All four ownership combinations feed the same consuming constructor.
macro_rules! binary {
    ($trait:ident,$method:ident,$lhs:ident,$rhs:ident,$output:ident,[$($bounds:tt)*],$combine:expr)=>{
        impl<L,T> $trait<$rhs<L,T>> for $lhs<L,T>
        where L:LocalStructure, $($bounds)* {
            type Output=$output<L,T>;
            fn $method(self,rhs:$rhs<L,T>)->Self::Output { ($combine)(self,rhs) }
        }
        impl<L,T> $trait<&$rhs<L,T>> for $lhs<L,T>
        where L:LocalStructure, $($bounds)* {
            type Output=$output<L,T>;
            fn $method(self,rhs:&$rhs<L,T>)->Self::Output { ($combine)(self,rhs.clone()) }
        }
        impl<L,T> $trait<$rhs<L,T>> for &$lhs<L,T>
        where L:LocalStructure, $($bounds)* {
            type Output=$output<L,T>;
            fn $method(self,rhs:$rhs<L,T>)->Self::Output { ($combine)(self.clone(),rhs) }
        }
        impl<L,T> $trait<&$rhs<L,T>> for &$lhs<L,T>
        where L:LocalStructure, $($bounds)* {
            type Output=$output<L,T>;
            fn $method(self,rhs:&$rhs<L,T>)->Self::Output { ($combine)(self.clone(),rhs.clone()) }
        }
    }
}

binary!(Add, add, LocalState, LocalState, LocalState, [], |a, b| {
    LocalState::sum([a, b])
});
binary!(
    Add,
    add,
    LocalOperator,
    LocalOperator,
    LocalOperator,
    [],
    |a, b| LocalOperator::sum([a, b])
);
binary!(
    Mul,
    mul,
    LocalOperator,
    LocalOperator,
    LocalOperator,
    [],
    |a: LocalOperator<L, T>, b| a.into_product(b)
);
binary!(
    Mul,
    mul,
    LocalOperator,
    LocalState,
    LocalState,
    [],
    |a, b: LocalState<L, T>| b.into_act(a)
);
binary!(Sub,sub,LocalState,LocalState,LocalState,[T:Zero+One+Neg<Output=T>,],|a,b:LocalState<L,T>| LocalState::sum([a,-b]));
binary!(Sub,sub,LocalOperator,LocalOperator,LocalOperator,[T:Zero+One+Neg<Output=T>,],|a,b:LocalOperator<L,T>| LocalOperator::sum([a,-b]));

macro_rules! scalar {
    ($expression:ident) => {
        impl<L: LocalStructure, T: Zero + One + Neg<Output = T>> Neg for $expression<L, T> {
            type Output = Self;
            fn neg(self) -> Self {
                self.into_scaled(Complex::new(-T::one(), T::zero()))
            }
        }
        impl<L: LocalStructure, T: Zero + One + Neg<Output = T>> Neg for &$expression<L, T> {
            type Output = $expression<L, T>;
            fn neg(self) -> Self::Output {
                -self.clone()
            }
        }
        impl<L: LocalStructure, T> Mul<Complex<T>> for $expression<L, T> {
            type Output = Self;
            fn mul(self, rhs: Complex<T>) -> Self {
                self.into_scaled(rhs)
            }
        }
        impl<L: LocalStructure, T> Mul<Complex<T>> for &$expression<L, T> {
            type Output = $expression<L, T>;
            fn mul(self, rhs: Complex<T>) -> Self::Output {
                self.clone().into_scaled(rhs)
            }
        }
        impl<L: LocalStructure, T> Mul<$expression<L, T>> for Complex<T> {
            type Output = $expression<L, T>;
            fn mul(self, rhs: $expression<L, T>) -> Self::Output {
                rhs.into_scaled(self)
            }
        }
        impl<L: LocalStructure, T> Mul<&$expression<L, T>> for Complex<T> {
            type Output = $expression<L, T>;
            fn mul(self, rhs: &$expression<L, T>) -> Self::Output {
                rhs.clone().into_scaled(self)
            }
        }
        impl<L: LocalStructure, T: Clone> Mul<&Complex<T>> for $expression<L, T> {
            type Output = Self;
            fn mul(self, rhs: &Complex<T>) -> Self {
                self.into_scaled(rhs.clone())
            }
        }
        impl<L: LocalStructure, T: Clone> Mul<&Complex<T>> for &$expression<L, T> {
            type Output = $expression<L, T>;
            fn mul(self, rhs: &Complex<T>) -> Self::Output {
                self.clone().into_scaled(rhs.clone())
            }
        }
        impl<L: LocalStructure, T: Clone> Mul<$expression<L, T>> for &Complex<T> {
            type Output = $expression<L, T>;
            fn mul(self, rhs: $expression<L, T>) -> Self::Output {
                rhs.into_scaled(self.clone())
            }
        }
        impl<L: LocalStructure, T: Clone> Mul<&$expression<L, T>> for &Complex<T> {
            type Output = $expression<L, T>;
            fn mul(self, rhs: &$expression<L, T>) -> Self::Output {
                rhs.clone().into_scaled(self.clone())
            }
        }
    };
}
scalar!(LocalState);
scalar!(LocalOperator);

impl<L: LocalStructure, T: Zero + One + Neg<Output = T>> LocalOperator<L, T> {
    /// Constructs the commutator `self * rhs - rhs * self`.
    pub fn commutator(&self, rhs: &Self) -> Self {
        self.product(rhs) - rhs.product(self)
    }
}
