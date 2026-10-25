use super::{Operator, State};
use crate::{Complex, ExpressionError, LocalStructure};
use core::ops::{Add, Mul, Neg, Sub};
use num_traits::{One, Zero};

macro_rules! binary {
    ($trait:ident,$method:ident,$lhs:ident,$rhs:ident,$output:ident,[$($bounds:tt)*],$combine:expr)=>{
        impl<I:Eq,L:LocalStructure,T> $trait<$rhs<I,L,T>> for $lhs<I,L,T>
        where $($bounds)* {
            type Output=Result<$output<I,L,T>,ExpressionError>;
            fn $method(self,rhs:$rhs<I,L,T>)->Self::Output { ($combine)(self,rhs) }
        }
        impl<I:Eq,L:LocalStructure,T> $trait<&$rhs<I,L,T>> for $lhs<I,L,T>
        where $($bounds)* {
            type Output=Result<$output<I,L,T>,ExpressionError>;
            fn $method(self,rhs:&$rhs<I,L,T>)->Self::Output { ($combine)(self,rhs.clone()) }
        }
        impl<I:Eq,L:LocalStructure,T> $trait<$rhs<I,L,T>> for &$lhs<I,L,T>
        where $($bounds)* {
            type Output=Result<$output<I,L,T>,ExpressionError>;
            fn $method(self,rhs:$rhs<I,L,T>)->Self::Output { ($combine)(self.clone(),rhs) }
        }
        impl<I:Eq,L:LocalStructure,T> $trait<&$rhs<I,L,T>> for &$lhs<I,L,T>
        where $($bounds)* {
            type Output=Result<$output<I,L,T>,ExpressionError>;
            fn $method(self,rhs:&$rhs<I,L,T>)->Self::Output { ($combine)(self.clone(),rhs.clone()) }
        }
    }
}
binary!(Add, add, State, State, State, [], |a: State<I, L, T>, b| a
    .into_add(b));
binary!(
    Add,
    add,
    Operator,
    Operator,
    Operator,
    [],
    |a: Operator<I, L, T>, b| a.into_add(b)
);
binary!(
    Mul,
    mul,
    Operator,
    Operator,
    Operator,
    [],
    |a: Operator<I, L, T>, b| a.into_product(b)
);
binary!(
    Mul,
    mul,
    Operator,
    State,
    State,
    [],
    |a, b: State<I, L, T>| b.into_act(a)
);
binary!(Sub,sub,State,State,State,[T:Zero+One+Neg<Output=T>,],|a:State<I,L,T>,b:State<I,L,T>|a.into_add(-b));
binary!(Sub,sub,Operator,Operator,Operator,[T:Zero+One+Neg<Output=T>,],|a:Operator<I,L,T>,b:Operator<I,L,T>|a.into_add(-b));

macro_rules! scalar {
    ($expression:ident) => {
        impl<I, L: LocalStructure, T: Zero + One + Neg<Output = T>> Neg for $expression<I, L, T> {
            type Output = Self;
            fn neg(self) -> Self {
                self.into_scaled(Complex::new(-T::one(), T::zero()))
            }
        }
        impl<I, L: LocalStructure, T: Zero + One + Neg<Output = T>> Neg for &$expression<I, L, T> {
            type Output = $expression<I, L, T>;
            fn neg(self) -> Self::Output {
                -self.clone()
            }
        }
        impl<I, L: LocalStructure, T> Mul<Complex<T>> for $expression<I, L, T> {
            type Output = Self;
            fn mul(self, rhs: Complex<T>) -> Self {
                self.into_scaled(rhs)
            }
        }
        impl<I, L: LocalStructure, T> Mul<Complex<T>> for &$expression<I, L, T> {
            type Output = $expression<I, L, T>;
            fn mul(self, rhs: Complex<T>) -> Self::Output {
                self.clone().into_scaled(rhs)
            }
        }
        impl<I, L: LocalStructure, T> Mul<$expression<I, L, T>> for Complex<T> {
            type Output = $expression<I, L, T>;
            fn mul(self, rhs: $expression<I, L, T>) -> Self::Output {
                rhs.into_scaled(self)
            }
        }
        impl<I, L: LocalStructure, T> Mul<&$expression<I, L, T>> for Complex<T> {
            type Output = $expression<I, L, T>;
            fn mul(self, rhs: &$expression<I, L, T>) -> Self::Output {
                rhs.clone().into_scaled(self)
            }
        }
        impl<I, L: LocalStructure, T: Clone> Mul<&Complex<T>> for $expression<I, L, T> {
            type Output = Self;
            fn mul(self, rhs: &Complex<T>) -> Self {
                self.into_scaled(rhs.clone())
            }
        }
        impl<I, L: LocalStructure, T: Clone> Mul<&Complex<T>> for &$expression<I, L, T> {
            type Output = $expression<I, L, T>;
            fn mul(self, rhs: &Complex<T>) -> Self::Output {
                self.clone().into_scaled(rhs.clone())
            }
        }
        impl<I, L: LocalStructure, T: Clone> Mul<$expression<I, L, T>> for &Complex<T> {
            type Output = $expression<I, L, T>;
            fn mul(self, rhs: $expression<I, L, T>) -> Self::Output {
                rhs.into_scaled(self.clone())
            }
        }
        impl<I, L: LocalStructure, T: Clone> Mul<&$expression<I, L, T>> for &Complex<T> {
            type Output = $expression<I, L, T>;
            fn mul(self, rhs: &$expression<I, L, T>) -> Self::Output {
                rhs.clone().into_scaled(self.clone())
            }
        }
    };
}
scalar!(State);
scalar!(Operator);

impl<I: Eq, L: LocalStructure, T: Zero + One + Neg<Output = T>> Operator<I, L, T> {
    /// Constructs the commutator `self * rhs - rhs * self` without evaluation.
    pub fn commutator(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        self.product(rhs)?.into_add(-rhs.product(self)?)
    }
}
