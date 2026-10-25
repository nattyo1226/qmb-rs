use core::ops::{Add, Mul, Neg, Sub};
use num_traits::{One, Zero};

use super::{LocalOperator, LocalState, Operator, State};
use crate::{Complex, ExpressionError, LocalStructure};

// Forward the three borrowed forms to the owned implementation. Only borrowed
// operands need Clone; the owned implementation determines Output and bounds.
macro_rules! borrowed_binary {
    ([$($generic:ident),+], $trait:ident, $method:ident, $lhs:ty, $rhs:ty, [$($bounds:tt)*]) => {
        impl<$($generic),+> $trait<&$rhs> for $lhs
        where
            $($bounds)*
            $lhs: $trait<$rhs>,
            $rhs: Clone,
        {
            type Output = <$lhs as $trait<$rhs>>::Output;
            fn $method(self, rhs: &$rhs) -> Self::Output {
                <$lhs as $trait<$rhs>>::$method(self, <$rhs as Clone>::clone(rhs))
            }
        }
        impl<$($generic),+> $trait<$rhs> for &$lhs
        where
            $($bounds)*
            $lhs: $trait<$rhs> + Clone,
        {
            type Output = <$lhs as $trait<$rhs>>::Output;
            fn $method(self, rhs: $rhs) -> Self::Output {
                <$lhs as $trait<$rhs>>::$method(<$lhs as Clone>::clone(self), rhs)
            }
        }
        impl<$($generic),+> $trait<&$rhs> for &$lhs
        where
            $($bounds)*
            $lhs: $trait<$rhs> + Clone,
            $rhs: Clone,
        {
            type Output = <$lhs as $trait<$rhs>>::Output;
            fn $method(self, rhs: &$rhs) -> Self::Output {
                <$lhs as $trait<$rhs>>::$method(
                    <$lhs as Clone>::clone(self),
                    <$rhs as Clone>::clone(rhs),
                )
            }
        }
    };
}

impl<L: LocalStructure, T> Add for LocalState<L, T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::sum([self, rhs])
    }
}

impl<L: LocalStructure, T> Add for LocalOperator<L, T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::sum([self, rhs])
    }
}

impl<L: LocalStructure, T> Sub for LocalState<L, T>
where
    T: Zero + One + Neg<Output = T>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::sum([self, -rhs])
    }
}

impl<L: LocalStructure, T> Sub for LocalOperator<L, T>
where
    T: Zero + One + Neg<Output = T>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::sum([self, -rhs])
    }
}

impl<L: LocalStructure, T> Mul for LocalOperator<L, T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.into_product(rhs)
    }
}

impl<L: LocalStructure, T> Mul<LocalState<L, T>> for LocalOperator<L, T> {
    type Output = LocalState<L, T>;
    fn mul(self, rhs: LocalState<L, T>) -> Self::Output {
        rhs.into_act(self)
    }
}

impl<I: Eq, L: LocalStructure, T> Add for State<I, L, T> {
    type Output = Result<Self, ExpressionError>;
    fn add(self, rhs: Self) -> Self::Output {
        self.into_add(rhs)
    }
}

impl<I: Eq, L: LocalStructure, T> Add for Operator<I, L, T> {
    type Output = Result<Self, ExpressionError>;
    fn add(self, rhs: Self) -> Self::Output {
        self.into_add(rhs)
    }
}

impl<I: Eq, L: LocalStructure, T> Sub for State<I, L, T>
where
    T: Zero + One + Neg<Output = T>,
{
    type Output = Result<Self, ExpressionError>;
    fn sub(self, rhs: Self) -> Self::Output {
        self.into_add(-rhs)
    }
}

impl<I: Eq, L: LocalStructure, T> Sub for Operator<I, L, T>
where
    T: Zero + One + Neg<Output = T>,
{
    type Output = Result<Self, ExpressionError>;
    fn sub(self, rhs: Self) -> Self::Output {
        self.into_add(-rhs)
    }
}

impl<I: Eq, L: LocalStructure, T> Mul for Operator<I, L, T> {
    type Output = Result<Self, ExpressionError>;
    fn mul(self, rhs: Self) -> Self::Output {
        self.into_product(rhs)
    }
}

impl<I: Eq, L: LocalStructure, T> Mul<State<I, L, T>> for Operator<I, L, T> {
    type Output = Result<State<I, L, T>, ExpressionError>;
    fn mul(self, rhs: State<I, L, T>) -> Self::Output {
        rhs.into_act(self)
    }
}

borrowed_binary!([L, T], Add, add, LocalState<L, T>, LocalState<L, T>, [L: LocalStructure,]);
borrowed_binary!([L, T], Add, add, LocalOperator<L, T>, LocalOperator<L, T>, [L: LocalStructure,]);
borrowed_binary!([L, T], Sub, sub, LocalState<L, T>, LocalState<L, T>, [L: LocalStructure,]);
borrowed_binary!([L, T], Sub, sub, LocalOperator<L, T>, LocalOperator<L, T>, [L: LocalStructure,]);
borrowed_binary!([L, T], Mul, mul, LocalOperator<L, T>, LocalOperator<L, T>, [L: LocalStructure,]);
borrowed_binary!([L, T], Mul, mul, LocalOperator<L, T>, LocalState<L, T>, [L: LocalStructure,]);
borrowed_binary!([I, L, T], Add, add, State<I, L, T>, State<I, L, T>, [L: LocalStructure,]);
borrowed_binary!([I, L, T], Add, add, Operator<I, L, T>, Operator<I, L, T>, [L: LocalStructure,]);
borrowed_binary!([I, L, T], Sub, sub, State<I, L, T>, State<I, L, T>, [L: LocalStructure,]);
borrowed_binary!([I, L, T], Sub, sub, Operator<I, L, T>, Operator<I, L, T>, [L: LocalStructure,]);
borrowed_binary!([I, L, T], Mul, mul, Operator<I, L, T>, Operator<I, L, T>, [L: LocalStructure,]);
borrowed_binary!([I, L, T], Mul, mul, Operator<I, L, T>, State<I, L, T>, [L: LocalStructure,]);

// Only scalar construction lives here. Borrowed multiplication uses the same
// forwarding macro as expression arithmetic, without a separate ownership policy.
macro_rules! scalar {
    ([$($generic:ident),+], $expression:ty, $scalar:ident, [$($bounds:tt)*]) => {
        impl<$($generic),+> Neg for $expression
        where
            $($bounds)*
            $scalar: Zero + One + Neg<Output = $scalar>,
        {
            type Output = Self;
            fn neg(self) -> Self {
                self.into_scaled(Complex::new(-$scalar::one(), $scalar::zero()))
            }
        }
        impl<$($generic),+> Neg for &$expression
        where
            $($bounds)*
            $expression: Neg<Output = $expression> + Clone,
        {
            type Output = $expression;
            fn neg(self) -> Self::Output {
                <$expression as Neg>::neg(<$expression as Clone>::clone(self))
            }
        }
        impl<$($generic),+> Mul<Complex<$scalar>> for $expression
        where $($bounds)*
        {
            type Output = Self;
            fn mul(self, rhs: Complex<$scalar>) -> Self {
                self.into_scaled(rhs)
            }
        }
        impl<$($generic),+> Mul<$expression> for Complex<$scalar>
        where $($bounds)*
        {
            type Output = $expression;
            fn mul(self, rhs: $expression) -> Self::Output {
                rhs.into_scaled(self)
            }
        }
        borrowed_binary!([$($generic),+], Mul, mul, $expression, Complex<$scalar>, [$($bounds)*]);
        borrowed_binary!([$($generic),+], Mul, mul, Complex<$scalar>, $expression, [$($bounds)*]);
    };
}

scalar!([L, T], LocalState<L, T>, T, [L: LocalStructure,]);
scalar!([L, T], LocalOperator<L, T>, T, [L: LocalStructure,]);
scalar!([I, L, T], State<I, L, T>, T, [L: LocalStructure,]);
scalar!([I, L, T], Operator<I, L, T>, T, [L: LocalStructure,]);

impl<L: LocalStructure, T: Zero + One + Neg<Output = T>> LocalOperator<L, T> {
    /// Constructs the commutator `self * rhs - rhs * self`.
    pub fn commutator(&self, rhs: &Self) -> Self {
        self.product(rhs) - rhs.product(self)
    }
}

impl<I: Eq, L: LocalStructure, T: Zero + One + Neg<Output = T>> Operator<I, L, T> {
    /// Constructs the commutator `self * rhs - rhs * self` without evaluation.
    pub fn commutator(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        self.product(rhs)?.into_add(-rhs.product(self)?)
    }
}
