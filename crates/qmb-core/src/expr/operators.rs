use core::ops::{Add, Mul, Neg, Sub};
use num_traits::One;

use super::{LocalOperator, LocalState, Operator, State};
use crate::{ExpressionError, LocalStructure};
use num_complex::Complex;

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

impl<L: LocalStructure, C> Add for LocalState<L, C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::sum([self, rhs])
    }
}

impl<L: LocalStructure, C> Add for LocalOperator<L, C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::sum([self, rhs])
    }
}

impl<L: LocalStructure, C> Sub for LocalState<L, C>
where
    C: One + Neg<Output = C>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::sum([self, -rhs])
    }
}

impl<L: LocalStructure, C> Sub for LocalOperator<L, C>
where
    C: One + Neg<Output = C>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::sum([self, -rhs])
    }
}

impl<L: LocalStructure, C> Mul for LocalOperator<L, C> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.into_product(rhs)
    }
}

impl<L: LocalStructure, C> Mul<LocalState<L, C>> for LocalOperator<L, C> {
    type Output = LocalState<L, C>;
    fn mul(self, rhs: LocalState<L, C>) -> Self::Output {
        rhs.into_act(self)
    }
}

impl<I: Eq, L: LocalStructure, C> Add for State<I, L, C> {
    type Output = Result<Self, ExpressionError>;
    fn add(self, rhs: Self) -> Self::Output {
        self.into_add(rhs)
    }
}

impl<I: Eq, L: LocalStructure, C> Add for Operator<I, L, C> {
    type Output = Result<Self, ExpressionError>;
    fn add(self, rhs: Self) -> Self::Output {
        self.into_add(rhs)
    }
}

impl<I: Eq, L: LocalStructure, C> Sub for State<I, L, C>
where
    C: One + Neg<Output = C>,
{
    type Output = Result<Self, ExpressionError>;
    fn sub(self, rhs: Self) -> Self::Output {
        self.into_add(-rhs)
    }
}

impl<I: Eq, L: LocalStructure, C> Sub for Operator<I, L, C>
where
    C: One + Neg<Output = C>,
{
    type Output = Result<Self, ExpressionError>;
    fn sub(self, rhs: Self) -> Self::Output {
        self.into_add(-rhs)
    }
}

impl<I: Eq, L: LocalStructure, C> Mul for Operator<I, L, C> {
    type Output = Result<Self, ExpressionError>;
    fn mul(self, rhs: Self) -> Self::Output {
        self.into_product(rhs)
    }
}

impl<I: Eq, L: LocalStructure, C> Mul<State<I, L, C>> for Operator<I, L, C> {
    type Output = Result<State<I, L, C>, ExpressionError>;
    fn mul(self, rhs: State<I, L, C>) -> Self::Output {
        rhs.into_act(self)
    }
}

borrowed_binary!([L, C], Add, add, LocalState<L, C>, LocalState<L, C>, [L: LocalStructure,]);
borrowed_binary!([L, C], Add, add, LocalOperator<L, C>, LocalOperator<L, C>, [L: LocalStructure,]);
borrowed_binary!([L, C], Sub, sub, LocalState<L, C>, LocalState<L, C>, [L: LocalStructure,]);
borrowed_binary!([L, C], Sub, sub, LocalOperator<L, C>, LocalOperator<L, C>, [L: LocalStructure,]);
borrowed_binary!([L, C], Mul, mul, LocalOperator<L, C>, LocalOperator<L, C>, [L: LocalStructure,]);
borrowed_binary!([L, C], Mul, mul, LocalOperator<L, C>, LocalState<L, C>, [L: LocalStructure,]);
borrowed_binary!([I, L, C], Add, add, State<I, L, C>, State<I, L, C>, [L: LocalStructure,]);
borrowed_binary!([I, L, C], Add, add, Operator<I, L, C>, Operator<I, L, C>, [L: LocalStructure,]);
borrowed_binary!([I, L, C], Sub, sub, State<I, L, C>, State<I, L, C>, [L: LocalStructure,]);
borrowed_binary!([I, L, C], Sub, sub, Operator<I, L, C>, Operator<I, L, C>, [L: LocalStructure,]);
borrowed_binary!([I, L, C], Mul, mul, Operator<I, L, C>, Operator<I, L, C>, [L: LocalStructure,]);
borrowed_binary!([I, L, C], Mul, mul, Operator<I, L, C>, State<I, L, C>, [L: LocalStructure,]);

// Only scalar construction lives here. Borrowed multiplication uses the same
// forwarding macro as expression arithmetic, without a separate ownership policy.
macro_rules! scalar {
    ([$($generic:ident),+], $expression:ty, $scalar:ident, [$($bounds:tt)*]) => {
        impl<$($generic),+> Neg for $expression
        where
            $($bounds)*
            $scalar: One + Neg<Output = $scalar>,
        {
            type Output = Self;
            fn neg(self) -> Self {
                self.into_scaled(-$scalar::one())
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
        impl<$($generic),+> Mul<$scalar> for $expression
        where $($bounds)*
        {
            type Output = Self;
            fn mul(self, rhs: $scalar) -> Self { self.into_scaled(rhs) }
        }
        borrowed_binary!([$($generic),+], Mul, mul, $expression, $scalar, [$($bounds)*]);
    };
}

scalar!([L, C], LocalState<L, C>, C, [L: LocalStructure,]);
scalar!([L, C], LocalOperator<L, C>, C, [L: LocalStructure,]);
scalar!([I, L, C], State<I, L, C>, C, [L: LocalStructure,]);
scalar!([I, L, C], Operator<I, L, C>, C, [L: LocalStructure,]);

impl<L: LocalStructure, C: One + Neg<Output = C>> LocalOperator<L, C> {
    /// Constructs the commutator `self * rhs - rhs * self`.
    pub fn commutator(&self, rhs: &Self) -> Self {
        self.product(rhs) - rhs.product(self)
    }
}

impl<I: Eq, L: LocalStructure, C: One + Neg<Output = C>> Operator<I, L, C> {
    /// Constructs the commutator `self * rhs - rhs * self` without evaluation.
    pub fn commutator(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        self.product(rhs)?.into_add(-rhs.product(self)?)
    }
}

// A blanket scalar-left impl violates the orphan rule. Standard complex scalars
// get concrete impls; arbitrary coefficients always have `scaled(c)` available.
macro_rules! left_scalar {
    ($real:ty, [$($generic:ident),+], $expression:ty) => {
        impl<$($generic),+> Mul<$expression> for Complex<$real>
        where L: LocalStructure,
        {
            type Output = $expression;
            fn mul(self, rhs: $expression) -> Self::Output { rhs.into_scaled(self) }
        }
        borrowed_binary!([$($generic),+], Mul, mul, Complex<$real>, $expression, [L: LocalStructure,]);
    };
}
macro_rules! complex_left {
    ($real:ty) => {
        left_scalar!($real, [L], LocalState<L, Complex<$real>>);
        left_scalar!($real, [L], LocalOperator<L, Complex<$real>>);
        left_scalar!($real, [I, L], State<I, L, Complex<$real>>);
        left_scalar!($real, [I, L], Operator<I, L, Complex<$real>>);
    };
}
complex_left!(f32);
complex_left!(f64);
