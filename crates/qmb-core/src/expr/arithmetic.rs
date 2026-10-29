use core::ops::{Add, Mul, Neg, Sub};
use num_traits::One;

use super::{LocalOperator, LocalState, Operator, State};
use crate::ExpressionError;
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

impl<S, O, C> Add for LocalState<S, O, C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::sum([self, rhs])
    }
}

impl<O, C> Add for LocalOperator<O, C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::sum([self, rhs])
    }
}

impl<S, O, C> Sub for LocalState<S, O, C>
where
    C: One + Neg<Output = C>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::sum([self, -rhs])
    }
}

impl<O, C> Sub for LocalOperator<O, C>
where
    C: One + Neg<Output = C>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::sum([self, -rhs])
    }
}

impl<O, C> Mul for LocalOperator<O, C> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.into_product(rhs)
    }
}

impl<S, O, C> Mul<LocalState<S, O, C>> for LocalOperator<O, C> {
    type Output = LocalState<S, O, C>;
    fn mul(self, rhs: LocalState<S, O, C>) -> Self::Output {
        rhs.into_act(self)
    }
}

impl<I: Eq, S, O, C> Add for State<I, S, O, C> {
    type Output = Result<Self, ExpressionError>;
    fn add(self, rhs: Self) -> Self::Output {
        self.into_add(rhs)
    }
}

impl<I: Eq, O, C> Add for Operator<I, O, C> {
    type Output = Result<Self, ExpressionError>;
    fn add(self, rhs: Self) -> Self::Output {
        self.into_add(rhs)
    }
}

impl<I: Eq, S, O, C> Sub for State<I, S, O, C>
where
    C: One + Neg<Output = C>,
{
    type Output = Result<Self, ExpressionError>;
    fn sub(self, rhs: Self) -> Self::Output {
        self.into_add(-rhs)
    }
}

impl<I: Eq, O, C> Sub for Operator<I, O, C>
where
    C: One + Neg<Output = C>,
{
    type Output = Result<Self, ExpressionError>;
    fn sub(self, rhs: Self) -> Self::Output {
        self.into_add(-rhs)
    }
}

impl<I: Eq, O, C> Mul for Operator<I, O, C> {
    type Output = Result<Self, ExpressionError>;
    fn mul(self, rhs: Self) -> Self::Output {
        self.into_product(rhs)
    }
}

impl<I: Eq, S, O, C> Mul<State<I, S, O, C>> for Operator<I, O, C> {
    type Output = Result<State<I, S, O, C>, ExpressionError>;
    fn mul(self, rhs: State<I, S, O, C>) -> Self::Output {
        rhs.into_act(self)
    }
}

borrowed_binary!([S, O, C], Add, add, LocalState<S, O, C>, LocalState<S, O, C>, []);
borrowed_binary!([O, C], Add, add, LocalOperator<O, C>, LocalOperator<O, C>, []);
borrowed_binary!([S, O, C], Sub, sub, LocalState<S, O, C>, LocalState<S, O, C>, []);
borrowed_binary!([O, C], Sub, sub, LocalOperator<O, C>, LocalOperator<O, C>, []);
borrowed_binary!([O, C], Mul, mul, LocalOperator<O, C>, LocalOperator<O, C>, []);
borrowed_binary!([S, O, C], Mul, mul, LocalOperator<O, C>, LocalState<S, O, C>, []);
borrowed_binary!([I, S, O, C], Add, add, State<I, S, O, C>, State<I, S, O, C>, []);
borrowed_binary!([I, O, C], Add, add, Operator<I, O, C>, Operator<I, O, C>, []);
borrowed_binary!([I, S, O, C], Sub, sub, State<I, S, O, C>, State<I, S, O, C>, []);
borrowed_binary!([I, O, C], Sub, sub, Operator<I, O, C>, Operator<I, O, C>, []);
borrowed_binary!([I, O, C], Mul, mul, Operator<I, O, C>, Operator<I, O, C>, []);
borrowed_binary!([I, S, O, C], Mul, mul, Operator<I, O, C>, State<I, S, O, C>, []);

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

scalar!([S, O, C], LocalState<S, O, C>, C, []);
scalar!([O, C], LocalOperator<O, C>, C, []);
scalar!([I, S, O, C], State<I, S, O, C>, C, []);
scalar!([I, O, C], Operator<I, O, C>, C, []);

impl<O, C: One + Neg<Output = C>> LocalOperator<O, C> {
    /// Constructs the commutator `self * rhs - rhs * self`.
    pub fn commutator(&self, rhs: &Self) -> Self {
        (self * rhs) - (rhs * self)
    }
}

impl<I: Eq, O, C: One + Neg<Output = C>> Operator<I, O, C> {
    /// Constructs the commutator `self * rhs - rhs * self` without evaluation.
    pub fn commutator(&self, rhs: &Self) -> Result<Self, ExpressionError> {
        (self * rhs)?.into_add(-(rhs * self)?)
    }
}

// A blanket scalar-left impl violates the orphan rule. Standard complex scalars
// get concrete impls; arbitrary coefficients always have `scaled(c)` available.
macro_rules! left_scalar {
    ($real:ty, [$($generic:ident),+], $expression:ty) => {
        impl<$($generic),+> Mul<$expression> for Complex<$real>
        {
            type Output = $expression;
            fn mul(self, rhs: $expression) -> Self::Output { rhs.into_scaled(self) }
        }
        borrowed_binary!([$($generic),+], Mul, mul, Complex<$real>, $expression, []);
    };
}
macro_rules! complex_left {
    ($real:ty) => {
        left_scalar!($real, [S, O], LocalState<S, O, Complex<$real>>);
        left_scalar!($real, [O], LocalOperator<O, Complex<$real>>);
        left_scalar!($real, [I, S, O], State<I, S, O, Complex<$real>>);
        left_scalar!($real, [I, O], Operator<I, O, Complex<$real>>);
    };
}
complex_left!(f32);
complex_left!(f64);
