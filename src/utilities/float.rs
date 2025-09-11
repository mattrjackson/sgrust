use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

pub trait Float : Into<f64> + Copy + Mul<Output = Self> + Div<Output = Self> + Add<Output = Self> + Sub<Output = Self> + AddAssign + SubAssign + MulAssign + DivAssign {    
    fn to_f64(&self) -> f64;
    fn from(value: f64) -> Self;
    fn zero() -> Self;
    fn exp(&self) -> Self;
    fn sqrt(&self) -> Self;
    fn ln(&self) -> Self;
    fn log10(&self) -> Self;
    fn powf(&self, power: Self) -> Self;
    fn max(&self, other: Self) -> Self;
}
impl Float for f64{
    fn to_f64(&self) -> f64 {
        *self
    }

    fn from(value: f64) -> Self {
        value
    }
    
    fn zero() -> Self {
        0.0
    }
    
    fn sqrt(&self) -> Self {
        f64::sqrt(*self)
    }
    
    fn ln(&self) -> Self {
        f64::ln(*self)
    }
    
    fn log10(&self) -> Self {
        f64::log10(*self)
    }
    
    fn powf(&self, power: Self) -> Self {
        f64::powf(*self, power)
    }
    
    fn max(&self, other: Self) -> Self {
        f64::max(*self, other)
    }
    
    fn exp(&self) -> Self {
        f64::exp(*self)
    }
}

impl Float for f32
{
    fn to_f64(&self) -> f64 {
        *self as f64
    }

    fn from(value: f64) -> Self {
        value as f32
    }

    fn zero() -> Self {
        0.0
    }
    fn sqrt(&self) -> Self {
        f32::sqrt(*self)
    }
    
    fn ln(&self) -> Self {
        f32::ln(*self)
    }
    
    fn log10(&self) -> Self {
        f32::log10(*self)
    }
    
    fn powf(&self, power: Self) -> Self {
        f32::powf(*self, power)
    }

    fn max(&self, other: Self) -> Self {
        f32::max(*self, other)
    }
    
    fn exp(&self) -> Self {
        f32::exp(*self)
    }
}
