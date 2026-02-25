use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

/// Marker trait to support rkyv serialization. Stub implementation when feature not enabled.
#[cfg(not(feature="rkyv"))]
pub trait RkyvMarker{}
/// Marker trait to support rkyv serialization
#[cfg(feature="rkyv")]
pub trait RkyvMarker : rkyv::Archive + rkyv::Archive<Archived: 
        rkyv::Deserialize<Self, rkyv::rancor::Strategy<rkyv::de::Pool, rkyv::rancor::Error>> +
        rkyv::Portable +
        for<'a> rkyv::bytecheck::CheckBytes<rkyv::rancor::Strategy<rkyv::validation::Validator<rkyv::validation::archive::ArchiveValidator<'a>, rkyv::validation::shared::SharedValidator>, rkyv::rancor::Error>>
    > +
    for<'a> rkyv::Serialize<rkyv::rancor::Strategy<rkyv::ser::Serializer<rkyv::util::AlignedVec, rkyv::ser::allocator::ArenaHandle<'a>, rkyv::ser::sharing::Share>, rkyv::rancor::Error>> where Self: Sized{}
impl RkyvMarker for f32{}
impl RkyvMarker for f64{}
pub trait Float : Into<f64> + Copy + Mul<Output = Self> + Div<Output = Self> + Add<Output = Self> + Sub<Output = Self> + AddAssign + SubAssign + MulAssign + DivAssign + RkyvMarker {
    fn to_f64(&self) -> f64;
    fn from(value: f64) -> Self;
    fn zero() -> Self;
    fn exp(&self) -> Self;
    fn sqrt(&self) -> Self;
    fn ln(&self) -> Self;
    fn log10(&self) -> Self;
    fn powf(&self, power: Self) -> Self;
    fn max(&self, other: Self) -> Self;
    /// Accumulates `result[i] += alpha[i] * scale` for all i.
    /// Default implementation is scalar; f64/f32 override with SIMD via the `wide` crate.
    fn accumulate(result: &mut [Self], alpha: &[Self], scale: Self) {
        let n = result.len();
        for i in 0..n {
            result[i] += alpha[i] * scale;
        }
    }
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
    #[inline]
    fn accumulate(result: &mut [f64], alpha: &[f64], scale: f64) {
        use wide::f64x4;
        let n = result.len();
        let scale4 = f64x4::splat(scale);
        let mut i = 0;
        while i + 4 <= n {
            let a = f64x4::from([alpha[i], alpha[i + 1], alpha[i + 2], alpha[i + 3]]);
            let r = f64x4::from([result[i], result[i + 1], result[i + 2], result[i + 3]]);
            let sum: [f64; 4] = (r + a * scale4).into();
            result[i..i + 4].copy_from_slice(&sum);
            i += 4;
        }
        let rem = n - i;
        if rem > 2 {
            // rem == 3: one zero-padded SIMD op beats 3 scalar multiply-adds
            let a = f64x4::from([alpha[i], alpha[i + 1], alpha[i + 2], 0.0]);
            let r = f64x4::from([result[i], result[i + 1], result[i + 2], 0.0]);
            let sum: [f64; 4] = (r + a * scale4).into();
            result[i] = sum[0]; result[i + 1] = sum[1]; result[i + 2] = sum[2];
        } else {
            while i < n {
                result[i] += alpha[i] * scale;
                i += 1;
            }
        }
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

    fn accumulate(result: &mut [f32], alpha: &[f32], scale: f32) {
        use wide::f32x4;
        let n = result.len();
        let scale4 = f32x4::splat(scale);
        let mut i = 0;
        while i + 4 <= n {
            let a = f32x4::from([alpha[i], alpha[i + 1], alpha[i + 2], alpha[i + 3]]);
            let r = f32x4::from([result[i], result[i + 1], result[i + 2], result[i + 3]]);
            let sum: [f32; 4] = (r + a * scale4).into();
            result[i..i + 4].copy_from_slice(&sum);
            i += 4;
        }
        let rem = n - i;
        if rem > 2 {
            // rem == 3: one zero-padded SIMD op beats 3 scalar multiply-adds
            let a = f32x4::from([alpha[i], alpha[i + 1], alpha[i + 2], 0.0]);
            let r = f32x4::from([result[i], result[i + 1], result[i + 2], 0.0]);
            let sum: [f32; 4] = (r + a * scale4).into();
            result[i] = sum[0]; result[i + 1] = sum[1]; result[i + 2] = sum[2];
        } else {
            while i < n {
                result[i] += alpha[i] * scale;
                i += 1;
            }
        }
    }
}
