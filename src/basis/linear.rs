use std::ops::AddAssign;

use num_traits::Float;

use crate::const_generic::storage::SparseGridData;

use super::base::Basis;

#[derive(Copy, Clone)]
pub struct LinearBasis;

impl Basis for LinearBasis
{
    fn quadrature_exactness(&self, level: u32) -> u32
    {
        if level == 0
        {
            1
        }
        else
        {
            2_u32.pow(level) - 1
        }
    }

    fn interpolation_exactness(&self, level: u32) -> u32
    {
        if level == 0
        {
            0
        }
        else
        {
            2_u32.pow(level) - 1
        }
    }

    #[inline]
    fn node(&self, level: u32, index: u32) -> f64
    {
        index as f64 / (1 << level) as f64        
    }
    #[inline]
    fn get_point(&self, level: u32, index: u32) -> f64
    {
        self.node(level, index)
    }
    #[inline]
    fn eval(&self, level: u32, index: u32, x: f64) -> f64 {
        if level == 0 {
            // Avoid branching on index by using a simple formula
            // index is 0 or 1, so: if 0 => 1.0 - x, if 1 => x
            (1.0 - x) * (1 - index) as f64 + x * index as f64
        } else {
            // Avoid repeated computation and function calls
            let scale = (1 << level) as f64;
            let center = index as f64 / scale;
            let dist = (x - center).abs();
            (1.0 - dist * scale).max(0.0)
        }
    }

    fn eval_deriv(&self, _level: u32, _index: u32, _x: f64) -> f64 {
        panic!("LinearBasis does not implement eval_deriv")
    }

    fn degree(&self) -> usize {
        1
    }
    #[inline]
    fn integral(&self, level: u32, _index: u32) -> f64 {
        if level == 0
        {
            0.5
        }
        else
        {
            1.0 / (1 << level) as f64
        }
    }
    
    fn basis_type(&self) -> super::base::BasisFunction {
        super::base::BasisFunction::Linear
    }
    #[inline]
    fn num_nodes(&self, level: u32) -> usize
    {
        (1 << level) as usize
    }
    #[inline]
    fn weights(&self, level: u32) -> Vec<f64> 
    {
        if level == 0
        {
            vec![0.5]
        }
        else
        {
            let n = self.num_nodes(level);
            vec![ 1.0 / (n - 1) as f64; self.num_nodes(level)]
        }
    }
    #[inline]
    fn nodes(&self, level: u32) -> Vec<f64> {
        if level == 0
        {
            vec![0.5]
        }
        else
        {
            let n = self.num_nodes(level);
            let mut r = vec![0.0; n];
            #[allow(clippy::needless_range_loop)]
            for index in 0..n
            {
                r[index] = index as f64 / (1 << level) as f64;
            }
            r
        }
    }
   
    
}

impl LinearBasis
{
    #[inline]
    pub fn eval_point<const D: usize, const DIM_OUT: usize, T: Float + AddAssign>(storage: &SparseGridData<D>, alpha: &[[T; DIM_OUT]]) -> [T; DIM_OUT]
    {        
        let volume = T::from(storage.bounding_box.volume()).unwrap();
        let mut integral = [T::zero(); DIM_OUT];
        for (i, point) in storage.nodes().iter().enumerate()
        {
            let mut pow = T::from( 1.0 / (1 << point.level_sum())  as f64).unwrap();
            if !point.is_inner_point()
            {
                let mut num_boundaries = 0;
                for d in 0..D
                {
                    if point.index[d] == 0 || (1 << point.level[d] as usize) == point.index[d] as usize
                    {
                        num_boundaries += 1;
                    }
                }
                pow = pow * T::from(1.0 / (1 << num_boundaries) as f64).unwrap();
            }                
            #[allow(clippy::needless_range_loop)]
            for dim in 0..DIM_OUT
            {
                integral[dim] += alpha[i][dim] * pow;
            }
        }
        #[allow(clippy::needless_range_loop)]
        for dim in 0..DIM_OUT
        {
            integral[dim] = integral[dim] * volume;
        }
        integral
    }

    #[inline]
    pub fn eval_point_dynamic<T: Float + AddAssign>(storage: &crate::dynamic::storage::SparseGridData, alpha: &[T], integral: &mut [T])
    {        
        integral.fill(T::zero());
        let volume = T::from(storage.bounding_box.volume()).unwrap();
        for (i, point) in storage.nodes().enumerate()
        {
            let mut pow = T::from( 1.0 / (1 << point.level_sum())  as f64).unwrap();
            if !point.flags.is_inner()
            {
                let mut num_boundaries = 0;
                for d in 0..storage.num_inputs
                {
                    if point.index[d] == 0 || (1 << point.level[d] as usize) == point.index[d] as usize
                    {
                        num_boundaries += 1;
                    }
                }
                pow = pow * T::from(1.0 / (1 << num_boundaries) as f64).unwrap();
            }                
            #[allow(clippy::needless_range_loop)]
            for dim in 0..storage.num_outputs
            {
                integral[dim] += alpha[i*storage.num_outputs+dim] * pow;
            }
        }
        #[allow(clippy::needless_range_loop)]
        for dim in 0..storage.num_outputs
        {
            integral[dim] = integral[dim] * volume;
        }
    }
}