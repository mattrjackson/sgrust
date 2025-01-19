use std::ops::AddAssign;

use num_traits::Float;

use crate::{algorithms::integration::{AnisotropicQuadrature, BasisAndQuadrature, IsotropicQuadrature, Quadrature}, storage::linear_grid::SparseGridData};

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
        if level == 0
        {        
            if index == 0 
            {
                1.0 - x
            } 
            else 
            {
                x
            }
        } 
        else 
        {
            0.0_f64.max(1.0-f64::abs((1 << level) as f64 * x - index as f64 ))
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

impl<T: Float + AddAssign, const D: usize, const DIM_OUT: usize> IsotropicQuadrature<T, D, DIM_OUT> for LinearBasis
{
    #[inline]
    fn eval(&self, storage: &SparseGridData<D>, alpha: &[[T; DIM_OUT]]) -> [T; DIM_OUT]
    {        
        let volume = T::from(storage.bounding_box.volume()).unwrap();
        let mut integral = [T::zero(); DIM_OUT];
        for (i, point) in storage.nodes().iter().enumerate()
        {
            let mut pow = T::from(2.0_f64.powf(-(point.level_sum() as f64))).unwrap();
            if !point.is_inner_point()
            {
                let mut num_boundaries = 0;
                for d in 0..D
                {
                    if point.index[d] == 0 || 2_usize.pow(point.level[d] as u32) == point.index[d] as usize
                    {
                        num_boundaries += 1;
                    }
                }
                pow = pow * T::from(2.0_f64.powi(-num_boundaries)).unwrap();
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
}



impl<const D: usize, const DIM_OUT: usize> AnisotropicQuadrature<D, DIM_OUT> for LinearBasis
{
    #[inline]
    fn eval(&self, storage: &SparseGridData<D>, index: usize, dim: usize) -> [f64; DIM_OUT]
    {   
        let mut integral_component = [1.0; DIM_OUT];
        let point = &storage[index];
        let mut weight = 2.0_f64.powf(-(point.level[dim] as f64));
        if !point.is_inner_point()
        {
            let mut num_boundaries = 0;
            if point.index[dim] == 0 || 2_usize.pow(point.level[dim] as u32) == point.index[dim] as usize
            {
                num_boundaries += 1;
            }
            weight *= 2.0_f64.powi(-num_boundaries);
        }        
        #[allow(clippy::needless_range_loop)]
        for i in 0..DIM_OUT
        {
            integral_component[i] = weight;
        }
        integral_component
    }
}

impl<const D: usize, const DIM_OUT: usize>  Quadrature<D, DIM_OUT> for LinearBasis{}
impl<const D: usize, const DIM_OUT: usize>  BasisAndQuadrature<D, DIM_OUT> for LinearBasis{}