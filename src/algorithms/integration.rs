use std::ops::AddAssign;

use num_traits::Float;

use crate::{ basis::base::Basis, storage::linear_grid::SparseGridData};



pub trait IsotropicQuadrature<T: Float + AddAssign, const D: usize, const DIM_OUT: usize>
{
    fn eval(&self, _storage: &SparseGridData<D>, _alpha: &[[T; DIM_OUT]]) -> [T; DIM_OUT]
    {
        panic!("Isotropic Quadrature not supported for this basis.");
    }
}

pub trait AnisotropicQuadrature<const D: usize, const DIM_OUT: usize>
{
    fn eval(&self, storage: &SparseGridData<D>, index: usize, dim: usize) -> [f64; DIM_OUT];
}

pub trait BasisAndQuadrature<const D: usize, const DIM_OUT: usize> : Quadrature<D, DIM_OUT> + Basis{}

pub trait Quadrature<const D: usize, const DIM_OUT: usize> : AnisotropicQuadrature<D, DIM_OUT> + AnisotropicQuadrature<D, DIM_OUT>
{
    
}