use crate::{ basis::base::Basis, storage::linear_grid::SparseGridStorage};



pub trait IsotropicQuadrature<const D: usize, const DIM_OUT: usize>
{
    fn eval(&self, _storage: &SparseGridStorage<D>, _alpha: &[[f64; DIM_OUT]]) -> [f64; DIM_OUT]
    {
        panic!("Isotropic Quadrature not supported for this basis.");
    }
}

pub trait AnisotropicQuadrature<const D: usize, const DIM_OUT: usize>
{
    fn eval(&self, storage: &SparseGridStorage<D>, index: usize, dim: usize) -> [f64; DIM_OUT];
}

pub trait BasisAndQuadrature<const D: usize, const DIM_OUT: usize> : Quadrature<D, DIM_OUT> + Basis{}

pub trait Quadrature<const D: usize, const DIM_OUT: usize> : AnisotropicQuadrature<D, DIM_OUT> + AnisotropicQuadrature<D, DIM_OUT>
{
    
}