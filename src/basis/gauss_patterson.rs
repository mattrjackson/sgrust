use std::ops::AddAssign;

use num_traits::Float;
use static_init::dynamic;

use crate::{algorithms::integration::{AnisotropicQuadrature, BasisAndQuadrature, IsotropicQuadrature, Quadrature}, tables};

use super::base::Basis;

static GP_MAX_LEVEL: u32 = 8;

#[derive(Clone)]
pub struct GaussPattersonCache{ nodes: Vec<Vec<f64>>, weights: Vec<Vec<f64>>}

impl GaussPattersonCache
{
    pub fn new(max_level: u32) -> Self
    {
        use tables::gauss_patterson_table::*;        
        let mut nodes = Vec::new();
        let mut weights = Vec::new();
        for level in 0..=max_level
        {
            nodes.push(gauss_patterson_abscissa(level));
            weights.push(gauss_patterson_weights(level));
        }        
        Self{ nodes, weights}
    }
    pub(crate) fn node(&self, level: u32, index: u32) -> f64
    {
        if level <= GP_MAX_LEVEL
        {
            GP_CACHE.nodes[level as usize][index as usize]
        }
        else
        {
            panic!("Gauss Patterson only supported for level <= 8!")
        }
    }
    pub fn len(&self) -> usize
    {
        self.nodes.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    pub fn weight(&self, level: u32, index: u32) -> f64
    {
        if level <= GP_MAX_LEVEL
        {
            GP_CACHE.weights[level as usize][index as usize]
        }
        else
        {
            panic!("Gauss Patterson only supported for level <= 8!")
        }
    }
    pub fn weights_for_level(&self, level: u32) -> Vec<f64>
    {
        if level <= GP_MAX_LEVEL
        {
            GP_CACHE.weights[level as usize].clone()
        }
        else
        {
            panic!("Gauss Patterson only supported for level <= 8!")
        }
    }

    pub fn nodes_for_level(&self, level: u32) -> Vec<f64>
    {
        if level <= GP_MAX_LEVEL
        {
            GP_CACHE.nodes[level as usize].clone()
        }
        else
        {
            panic!("Gauss Patterson only supported for level <= 8!")       
        }
    }
}

#[dynamic]
static GP_CACHE: GaussPattersonCache = GaussPattersonCache::new(GP_MAX_LEVEL) ;

#[derive(Clone, Copy)]
pub struct GaussPatterson;

impl Basis for GaussPatterson
{
    fn eval(&self, level: u32, index: u32, x: f64) -> f64 {
        let xl = GP_CACHE.node(level, index  - 1);        
        let xr = GP_CACHE.node(level, index + 1);
        let xm = GP_CACHE.node(level, index);
        if x < xm
        {
            f64::max(0.0, 1.0-(xm-x)/(xm-xl))
        }
        else
        {
            f64::max(0.0, (xr-x)/(xr-xm))
        }
    }
    fn get_point(&self, level: u32, index: u32) -> f64
    {
        GP_CACHE.node(level, index)
    }

    fn eval_deriv(&self, _level: u32, _index: u32, _x: f64) -> f64 {
        panic!("Derivative not implemented for GaussPatterson Basis.")
    }

    fn degree(&self) -> usize {
        1
    }

    fn integral(&self, level: u32, index: u32) -> f64 {
        GP_CACHE.weight(level, index)
    }

    fn basis_type(&self) -> super::base::BasisFunction {
        super::base::BasisFunction::GaussPatterson
    }

    fn num_nodes(&self, level: u32) -> usize
    {
        crate::tables::gauss_patterson_table::num_nodes(level)         
    }
    
    fn node(&self, level: u32, index: u32) -> f64 {
        GP_CACHE.node(level, index)
    }
    
    fn weights(&self, level: u32) -> Vec<f64> {
        GP_CACHE.weights_for_level(level)
    }
    
    fn nodes(&self, level: u32) -> Vec<f64> {
        GP_CACHE.nodes_for_level(level)
    }
    fn quadrature_exactness(&self, level: u32) -> u32
    {
        if level == 0
        {
            1
        }
        else
        {
            3*2_u32.pow(level) - 1
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
            2_u32.pow(level+1) - 2
        }
    }

}

impl<const D: usize, const DIM_OUT: usize> AnisotropicQuadrature<D, DIM_OUT> for GaussPatterson
{
    fn eval(&self, storage: &crate::storage::linear_grid::SparseGridData<D>, index: usize, dim: usize) -> [f64; DIM_OUT] {
        
        let mut integral_component = [1.0; DIM_OUT];
        let point = &storage[index];
        let weight = self.integral(point.level[dim] as u32, point.index[dim]);        
        #[allow(clippy::needless_range_loop)]
        for d in 0..DIM_OUT
        {
            integral_component[d] = weight;
        }
        integral_component
    }
}

impl<T: Float + AddAssign, const D: usize, const DIM_OUT: usize>  IsotropicQuadrature<T, D, DIM_OUT> for GaussPatterson
{}

impl<const D: usize, const DIM_OUT: usize>  Quadrature<D, DIM_OUT> for GaussPatterson{}

impl<const D: usize, const DIM_OUT: usize> BasisAndQuadrature<D, DIM_OUT> for GaussPatterson{}

#[test]
fn print_gauss_patterson_nodes_and_weights()
{
    let cc = GaussPatterson;
    
    for level in 0..=GP_MAX_LEVEL
    {
        let weights = crate::tables::gauss_patterson_table::gauss_patterson_weights(level);
        for i in 0..cc.num_nodes(level) as u32
        {
            let integral = cc.integral(level, i);
            let weight = weights[i as usize];
            assert_eq!(integral, weight);
            // ensure we have the expected (0,1) interval
            assert!(integral >= 0.0);
            assert!(integral <= 1.0)
        }
    } 
}