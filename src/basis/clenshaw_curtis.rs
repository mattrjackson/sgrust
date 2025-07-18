use static_init::dynamic;

use super::base::Basis;

static CC_MAX_LEVEL: u32 = 12;
use crate::tables::clenshaw_curtis_table::*;

#[derive(Clone)]
pub struct ClenshawCurtisCache{ nodes: Vec<Vec<f64>>, weights: Vec<Vec<f64>>}

impl ClenshawCurtisCache
{
    pub fn new(max_level: u32) -> Self
    {        
        let mut nodes = Vec::new();
        let mut weights = Vec::new();
        for level in 0..=max_level
        {
            nodes.push(cc_nodes(level));
            weights.push(cc_weights(level));
        }        
        Self{ nodes, weights}
    }
    pub(crate) fn node(&self, level: u32, index: u32) -> f64
    {
        if level <= CC_MAX_LEVEL
        {
            CC_CACHE.nodes[level as usize][index as usize]
        }
        else
        {
            cc_nodes(level)[index as usize]
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
        if level <= CC_MAX_LEVEL
        {
            CC_CACHE.weights[level as usize][index as usize]
        }
        else
        {
            cc_weights(level)[index as usize]
        }
    }
    pub fn weights_for_level(&self, level: u32) -> Vec<f64>
    {
        if level <= CC_MAX_LEVEL
        {
            CC_CACHE.weights[level as usize].clone()
        }
        else
        {
            cc_weights(level)
        }
    }

    pub fn nodes_for_level(&self, level: u32) -> Vec<f64>
    {
        if level <= CC_MAX_LEVEL
        {
            CC_CACHE.nodes[level as usize].clone()
        }
        else
        {
            cc_nodes(level)            
        }
    }
}

#[dynamic]
pub(crate) static CC_CACHE: ClenshawCurtisCache = ClenshawCurtisCache::new(CC_MAX_LEVEL) ;

#[derive(Clone, Copy)]
pub struct ClenshawCurtis;

impl Basis for ClenshawCurtis
{
    fn eval(&self, level: u32, index: u32, x: f64) -> f64 {
        let xl = CC_CACHE.node(level, index  - 1);        
        let xr = CC_CACHE.node(level, index + 1);
        let xm = CC_CACHE.node(level, index);
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
        CC_CACHE.node(level, index)
    }

    fn eval_deriv(&self, _level: u32, _index: u32, _x: f64) -> f64 {
        panic!("Derivative not implemented for CurtisClenshaw Basis.")
    }

    fn degree(&self) -> usize {
        1
    }

    fn integral(&self, level: u32, index: u32) -> f64 {
        CC_CACHE.weight(level, index)
    }

    fn basis_type(&self) -> super::base::BasisFunction {
        super::base::BasisFunction::ClenshawCurtis
    }

    fn num_nodes(&self, level: u32) -> usize
    {
        cc_num_nodes(level) as usize 
    }
    
    fn node(&self, level: u32, index: u32) -> f64 {
        CC_CACHE.node(level, index)
    }
    
    fn weights(&self, level: u32) -> Vec<f64> {
        CC_CACHE.weights_for_level(level)
    }
    
    fn nodes(&self, level: u32) -> Vec<f64> {
        CC_CACHE.nodes_for_level(level)
    }

    fn quadrature_exactness(&self, level: u32) -> u32
    {
        if level == 0
        {
            1
        }
        else
        {
            2_u32.pow(level) + 1
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
            2_u32.pow(level)
        }
    }

}
