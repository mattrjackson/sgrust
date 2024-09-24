use core::f64;
use f64::consts::PI;
use crate::{rules::SparseGridRule, tables};

pub struct OneDimensionalNodes
{
    pub x: Vec<f64>,
    pub weights: Vec<f64>,    
}
impl OneDimensionalNodes
{
    pub fn new(num_points: usize) -> Self
    {
        Self { x: vec![0.0; num_points], weights: vec![0.0; num_points] }
    }
}

#[inline]
fn pow2(level: u32) -> u32
{
    1 << level
}

pub fn get_num_points(level: u32, rule: SparseGridRule) -> u32
{
    match rule
    {        
        SparseGridRule::GaussPatterson => pow2(level+1) - 1,
        SparseGridRule::ClenshawCurtis => if level == 0 {1} else { pow2(level) + 1 }
    }
}

pub trait OneDimensionalRule
{
    fn num_nodes(&self, level: u32) -> usize;
    fn nodes(&self, level: u32) -> Vec<f64>;
    fn node(&self, level: u32, index: u32) -> f64;
    fn weights(&self, level: u32) -> Vec<f64>;
    fn weight(&self, level: u32, index: u32) -> f64;
    
}

pub struct ClenshawCurtis;

impl ClenshawCurtis
{
    #[inline(always)]
    fn calculate_node(h: f64, index: u32) -> f64
    {
        0.5*(f64::cos(PI * (1.0 - index as f64 * h)) + 1.0)
    }
    
}

impl OneDimensionalRule for ClenshawCurtis
{
    fn num_nodes(&self, level: u32) -> usize {
         if level == 0 {1} else { (1 << level) + 1 }
    }
    fn nodes(&self, level: u32) -> Vec<f64> 
    {
        clenshaw_curtis_nodes(level)
    }

    fn node(&self, level: u32, index: u32) -> f64 {
        let h = 2.0_f64.powi(level as i32);
        Self::calculate_node(h, index)
    }
    fn weights(&self, level: u32) -> Vec<f64>
    {
        tables::clenshaw_curtis_table::cc_weights(level)
    }   

    // returns weight for (0,1) interval
    fn weight(&self, level: u32, index: u32) -> f64 {
        clenshaw_curtis_weight(level, index) * 0.5
    }
    
}

/// Return Clenshaw-Curtis nodes over a (0,1) interval
pub fn clenshaw_curtis_nodes(level: u32) -> Vec<f64>
{
    let n = get_num_points(level, SparseGridRule::ClenshawCurtis) ;
    let mut nodes = vec![0.0; n as usize];    
    if level > 0
    {
        #[allow(clippy::needless_range_loop)]
        for i in 0..n as usize
        {
            nodes[i] = 0.5*(1.0 + f64::cos(PI*(n as usize - 1 - i) as f64 / (n - 1) as f64));
        }
    }
    else 
    {
        nodes[0] = 0.5;
    }
    nodes
}

/// Return Clenshaw-Curtis weight over a (0,1) interval
pub fn clenshaw_curtis_weight(level: u32, point: u32) -> f64
{
    crate::tables::clenshaw_curtis_table::cc_weights(level)[point as usize]
}

#[test]
fn test_clenshaw_curtis()
{
    let level =2;
    let nodes = clenshaw_curtis_nodes(level);
    #[allow(clippy::needless_range_loop)]
    for i in 0..nodes.len()
    {
        println!("{},{}",nodes[i], clenshaw_curtis_weight(level, i as u32));
    }
}