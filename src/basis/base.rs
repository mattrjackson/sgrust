use serde::{Deserialize, Serialize};

pub trait Basis
{
    fn eval(&self, level: u32, index: u32, x: f64) -> f64;
    fn eval_deriv(&self, level: u32, index: u32, x: f64) -> f64;
    fn degree(&self) -> usize;
    fn integral(&self, level: u32, index: u32) -> f64; 
    fn basis_type(&self) -> BasisFunction;
    fn node(&self, level: u32, index: u32) -> f64;
    fn num_nodes(&self, level: u32) -> usize;
    fn get_point(&self, level: u32, index: u32) -> f64;
    fn weights(&self, level: u32) -> Vec<f64>;
    fn nodes(&self, level: u32) -> Vec<f64>;
    fn interpolation_exactness(&self, level: u32) -> u32;
    fn quadrature_exactness(&self, level: u32) -> u32;
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum BasisFunction
{
    Linear = 0,
    ClenshawCurtis = 1,
    GaussPatterson = 2,
    CustomBasis = 3,
    GaussLegendre = 4,
}
