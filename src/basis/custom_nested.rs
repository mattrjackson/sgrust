use std::ops::{Deref, Index};

use serde::{Deserialize, Serialize};

use super::base::Basis;
#[derive(Clone, Serialize, Deserialize, Default)]
pub struct NodesForLevel(Vec<f64>);

impl Index<u32> for NodesForLevel
{
    type Output=f64;

    fn index(&self, index: u32) -> &Self::Output {
        &self.0[index as usize]
    }
}
impl Deref for NodesForLevel
{
    type Target=Vec<f64>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Clone, Serialize, Deserialize, Default)]
pub struct WeightsForLevel(Vec<f64>);

impl Index<u32> for WeightsForLevel
{
    type Output=f64;

    fn index(&self, index: u32) -> &Self::Output {
        &self.0[index as usize]
    }
}

impl Deref for WeightsForLevel
{
    type Target=Vec<f64>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Clone, Serialize, Deserialize, Default)]
pub struct CustomNested
{
    nodes: Vec<NodesForLevel>,
    weights: Vec<WeightsForLevel>,
    quadrature_exactness: Vec<u32>,
    interpolation_exactness: Vec<u32>,
}

impl Basis for CustomNested
{
    fn eval(&self, level: u32, index: u32, x: f64) -> f64 {
        let xl = self.nodes[level as usize][index  - 1];        
        let xr = self.nodes[level as usize][index  + 1];
        let xm = self.nodes[level as usize][index];
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
        self.nodes[level as usize][index]
    }

    fn eval_deriv(&self, _level: u32, _index: u32, _x: f64) -> f64 {
        panic!("Derivative not implemented for CurtisClenshaw Basis.")
    }

    fn degree(&self) -> usize {
        1
    }

    fn integral(&self, level: u32, index: u32) -> f64 {
        self.weights[level as usize][index]
    }

    fn basis_type(&self) -> super::base::BasisFunction {
        super::base::BasisFunction::CustomBasis
    }

    fn num_nodes(&self, level: u32) -> usize
    {
        self.nodes[level as usize].len()
    }
    
    fn node(&self, level: u32, index: u32) -> f64 {
        self.nodes[level as usize][index]
    }
    
    fn weights(&self, level: u32) -> Vec<f64> {
        self.weights[level as usize].0.clone()
    }
    
    fn nodes(&self, level: u32) -> Vec<f64> {
        self.nodes[level as usize].0.clone()
    }
    
    fn quadrature_exactness(&self, level: u32) -> u32
    {
        self.quadrature_exactness[level as usize]
    }

    fn interpolation_exactness(&self, level: u32) -> u32
    {
        self.interpolation_exactness[level as usize]
    }

}

