use serde::{Deserialize, Serialize};

use super::{base::Basis, clenshaw_curtis::ClenshawCurtis, custom_nested::CustomNested, gauss_patterson::GaussPatterson, linear::LinearBasis};

#[derive(Serialize, Deserialize)]
#[derive(Hash, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum GlobalBasisType
{
    Linear = 0,
    ClenshawCurtis = 1,
    GaussPatterson = 2,
    CustomBasis = 3,
}

#[derive(Serialize, Deserialize)]
#[derive(Clone)]
pub struct GlobalBasis
{
    pub basis_type: GlobalBasisType,
    pub custom_rule: Option<CustomNested>,
}
impl GlobalBasis
{
    pub fn new(basis_type: GlobalBasisType) -> Self
    {
        Self { basis_type, custom_rule: None }
    }
}
impl From<GlobalBasisType> for GlobalBasis
{
    fn from(value: GlobalBasisType) -> Self {
        GlobalBasis { basis_type: value, custom_rule: None }
    }
}
impl Basis for GlobalBasis
{
    fn eval(&self, level: u32, index: u32, x: f64) -> f64 {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.eval(level, index, x),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.eval(level, index, x),
            GlobalBasisType::GaussPatterson => GaussPatterson.eval(level, index, x),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.eval(level, index, x)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn eval_deriv(&self, level: u32, index: u32, x: f64) -> f64 {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.eval_deriv(level, index, x),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.eval_deriv(level, index, x),
            GlobalBasisType::GaussPatterson => GaussPatterson.eval_deriv(level, index, x),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.eval_deriv(level, index, x)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn degree(&self) -> usize {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.degree(),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.degree(),
            GlobalBasisType::GaussPatterson => GaussPatterson.degree(),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.degree()
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn integral(&self, level: u32, index: u32) -> f64 {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.integral(level, index),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.integral(level, index),
            GlobalBasisType::GaussPatterson => GaussPatterson.integral(level, index),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.integral(level, index)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn basis_type(&self) -> crate::basis::base::BasisFunction {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.basis_type(),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.basis_type(),
            GlobalBasisType::GaussPatterson => GaussPatterson.basis_type(),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.basis_type()
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn node(&self, level: u32, index: u32) -> f64 {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.node(level, index),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.node(level, index),
            GlobalBasisType::GaussPatterson => GaussPatterson.node(level, index),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.node(level, index)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn num_nodes(&self, level: u32) -> usize {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.num_nodes(level),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.num_nodes(level),
            GlobalBasisType::GaussPatterson => GaussPatterson.num_nodes(level),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.num_nodes(level)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn get_point(&self, level: u32, index: u32) -> f64 {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.get_point(level, index),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.get_point(level, index),
            GlobalBasisType::GaussPatterson => GaussPatterson.get_point(level, index),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.get_point(level, index)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn weights(&self, level: u32) -> Vec<f64> {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.weights(level),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.weights(level),
            GlobalBasisType::GaussPatterson => GaussPatterson.weights(level),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.weights(level)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn nodes(&self, level: u32) -> Vec<f64> {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.nodes(level),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.nodes(level),
            GlobalBasisType::GaussPatterson => GaussPatterson.nodes(level),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.nodes(level)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }
    
    fn quadrature_exactness(&self, level: u32) -> u32 {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.quadrature_exactness(level),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.quadrature_exactness(level),
            GlobalBasisType::GaussPatterson => GaussPatterson.quadrature_exactness(level),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.quadrature_exactness(level)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }

    fn interpolation_exactness(&self, level: u32) -> u32 {
        match self.basis_type
        {
            GlobalBasisType::Linear => LinearBasis.interpolation_exactness(level),
            GlobalBasisType::ClenshawCurtis => ClenshawCurtis.interpolation_exactness(level),
            GlobalBasisType::GaussPatterson => GaussPatterson.interpolation_exactness(level),
            GlobalBasisType::CustomBasis => if let Some(rule) = &self.custom_rule
            {
                rule.interpolation_exactness(level)
            }
            else
            {
                panic!("CustomBasis selected but no rule defined.")
            },
        }
    }
}