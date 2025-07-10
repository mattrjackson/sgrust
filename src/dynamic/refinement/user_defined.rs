use crate::dynamic::{algorithms::refinement::RefinementFunctor, storage::PointIterator};

///
/// A function that defines how refinement is performed.
/// 
/// # Arguments
/// - `storage`: Storage of sparse grid.
/// - `alpha`: Surplus Coefficients
/// 
pub type UserRefinementFunction = dyn Fn(&[f64], &[f64]) -> f64 + Send + Sync;
pub struct UserDefinedRefinement<'a>
{
    pub fun_eval: &'a UserRefinementFunction,
    pub num_inputs: usize,
    pub num_outputs: usize,
}

impl RefinementFunctor for UserDefinedRefinement<'_>
{
    fn eval(&self, _points: PointIterator, alpha: &[f64], values: &[f64]) -> Vec<f64> {
        
        alpha.chunks_exact(self.num_outputs()).zip(values.chunks_exact(self.num_outputs())).map(|(alpha_i, values_i)|
        {
            (self.fun_eval)(alpha_i, values_i)
        }).collect()    
    }
    
    fn num_outputs(&self) -> usize {
        self.num_outputs
    }
    
    fn num_inputs(&self) -> usize {
        self.num_inputs
    }
}