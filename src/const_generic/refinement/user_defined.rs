use crate::const_generic::{algorithms::refinement::RefinementFunctor, storage::PointIterator};

///
/// A function that defines how refinement is performed.
/// 
/// # Arguments
/// - `storage`: Storage of sparse grid.
/// - `alpha`: Surplus Coefficients
/// 
pub type UserRefinementFunction<const D: usize, const DIM_OUT: usize> = dyn Fn(&[f64; D], &[f64; DIM_OUT], &[f64; DIM_OUT]) -> f64 + Send + Sync;
pub struct UserDefinedRefinement<'a, const D: usize, const DIM_OUT: usize>
{
    pub fun_eval: &'a UserRefinementFunction<D, DIM_OUT>,
}

impl<const D: usize, const DIM_OUT: usize> RefinementFunctor<D, DIM_OUT> for UserDefinedRefinement<'_, D, DIM_OUT>
{
    fn eval(&self, points: PointIterator<D>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]]) -> Vec<f64> {
        points
            .zip(alpha.iter().zip(values.iter()))
            .map(|(point, (alpha_i, values_i))| {
                (self.fun_eval)(&point, alpha_i, values_i)
            })
            .collect()
    }
}