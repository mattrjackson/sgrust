use crate::algorithms::refinement::RefinementFunctor;

///
/// A function that defines how refinement is performed.
/// 
/// # Arguments
/// - `storage`: Storage of sparse grid.
/// - `alpha`: Surplus Coefficients
/// 
pub type UserRefinementFunction<const D: usize, const DIM_OUT: usize> = dyn Fn(&crate::storage::linear_grid::SparseGridStorage<D>,  &[[f64; DIM_OUT]], &[[f64; DIM_OUT]], usize) -> f64 + Send + Sync;
pub struct UserDefinedRefinement<'a, const D: usize, const DIM_OUT: usize>
{
    pub fun_eval: &'a UserRefinementFunction<D, DIM_OUT>,
    pub threshold: f64,
}

impl<const D: usize, const DIM_OUT: usize> RefinementFunctor<D, DIM_OUT> for UserDefinedRefinement<'_, D, DIM_OUT>
{
    fn eval(&self, storage: &crate::storage::linear_grid::SparseGridStorage<D>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]], seq: usize) -> f64 {
        (self.fun_eval)(storage, alpha, values, seq)
    }

    fn threshold(&self) -> f64 {
        self.threshold
    }
}