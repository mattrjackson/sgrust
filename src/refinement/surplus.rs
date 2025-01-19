use crate::algorithms::refinement::RefinementFunctor;

pub struct SurplusRefinement<const D: usize, const DIM_OUT: usize>(pub f64);

impl<const D: usize, const DIM_OUT: usize> RefinementFunctor<D, DIM_OUT> for SurplusRefinement<D, DIM_OUT>
{
    fn eval(&self, _storage: &crate::storage::linear_grid::SparseGridData<D>, alpha: &[[f64; DIM_OUT]], _values: &[[f64; DIM_OUT]], seq: usize) -> f64 
    {
        let mut max = -1.0_f64;
        alpha[seq].iter().for_each(|&val| max = max.max(val.abs()));
        max
    }

    fn threshold(&self) -> f64 {
        self.0
    }
}