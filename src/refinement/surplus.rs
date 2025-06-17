use crate::{algorithms::refinement::RefinementFunctor, storage::linear_grid::PointIterator};

pub struct SurplusRefinement<const D: usize, const DIM_OUT: usize>;

impl<const D: usize, const DIM_OUT: usize> RefinementFunctor<D, DIM_OUT> for SurplusRefinement<D, DIM_OUT>
{
    fn eval(&self, _points: PointIterator<D>, alpha: &[[f64; DIM_OUT]], _values: &[[f64; DIM_OUT]]) -> Vec<f64>
    {
        alpha.iter().map(|alpha_i|
        {
            let mut max = -1.0_f64;
            alpha_i.iter().for_each(|&val| max = max.max(val.abs()));
            max
        }).collect()        
    }
}