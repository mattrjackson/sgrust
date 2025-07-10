use crate::dynamic::{algorithms::refinement::RefinementFunctor, storage::PointIterator};

pub struct SurplusRefinement(pub usize, pub usize);

impl RefinementFunctor for SurplusRefinement
{
    fn eval(&self, _points: PointIterator, alpha: &[f64], _values: &[f64]) -> Vec<f64>
    {
        alpha.chunks_exact(self.num_outputs()).map(|alpha_i|
        {
            let mut max = -1.0_f64;
            alpha_i.iter().for_each(|&val| max = max.max(val.abs()));
            max
        }).collect()        
    }
    #[inline]
    fn num_outputs(&self) -> usize {
        self.1
    }

    #[inline]    
    fn num_inputs(&self) -> usize {
        self.0
    }
}