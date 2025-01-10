use crate::{ basis::base::Basis, errors::SGError, iterators::grid_iterator_cache::GridIteratorWithCache};
use super::{basis_evaluation::BasisEvaluation, basis_evalution_with_boundary::eval_boundary};

pub(crate) struct InterpolationOperation<'a, const D: usize, const DIM_OUT: usize, BASIS: Basis>(pub bool, pub BasisEvaluation<'a, D, DIM_OUT, BASIS>);


impl<const D: usize, const DIM_OUT: usize, BASIS: Basis> InterpolationOperation<'_, D, DIM_OUT, BASIS>
{

    #[inline]
    pub(crate) fn interpolate(&self, x: [f64; D], alpha: &[[f64; DIM_OUT]], iterator: &mut GridIteratorWithCache<D>) -> Result<[f64; DIM_OUT], SGError>
    {
        match self.0
        {
            true =>
            {
                let mut result = [0.0; DIM_OUT];                            
                iterator.reset_to_level_zero();
                let xscaled = if let Some(bbox) = self.1.0.bounding_box()
                {
                    bbox.to_unit_coordinate(&x)
                }
                else
                {
                    x
                };
                eval_boundary(self.1.0, &self.1.1, &xscaled, 0, 1.0, iterator, alpha, &mut result);    
                Ok(result)
                },
            false =>  self.1.eval(x, alpha),
        }       
    }
}

