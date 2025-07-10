use num_traits::Float;

use crate::dynamic::iterators::dynamic_grid_iterator::GridIteratorT;
use super::{basis_evaluation::BasisEvaluation, basis_evalution_with_boundary::eval_boundary};
use crate::errors::SGError;
pub(crate) struct InterpolationOperation<'a>(pub bool, pub BasisEvaluation<'a>);


impl InterpolationOperation<'_>
{

    #[inline]
    pub(crate) fn interpolate<T: Float  + std::ops::AddAssign, Iterator: GridIteratorT>(&self, x: &[f64], alpha: &[T], iterator: &mut Iterator, result: &mut [T]) -> Result<(), SGError>
    {
        match self.0
        {
            true =>
            {
                iterator.reset_to_level_zero();
                let xscaled = self.1.0.bounding_box.to_unit_coordinate(&x);
                eval_boundary(self.1.0, &self.1.1, &xscaled, 0, T::from(1.0).unwrap(), iterator, alpha, result, x.len(), result.len());    
                Ok(())
                },
            false =>  self.1.eval(x, alpha, iterator, result),
        }       
    }
}

