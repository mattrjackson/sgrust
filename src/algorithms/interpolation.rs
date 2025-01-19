use num_traits::Float;

use crate::{ basis::base::Basis, errors::SGError, iterators::{grid_iterator::GridIteratorT, grid_iterator_cache::AdjacencyGridIterator}};
use super::{basis_evaluation::BasisEvaluation, basis_evalution_with_boundary::eval_boundary};

pub(crate) struct InterpolationOperation<'a, const D: usize, const DIM_OUT: usize, BASIS: Basis>(pub bool, pub BasisEvaluation<'a, D, DIM_OUT, BASIS>);


impl<const D: usize, const DIM_OUT: usize, BASIS: Basis> InterpolationOperation<'_, D, DIM_OUT, BASIS>
{

    #[inline]
    pub(crate) fn interpolate<T: Float  + std::ops::AddAssign, Iterator: GridIteratorT<D>>(&self, x: [f64; D], alpha: &[[T; DIM_OUT]], iterator: &mut Iterator) -> Result<[T; DIM_OUT], SGError>
    {
        match self.0
        {
            true =>
            {
                let mut result = [T::zero(); DIM_OUT];                            
                iterator.reset_to_level_zero();
                let xscaled = self.1.0.bounding_box.to_unit_coordinate(&x);
                eval_boundary(self.1.0, &self.1.1, &xscaled, 0, T::from(1.0).unwrap(), iterator, alpha, &mut result);    
                Ok(result)
                },
            false =>  self.1.eval(x, alpha, iterator),
        }       
    }
}

