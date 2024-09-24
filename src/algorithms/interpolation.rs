use crate::{ basis::base::Basis, errors::SGError, iterators::grid_iterator_cache::GridIteratorWithCache};

use super::basis_evaluation::BasisEvaluation;

pub(crate) struct InterpolationOperation<'a, const D: usize, const DIM_OUT: usize, BASIS: Basis>(pub bool, pub BasisEvaluation<'a, D, DIM_OUT, BASIS>);


impl<'a, const D: usize, const DIM_OUT: usize, BASIS: Basis> InterpolationOperation<'a, D, DIM_OUT, BASIS>
{
    #[inline]
    fn process_affected_basis_functions(&self, alpha: &[[f64; DIM_OUT]], results: Vec<(usize, f64)>) -> Result<[f64; DIM_OUT], SGError>
    {
        let mut r = [0.0; DIM_OUT];
        for (seq, value) in results
        {
            (0..DIM_OUT).for_each(|d| {
                r[d] += alpha[seq][d] * value;
            });
        }
        Ok(r)
    }
    #[inline]
    pub(crate) fn interpolate(&self, x: [f64; D], alpha: &[[f64; DIM_OUT]], iterator: &mut GridIteratorWithCache<D>) -> Result<[f64; DIM_OUT], SGError>
    {
        match self.0
        {
            true => self.process_affected_basis_functions(alpha, 
                crate::algorithms::affected_basis_functions::get_affected_basis_functions(x, &self.1.1, self.1.0, iterator)),
            false =>  self.1.eval(x, alpha),
        }       
    }
}

