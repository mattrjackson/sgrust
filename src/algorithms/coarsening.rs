use std::collections::HashSet;

use crate::storage::linear_grid::SparseGridStorage;

use super::refinement::RefinementFunctor;

pub(crate) fn coarsen<const D: usize, const DIM_OUT: usize>(storage: &mut SparseGridStorage<D>, functor: &dyn RefinementFunctor<D, DIM_OUT>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]]) -> HashSet<usize>
{
    let mut removed_points = HashSet::new();
    for (seq, point) in storage.iter().enumerate()
    {
        if point.is_leaf()
        {
            let r = functor.eval(storage, alpha, values, seq);
            if r < functor.threshold()
            {
                removed_points.insert(seq);
            }
        }
    }
    storage.remove(&removed_points);
    removed_points
}