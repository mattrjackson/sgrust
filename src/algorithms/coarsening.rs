use indexmap::IndexSet;

use crate::storage::linear_grid::SparseGridData;

use super::refinement::RefinementFunctor;

pub(crate) fn coarsen<const D: usize, const DIM_OUT: usize>(storage: &mut SparseGridData<D>, functor: &dyn RefinementFunctor<D, DIM_OUT>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]]) -> IndexSet<usize>
{
    let mut kept_points = IndexSet::default();
    let zero_index = storage.adjacency_data.zero_index;
    for (seq, point) in storage.nodes().iter().enumerate()
    {
        if point.is_leaf() && seq != zero_index        
        {
            let r = functor.eval(storage, alpha, values, seq);
            if r >= functor.threshold()
            {
                kept_points.insert(seq);
            }
        }
        else
        {
            kept_points.insert(seq);
        }
    }
    storage.remove(&kept_points);
    kept_points
}