use indexmap::IndexSet;

use crate::storage::linear_grid::SparseGridData;

use super::refinement::RefinementFunctor;

pub(crate) fn coarsen<const D: usize, const DIM_OUT: usize>(storage: &mut SparseGridData<D>, functor: &dyn RefinementFunctor<D, DIM_OUT>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]], threshold: f64) -> IndexSet<usize>
{
    let mut kept_points = IndexSet::default();
    let zero_index = storage.adjacency_data.zero_index;
    let values = functor.eval(alpha, values);
    for (seq, (point, r)) in storage.nodes().iter().zip(values.into_iter()).enumerate()
    {
        if point.is_leaf() && seq != zero_index        
        {            
            if r >= threshold
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