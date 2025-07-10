use indexmap::IndexSet;

use crate::dynamic::storage::SparseGridData;

use super::refinement::RefinementFunctor;

pub(crate) fn coarsen(storage: &mut SparseGridData, functor: &dyn RefinementFunctor, alpha: &[f64], values: &[f64], threshold: f64) -> IndexSet<usize>
{
    let mut kept_points = IndexSet::default();
    let zero_index = storage.adjacency_data.zero_index;
    let values = functor.eval(storage.points(), alpha, values);
    for (seq, (point, r)) in storage.nodes().zip(values.into_iter()).enumerate()
    {
        if point.flags.is_leaf() && seq != zero_index        
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