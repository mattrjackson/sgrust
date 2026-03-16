use crate::dynamic::{iterators::dynamic_grid_iterator::DynamicHashMapGridIterator, storage::SparseGridData};

use super::refinement::RefinementFunctor;

/// Expands the kept-node set to include the level-zero boundary closure required by
/// boundary hierarchization sweeps.
///
/// For every kept node and every dimension, both left/right level-zero projections in that
/// dimension must exist on the active fiber. Applying this recursively yields the minimal
/// boundary skeleton, e.g. `2^D` corners for a constant boundary grid.
fn add_required_boundary_closure(iterator: &mut DynamicHashMapGridIterator, storage: &SparseGridData, keep: &mut [bool])
{
    use crate::dynamic::iterators::dynamic_grid_iterator::GridIteratorT;
    let mut changed = true;
    while changed
    {
        changed = false;
        for seq in 0..storage.len()
        {
            if !keep[seq]
            {
                continue;
            }
            let point = storage.point(seq);
            for dim in 0..storage.num_inputs()
            {
                iterator.set_index(point.clone());
                if iterator.reset_to_left_level_zero(dim)
                {
                    let boundary_seq = iterator.index().unwrap();
                    if !keep[boundary_seq]
                    {
                        keep[boundary_seq] = true;
                        changed = true;
                    }
                }
                iterator.set_index(point.clone());
                if iterator.reset_to_right_level_zero(dim)
                {
                    let boundary_seq = iterator.index().unwrap();
                    if !keep[boundary_seq]
                    {
                        keep[boundary_seq] = true;
                        changed = true;
                    }
                }
            }
        }
    }
}

fn keep_mask_to_indices(keep: &[bool]) -> Vec<usize>
{
    let mut kept_points = Vec::with_capacity(keep.iter().filter(|&&is_kept| is_kept).count());
    for (seq, &is_kept) in keep.iter().enumerate()
    {
        if is_kept
        {
            kept_points.push(seq);
        }
    }
    kept_points
}


/// Coarsen the grid by removing leaf nodes with low error estimates.
/// 
/// # Arguments
/// * `storage` - The sparse grid storage to coarsen
/// * `functor` - Refinement functor that computes error estimates
/// * `alpha` - Hierarchical surplus coefficients (flattened)
/// * `values` - Function values at grid points (flattened)
/// * `threshold` - Nodes with error below this are candidates for removal
/// * `remove_boundary` - If true, also remove boundary nodes that are no longer needed
///
/// # Returns
/// Vec containing the indices of points that were kept, in original storage order.
pub(crate) fn coarsen(storage: &mut SparseGridData, functor: &dyn RefinementFunctor, alpha: &[f64], values: &[f64], threshold: f64, remove_boundary: bool) -> Vec<usize>
{
    let mut iterator = DynamicHashMapGridIterator::new(storage);
    let mut keep = vec![false; storage.len()];
    let zero_index = storage.adjacency_data.zero_index;
    let values = functor.eval(storage.points(), alpha, values);
    
    // First pass: determine which inner nodes to keep
    for (seq, r) in values.iter().enumerate()
    {
        let is_leaf = storage.is_leaf(seq);
        let is_inner = storage.is_inner_point(seq);
        
        if is_inner
        {
            if is_leaf && seq != zero_index
            {
                if *r >= threshold
                {
                    keep[seq] = true;
                }
            }
            else
            {
                // Keep non-leaf inner nodes and zero index
                keep[seq] = true;
            }
        }
    }
    
    if remove_boundary
    {
        // Keep structurally significant boundary nodes before adding the required closure.
        for (seq, r) in values.iter().enumerate()
        {
            if !storage.is_inner_point(seq)
            {
                let is_leaf = storage.is_leaf(seq);
                
                // Always keep zero index
                if seq == zero_index
                {
                    keep[seq] = true;
                }
                // Keep if not a leaf (has children that are kept)
                else if !is_leaf
                {
                    keep[seq] = true;
                }
                // Keep boundary leaf nodes with significant surplus
                else if *r >= threshold
                {
                    keep[seq] = true;
                }
                // Otherwise, this boundary node can be removed
            }
        }
        add_required_boundary_closure(&mut iterator, storage, &mut keep);
    }
    else
    {
        // Keep all boundary nodes
        for seq in 0..storage.len()
        {
            if !storage.is_inner_point(seq)
            {
                keep[seq] = true;
            }
        }
    }
    
    let kept_points = keep_mask_to_indices(&keep);
    storage.remove(&kept_points);
    kept_points
}
