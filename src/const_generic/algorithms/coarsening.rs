use crate::const_generic::{iterators::grid_iterator::{GridIteratorT, HashMapGridIterator}, storage::SparseGridData};

use super::refinement::RefinementFunctor;

/// Expands the kept-node set to include the level-zero boundary closure required by
/// boundary hierarchization sweeps.
///
/// For every kept node and every dimension, both left/right level-zero projections in that
/// dimension must exist on the active fiber. Applying this recursively yields the minimal
/// boundary skeleton, e.g. `2^D` corners for a constant boundary grid.
fn add_required_boundary_closure<const D: usize>(iterator: &mut HashMapGridIterator<D>, storage: &SparseGridData<D>, keep: &mut [bool])
{
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
            let point = storage.nodes()[seq];
            for dim in 0..D
            {
                iterator.set_index(point);
                if iterator.reset_to_left_level_zero(dim)
                {
                    let boundary_seq = iterator.index().unwrap();
                    if !keep[boundary_seq]
                    {
                        keep[boundary_seq] = true;
                        changed = true;
                    }
                }
                iterator.set_index(point);
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
    let mut kept_points = Vec::with_capacity(keep.iter().filter(|&&keep| keep).count());
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
/// * `alpha` - Hierarchical surplus coefficients
/// * `values` - Function values at grid points
/// * `threshold` - Nodes with error below this are candidates for removal
/// * `remove_boundary` - If true, also remove boundary nodes that are no longer needed
///
/// # Returns
/// Vec containing the indices of points that were kept, in original storage order.
pub(crate) fn coarsen<const D: usize, const DIM_OUT: usize>( storage: &mut SparseGridData<D>, functor: &dyn RefinementFunctor<D, DIM_OUT>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]], threshold: f64, remove_boundary: bool) -> Vec<usize>
{
    let mut iterator = HashMapGridIterator::new(storage);
    let mut keep = vec![false; storage.len()];
    iterator.reset_to_level_zero();
    let zero_index = iterator.index().unwrap();
    let values = functor.eval(storage.points(), alpha, values);
    
    // First pass: determine which inner nodes to keep
    for (seq, (point, r)) in storage.nodes().iter().zip(values.iter()).enumerate()
    {
        if point.is_inner_point()
        {
            if point.is_leaf() && seq != zero_index
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
        for (seq, (point, r)) in storage.nodes().iter().zip(values.iter()).enumerate()
        {
            if !point.is_inner_point()
            {
                // Always keep zero index
                if seq == zero_index
                {
                    keep[seq] = true;
                }
                // Keep if not a leaf (has children that are kept)
                else if !point.is_leaf()
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
        for (seq, point) in storage.nodes().iter().enumerate()
        {
            if !point.is_inner_point()
            {
                keep[seq] = true;
            }
        }
    }
    let kept_points = keep_mask_to_indices(&keep);
    storage.remove(&kept_points);
    kept_points
}