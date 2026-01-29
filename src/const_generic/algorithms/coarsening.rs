use indexmap::IndexSet;

use crate::const_generic::storage::SparseGridData;

use super::refinement::RefinementFunctor;

/// Determines which boundary nodes are required by kept inner nodes.
/// A boundary node is required if any kept inner node references it via left_zero or right_zero.
fn find_required_boundary_nodes<const D: usize>(storage: &SparseGridData<D>, kept_points: &IndexSet<usize>) -> IndexSet<usize>
{
    let mut required_boundary = IndexSet::new();
    
    for &seq in kept_points.iter()
    {
        let point = &storage.nodes()[seq];
        
        // Only inner nodes have left_zero/right_zero dependencies
        if !point.is_inner_point()
        {
            continue;
        }
        
        for dim in 0..D
        {
            let offset = dim * storage.len();
            let left_zero_idx = storage.adjacency_data.left_zero[offset + seq];
            let right_zero_idx = storage.adjacency_data.right_zero[offset + seq];
            
            if left_zero_idx != u32::MAX
            {
                required_boundary.insert(left_zero_idx as usize);
            }
            if right_zero_idx != u32::MAX
            {
                required_boundary.insert(right_zero_idx as usize);
            }
        }
    }
    
    required_boundary
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
/// IndexSet containing the indices of points that were kept
pub(crate) fn coarsen<const D: usize, const DIM_OUT: usize>(storage: &mut SparseGridData<D>, functor: &dyn RefinementFunctor<D, DIM_OUT>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]], threshold: f64, remove_boundary: bool) -> IndexSet<usize>
{
    let mut kept_points = IndexSet::default();
    let zero_index = storage.adjacency_data.zero_index;
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
                    kept_points.insert(seq);
                }
            }
            else
            {
                // Keep non-leaf inner nodes and zero index
                kept_points.insert(seq);
            }
        }
    }
    
    if remove_boundary
    {
        // Find which boundary nodes are required by kept inner nodes
        let required_boundary = find_required_boundary_nodes(storage, &kept_points);
        
        // Add boundary nodes based on various criteria
        for (seq, (point, r)) in storage.nodes().iter().zip(values.iter()).enumerate()
        {
            if !point.is_inner_point()
            {
                // Always keep zero index
                if seq == zero_index
                {
                    kept_points.insert(seq);
                }
                // Keep if required by a kept inner node
                else if required_boundary.contains(&seq)
                {
                    kept_points.insert(seq);
                }
                // Keep if not a leaf (has children that are kept)
                else if !point.is_leaf()
                {
                    kept_points.insert(seq);
                }
                // Keep boundary leaf nodes with significant surplus
                else if *r >= threshold
                {
                    kept_points.insert(seq);
                }
                // Otherwise, this boundary node can be removed
            }
        }
    }
    else
    {
        // Keep all boundary nodes
        for (seq, point) in storage.nodes().iter().enumerate()
        {
            if !point.is_inner_point()
            {
                kept_points.insert(seq);
            }
        }
    }
    
    storage.remove(&kept_points);
    kept_points
}