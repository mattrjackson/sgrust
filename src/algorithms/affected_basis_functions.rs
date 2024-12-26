use crate::{basis::base::Basis, iterators::grid_iterator_cache::GridIteratorWithCache, storage::linear_grid::SparseGridStorage};

#[inline]
fn rec_linear<const D: usize, BASIS: Basis>(storage: &SparseGridStorage<D>, basis: &[BASIS; D], x: &[f64; D], 
    dim: usize, value: f64, iterator: &mut GridIteratorWithCache<D>, result: &mut Vec<(usize, f64)>)
{
    if value.abs() < 1e-14
    {
        return;
    }
    let mut level = 0;
    loop
    {
        let work_index = storage[iterator.index].index[dim];        
        if level > 0
        {
            let new_value = basis[dim].eval(level, work_index, x[dim]);
            if dim == D - 1
            {
                result.push((iterator.index, value * new_value));
            }
            else 
            {
                rec_linear(storage, basis, x, dim + 1, value * new_value, iterator, result);    
            }
        }
        else
        {
            // handle boundaries if we are on level 0
            // level 0, index 0
            // reset_to_left_level_zero now checks if the node exists - after grid coarsening some boundary nodes are removed.            
            if iterator.reset_to_left_level_zero(dim, storage)
            {
                let seq_l = iterator.index;
                let new_value_l = basis[dim].eval(0, 0, x[dim]);
                if dim == D - 1
                {
                    result.push((seq_l, value * new_value_l));
                }
                else 
                {
                    rec_linear(storage, basis, x, dim + 1, value * new_value_l, iterator, result);
                }
            }
            // reset_to_right_level_zero now checks if the node exists - after grid coarsening some boundary nodes are removed.
            if iterator.reset_to_right_level_zero(dim, storage)
            {
                let seq_r = iterator.index;
                let new_value_r = basis[dim].eval(0, 1, x[dim]);
                if dim == D - 1
                {
                    result.push((seq_r, value * new_value_r));
                }
                else 
                {
                    rec_linear(storage, basis, x, dim + 1, value * new_value_r, iterator, result);
                }
            }
        }
        // Can't descend any further.
        if iterator.is_leaf(storage)
        {
            break;
        }
        // this decides in which direction we should descend by evaluating
        // the corresponding bit
        // the bits are coded from left to right starting with level 1
        // being in position max_level
        let x = x[dim];
        if level > 0
        {
            let h = 1 << level; 
            let x_coord = work_index as f64 / h as f64;            
            if (x - x_coord).abs() < 1e-15
            {
                break;
            }
            if x > x_coord && !iterator.right_child(dim, storage)
            {                
                break;
            }             
            if x <= x_coord && !iterator.left_child(dim, storage)
            {                 
                break;
            }
        }
        else 
        {
            if x.abs() < 1e-15 || (x-1.0).abs() < 1e-15
            {
                break;
            }            
            if !iterator.reset_to_level_one(dim, storage)
            {
                break;
            }            
        }
        level += 1;
    }
    iterator.reset_to_left_level_zero(dim, storage);    
}
#[inline]
pub(crate) fn get_affected_basis_functions<const D: usize, BASIS: Basis>(x: [f64; D], basis: &[BASIS; D], storage: &SparseGridStorage<D>,  iterator: &mut GridIteratorWithCache<D>) -> Vec<(usize, f64)>
{
    iterator.reset_to_level_zero();
    let mut result = Vec::with_capacity(2000);
    let xscaled = if let Some(bbox) = storage.bounding_box()
    {
        bbox.to_unit_coordinate(&x)
    }
    else
    {
        x
    };
    rec_linear(storage, basis, &xscaled, 0, 1.0, iterator, &mut result);    
    result
}