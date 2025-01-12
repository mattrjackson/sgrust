use num_traits::Float;

use crate::{basis::base::Basis, iterators::{grid_iterator::GridIteratorT, grid_iterator_cache::GridIteratorWithCache}, storage::linear_grid::SparseGridData};

#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) fn eval_boundary<const D: usize, const DIM_OUT: usize, BASIS: Basis, T: Float +std::ops::AddAssign>(storage: &SparseGridData<D>, basis: &[BASIS; D], x: &[f64; D], 
    dim: usize, value: T, iterator: &mut GridIteratorWithCache<D>, alpha: &[[T; DIM_OUT]], result: &mut [T; DIM_OUT])
{
    if value.abs().to_f64().unwrap() < 1e-14
    {
        return;
    }
    let mut level = 0;
    loop
    {
        let work_index = storage[iterator.index].index[dim];        
        if level > 0
        {
            let new_value = T::from(basis[dim].eval(level, work_index, x[dim])).unwrap();
            if dim == D - 1
            {
                #[allow(clippy::needless_range_loop)]
                for i in 0..DIM_OUT
                {
                    result[i] += alpha[iterator.index][i] * value * new_value;
                }
            }
            else 
            {
                eval_boundary(storage, basis, x, dim + 1, value * new_value, iterator, alpha, result);    
            }
        }
        else
        {
            // handle boundaries if we are on level 0
            // level 0, index 0
            // reset_to_left_level_zero now checks if the node exists - after grid coarsening some boundary nodes are removed.            
            if iterator.reset_to_left_level_zero(dim)
            {
                let seq_l = iterator.index;
                let new_value_l = T::from(basis[dim].eval(0, 0, x[dim])).unwrap();
                if dim == D - 1
                {
                    #[allow(clippy::needless_range_loop)]
                    for i in 0..DIM_OUT
                    {
                        result[i] += alpha[seq_l][i] * value * new_value_l;
                    }
                }
                else 
                {
                    eval_boundary(storage, basis, x, dim + 1, value * new_value_l, iterator, alpha, result);
                }
            }
            // reset_to_right_level_zero now checks if the node exists - after grid coarsening some boundary nodes are removed.
            if iterator.reset_to_right_level_zero(dim)
            {
                let seq_r = iterator.index;
                let new_value_r = T::from(basis[dim].eval(0, 1, x[dim])).unwrap();
                if dim == D - 1
                {
                    #[allow(clippy::needless_range_loop)]
                    for i in 0..DIM_OUT
                    {
                        result[i] += alpha[seq_r][i] * value * new_value_r;
                    }
                }
                else 
                {
                    eval_boundary(storage, basis, x, dim + 1, value * new_value_r, iterator, alpha, result);
                }
            }
        }
        // Can't descend any further.
        if iterator.is_leaf()
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
            if x > x_coord && !iterator.right_child(dim)
            {                
                break;
            }             
            if x <= x_coord && !iterator.left_child(dim)
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
            if !iterator.reset_to_level_one(dim)
            {
                break;
            }            
        }
        level += 1;
    }
    iterator.reset_to_left_level_zero(dim);    
}
