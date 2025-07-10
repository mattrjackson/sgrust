use num_traits::Float;

use crate::{basis::{base::Basis, linear::LinearBasis}, dynamic::{iterators::dynamic_grid_iterator::GridIteratorT, storage::SparseGridData}};

#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) fn eval_boundary<T: Float +std::ops::AddAssign, Iterator: GridIteratorT>(storage: &SparseGridData, basis: &[LinearBasis], x: &[f64], 
    dim: usize, value: T, iterator: &mut Iterator, alpha: &[T], result: &mut [T], ndim: usize, num_outputs: usize)
{
    let mut level = 0;
    loop
    {
        let node_index = iterator.index().unwrap();
        let work_index = iterator.point_index(dim);        
        if level > 0
        {
            let new_value = T::from(basis[dim].eval(level, work_index, x[dim])).unwrap();
            if dim == ndim - 1
            {
                #[allow(clippy::needless_range_loop)]
                for i in 0..num_outputs
                {
                    result[i] += alpha[node_index*num_outputs+i] * value * new_value;
                }
            }
            else 
            {
                eval_boundary(storage, basis, x, dim + 1, value * new_value, iterator, alpha, result, ndim, num_outputs);    
            }
        }
        else
        {
            // handle boundaries if we are on level 0
            // level 0, index 0
            // reset_to_left_level_zero now checks if the node exists - after grid coarsening some boundary nodes are removed.            
            if iterator.reset_to_left_level_zero(dim)
            {
                let seq_l = iterator.index().unwrap();
                let new_value_l = T::from(basis[dim].eval(0, 0, x[dim])).unwrap();
                if dim == ndim - 1
                {
                    #[allow(clippy::needless_range_loop)]
                    for i in 0..num_outputs
                    {
                        result[i] += alpha[seq_l*num_outputs+i] * value * new_value_l;
                    }
                }
                else 
                {
                    eval_boundary(storage, basis, x, dim + 1, value * new_value_l, iterator, alpha, result, ndim, num_outputs);
                }
            }
            // reset_to_right_level_zero now checks if the node exists - after grid coarsening some boundary nodes are removed.
            if iterator.reset_to_right_level_zero(dim)
            {
                let seq_r = iterator.index().unwrap();
                let new_value_r = T::from(basis[dim].eval(0, 1, x[dim])).unwrap();
                if dim == ndim - 1
                {
                    #[allow(clippy::needless_range_loop)]
                    for i in 0..num_outputs
                    {
                        result[i] += alpha[seq_r*num_outputs+i] * value * new_value_r;
                    }
                }
                else 
                {
                    eval_boundary(storage, basis, x, dim + 1, value * new_value_r, iterator, alpha, result, ndim, num_outputs);
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
