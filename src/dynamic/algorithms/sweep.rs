use crate::{dynamic::{iterators::dynamic_grid_iterator::{DynamicHashMapGridIterator, GridIteratorT}, storage::SparseGridData}, errors::SGError};

pub trait SweepFunction<T, R>
{
    fn execute_in_place(&mut self, values: &mut [R], iterator: &mut DynamicHashMapGridIterator, storage: &SparseGridData, dimension: usize) -> Result<(), SGError>;
}

#[allow(unused)]
fn sweep_recursive_in_place<T, R, F: SweepFunction<T, R>>(function: &mut F, storage: &SparseGridData, values: &mut [R], iterator: &mut DynamicHashMapGridIterator, 
    dim_list: &Vec<usize>, dim_rem: usize, dim_sweep: usize)
    {
        function.execute_in_place(values, iterator, storage, dim_sweep);
        for d in 0..dim_rem
        {
            let cur_dim = dim_list[d];
            if iterator.is_leaf()
            {
                continue;
            }
            iterator.left_child(cur_dim);
            if iterator.seq().is_some()
            {
                sweep_recursive_in_place(function, storage, values, iterator, dim_list, d + 1, dim_sweep);
            }
            iterator.step_right(cur_dim);
            if iterator.seq().is_some()
            {
                sweep_recursive_in_place(function, storage, values, iterator, dim_list, d + 1, dim_sweep);
            }
            iterator.up(cur_dim);
        }
    }


pub(crate) fn sweep_boundary_recursive_in_place<T, R, F: SweepFunction<T, R>>(function: &mut F, storage: &SparseGridData, values: &mut [R], iterator: &mut DynamicHashMapGridIterator, 
    dim_list: &Vec<usize>, dim_rem: usize, dim_sweep: usize) -> Result<(), SGError>
    {
        if dim_rem == 0
        {
            function.execute_in_place(values, iterator, storage, dim_sweep)?;
        }
        else
        {
            let current_index = iterator.point();
            let d = dim_list[dim_rem - 1];
            let current_level = current_index.level[d];
            if current_level > 0
            {
                sweep_boundary_recursive_in_place(function, storage, values, iterator, dim_list, dim_rem - 1, dim_sweep)?;
                if !iterator.is_leaf()
                {
                    iterator.left_child(d);
                    if iterator.seq().is_some()
                    {
                        sweep_boundary_recursive_in_place(function, storage, values, iterator, dim_list, dim_rem, dim_sweep)?;
                    }
                    iterator.step_right(d);
                    if iterator.seq().is_some()
                    {
                        sweep_boundary_recursive_in_place(function, storage, values, iterator, dim_list, dim_rem, dim_sweep)?;
                    }
                    iterator.up(d);
                }
            }
            else 
            {
                sweep_boundary_recursive_in_place(function, storage, values, iterator, dim_list, dim_rem - 1, dim_sweep)?;
                iterator.reset_to_right_level_zero(d);
                sweep_boundary_recursive_in_place(function, storage, values, iterator, dim_list, dim_rem - 1, dim_sweep)?;
                if !iterator.is_leaf()
                {
                    iterator.reset_to_level_one(d);
                    if iterator.seq().is_some()
                    {
                        sweep_boundary_recursive_in_place(function, storage, values, iterator, dim_list, dim_rem, dim_sweep)?;
                    }
                }
                iterator.reset_to_left_level_zero(d);
            }                   
        }
        Ok(())
    }
pub(crate) fn sweep_1d_in_place<T, R, F: SweepFunction<T, R>>(function: &mut F, storage: &SparseGridData, values: &mut [R],
    dim_sweep: usize, ndim: usize)
{
    let mut dim_list = vec![0; ndim-1];
    let mut idx = 0;
    for i in 0..ndim
    {
        if i != dim_sweep
        {
            dim_list[idx] = i;
            idx += 1;
        }
    }
    let mut iterator = DynamicHashMapGridIterator::new(storage);
    sweep_recursive_in_place(function, storage, values, &mut iterator, &dim_list, ndim-1, dim_sweep);
}

pub(crate) fn sweep_1d_boundary_in_place<T, R, F: SweepFunction<T, R>>(function: &mut F, storage: &SparseGridData, values: &mut [R],
    dim_sweep: usize, ndim: usize) -> Result<(), SGError>
{
    let mut dim_list = vec![0; ndim-1];
    let mut idx = 0;
    for i in 0..ndim
    {
        if i != dim_sweep
        {
            dim_list[idx] = i;
            idx += 1;
        }
    }
    let mut iterator = DynamicHashMapGridIterator::new(storage);
    iterator.reset_to_level_zero();
    sweep_boundary_recursive_in_place(function, storage, values, &mut iterator, &dim_list, ndim-1, dim_sweep)
}