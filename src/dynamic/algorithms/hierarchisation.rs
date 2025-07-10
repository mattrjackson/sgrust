use crate::dynamic::{algorithms::sweep::SweepFunction, iterators::dynamic_grid_iterator::{DynamicHashMapGridIterator, GridIteratorT}, storage::SparseGridData};


pub trait HierarchisationOperation : Copy
{
    fn hierarchize(&self, node_values: &mut [f64], storage: &SparseGridData); 
    fn dehierarchize(&self, alpha: &mut [f64], storage: &SparseGridData); 
}

pub struct LinearHierarchisation;
impl LinearHierarchisation
{
    fn recurse(values: &mut [f64], iterator: &mut DynamicHashMapGridIterator, dimension: usize, left_value: &[f64], right_value: &[f64], ndim: usize)
    {
        if let Some(seq) = iterator.seq()
        {
            let mid_value = Vec::from_iter(values.chunks_exact(left_value.len()).nth(seq).unwrap().iter().copied());
            if !iterator.is_leaf()
            {            
                iterator.left_child(dimension);
                if let Some(_seq) = iterator.seq()
                {                    
                    Self::recurse(values, iterator, dimension, left_value, &mid_value, ndim);                    
                }
                iterator.step_right(dimension);
                if let Some(_seq) = iterator.seq()                
                {
                    Self::recurse(values, iterator, dimension, &mid_value, right_value, ndim);
                }
                iterator.up(dimension);
            }
            for i in 0..left_value.len()
            {
                values[seq*right_value.len()+i] = mid_value[i] - 0.5 * (left_value[i] + right_value[i]);
            }
        }
    }
}

impl SweepFunction<&[f64], f64> for LinearHierarchisation
{
    fn execute_in_place(&mut self, values: &mut [f64], iterator: &mut DynamicHashMapGridIterator, storage: &SparseGridData, dimension: usize) {
        let left = vec![0.0; storage.num_outputs];
        let right = vec![0.0; storage.num_outputs];
        Self::recurse(values, iterator, dimension, &left, &right, storage.num_inputs);
    }
}


pub struct LinearBoundaryHierarchisation(LinearHierarchisation);

impl SweepFunction<&[f64], f64> for LinearBoundaryHierarchisation
{
    #[inline]
    fn execute_in_place(&mut self, values: &mut [f64], iterator: &mut DynamicHashMapGridIterator, storage: &SparseGridData, dimension: usize) {
        // left boundary
        iterator.reset_to_left_level_zero(dimension);
        let seq_left = iterator.seq().unwrap();
        let left_boundary: Vec<f64> = values.chunks_exact(storage.num_outputs).nth(seq_left).unwrap().iter().copied().collect();
        iterator.reset_to_right_level_zero(dimension);
        let seq_right = iterator.seq().unwrap();       
        let right_boundary: Vec<f64> = values.chunks_exact(storage.num_outputs).nth(seq_right).unwrap().iter().copied().collect();
        if !iterator.is_leaf()
        {
            iterator.reset_to_level_one(dimension);
            if iterator.seq().is_some()
            {
                LinearHierarchisation::recurse(values, iterator, dimension, &left_boundary, &right_boundary, storage.num_inputs);
            }
            iterator.reset_to_left_level_zero(dimension);
        }
    }
}

pub struct LinearDehierarchisation;

impl LinearDehierarchisation
{
    #[inline]
    fn recurse(values: &mut [f64], iterator: &mut DynamicHashMapGridIterator, dimension: usize, left_value: &[f64], right_value: &[f64])
    {
        if let Some(seq) = iterator.seq()
        {
            let mut mid_value: Vec<f64> = values.chunks_exact(left_value.len()).nth(seq).unwrap().iter().copied().collect();
            for i in 0..left_value.len()
            {
                mid_value[i] += (left_value[i] + right_value[i]) * 0.5;
            }
            values.chunks_exact_mut(left_value.len()).nth(seq).unwrap().copy_from_slice(&mid_value);
            if !iterator.is_leaf()
            {
                iterator.left_child(dimension);
                if iterator.seq().is_some()
                {
                    Self::recurse(values, iterator, dimension, left_value, &mid_value);
                }
                iterator.step_right(dimension);
                if iterator.seq().is_some()
                {
                    Self::recurse(values, iterator, dimension, &mid_value, right_value);
                }
                iterator.up(dimension);
            }
        }
    }
}
impl SweepFunction<&[f64], f64> for LinearDehierarchisation
{
    #[inline]
    fn execute_in_place(&mut self, values: &mut [f64], iterator: &mut DynamicHashMapGridIterator, storage: &SparseGridData, dimension: usize) {
        let dummy = vec![0.0; storage.num_outputs];
        Self::recurse(values, iterator, dimension, &dummy, &dummy);
    }
}

pub struct LinearBoundaryDehierarchisation(LinearDehierarchisation);

impl SweepFunction<&[f64], f64> for LinearBoundaryDehierarchisation
{
    #[inline]
    fn execute_in_place(&mut self, values: &mut [f64], iterator: &mut DynamicHashMapGridIterator, storage: &SparseGridData, dimension: usize) {
        // left boundary
        iterator.reset_to_left_level_zero(dimension);
        let seq_left = iterator.seq().unwrap();
        let left_boundary: Vec<f64> = values.chunks_exact(storage.num_outputs).nth(seq_left).unwrap().iter().copied().collect();
        iterator.reset_to_right_level_zero(dimension);
        let seq_right = iterator.seq().unwrap();
        let right_boundary: Vec<f64> = values.chunks_exact(storage.num_outputs).nth(seq_right).unwrap().iter().copied().collect();
        if !iterator.is_leaf()
        {
            iterator.reset_to_level_one(dimension);
            if iterator.seq().is_some()
            {
               LinearDehierarchisation::recurse(values, iterator, dimension, &left_boundary, &right_boundary);
            }
            iterator.reset_to_left_level_zero(dimension);
        }
    }
}

#[derive(Clone, Copy)]
pub struct LinearHierarchisationOperation;

impl HierarchisationOperation for LinearHierarchisationOperation
{        
    #[inline]
    fn hierarchize(&self, node_values: &mut [f64], storage: &SparseGridData) {
        use crate::dynamic::algorithms::sweep;
        let mut func = LinearHierarchisation;
        for d in 0..storage.num_inputs
        {
            sweep::sweep_1d_in_place(&mut func, storage, node_values, d, storage.num_inputs);
        }
    }
    #[inline]
    fn dehierarchize(&self, alpha: &mut [f64], storage: &SparseGridData) {
        use crate::dynamic::algorithms::sweep;
        let mut func = LinearDehierarchisation;
        for d in 0..storage.num_inputs
        {
            sweep::sweep_1d_in_place(&mut func, storage, alpha, d, storage.num_inputs);
        }
    }
}

#[derive(Clone, Copy)]
pub struct LinearBoundaryHierarchisationOperation;

impl HierarchisationOperation for LinearBoundaryHierarchisationOperation
{    
    #[inline]
    fn hierarchize(&self, node_values: &mut [f64], storage: &SparseGridData) {
        use crate::dynamic::algorithms::sweep;
        let mut func = LinearBoundaryHierarchisation(LinearHierarchisation);
        for d in 0..storage.num_inputs
        {
            sweep::sweep_1d_boundary_in_place(&mut func, storage, node_values, d, storage.num_inputs);
        }
    }
    #[inline]
    fn dehierarchize(&self, alpha: &mut [f64], storage: &SparseGridData) {
        use crate::dynamic::algorithms::sweep;
        let mut func = LinearDehierarchisation;
        for d in 0..storage.num_inputs
        {
            sweep::sweep_1d_in_place(&mut func, storage, alpha, d, storage.num_inputs);
        }
    }
}
