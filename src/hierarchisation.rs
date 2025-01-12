use crate::{algorithms::sweep::SweepFunction, iterators::grid_iterator::{GridIterator, GridIteratorT}, storage::linear_grid::SparseGridStorage};


pub trait HierarchisationOperation<const D: usize, const DIM_OUT: usize> : Copy
{
    fn hierarchize(&self, node_values: &mut [[f64; DIM_OUT]], storage: &SparseGridStorage<D>); 
    fn dehierarchize(&self, alpha: &mut [[f64; DIM_OUT]], storage: &SparseGridStorage<D>); 
}

pub struct LinearHierarchisation<const D: usize, const DIM_OUT: usize>;
impl<const D: usize, const DIM_OUT: usize> LinearHierarchisation<D, DIM_OUT>
{
    fn recurse(values: &mut [[f64; DIM_OUT]], iterator: &mut GridIterator<D>, dimension: usize, left_value: [f64; DIM_OUT], right_value: [f64; DIM_OUT])
    {
        if let Some(seq) = iterator.seq()
        {
            let mid_value = values[seq];
            if !iterator.is_leaf()
            {            
                iterator.left_child(dimension);
                if let Some(_seq) = iterator.seq()
                {                    
                    Self::recurse(values, iterator, dimension, left_value, mid_value);                    
                }
                iterator.step_right(dimension);
                if let Some(_seq) = iterator.seq()                
                {
                    Self::recurse(values, iterator, dimension, mid_value, right_value);
                }
                iterator.up(dimension);
            }
            for i in 0..DIM_OUT
            {
                values[seq][i] = mid_value[i] - 0.5 * (left_value[i] + right_value[i]);
            }
        }
    }
}

impl< const D: usize, const DIM_OUT: usize> SweepFunction<D, [f64; D], [f64; DIM_OUT]> for LinearHierarchisation<D, DIM_OUT>
{
    fn execute_in_place(&mut self, values: &mut [[f64; DIM_OUT]], iterator: &mut GridIterator<D>, _storage: &SparseGridStorage<D>, dimension: usize) {
        Self::recurse(values, iterator, dimension, [0.0; DIM_OUT], [0.0; DIM_OUT]);
    }
}


pub struct LinearBoundaryHierarchisation<const D: usize, const DIM_OUT: usize>(LinearHierarchisation<D, DIM_OUT>);

impl<const D: usize, const DIM_OUT: usize> SweepFunction<D, [f64; D], [f64; DIM_OUT]> for LinearBoundaryHierarchisation<D, DIM_OUT>
{
    #[inline]
    fn execute_in_place(&mut self, values: &mut [[f64; DIM_OUT]], iterator: &mut crate::iterators::grid_iterator::GridIterator<D>, _storage: &SparseGridStorage<D>, dimension: usize) {
        // left boundary
        iterator.reset_to_left_level_zero(dimension);
        let seq_left = iterator.seq().unwrap();
        let left_boundary = values[seq_left];
        iterator.reset_to_right_level_zero(dimension);
        let seq_right = iterator.seq().unwrap();       
        let right_boundary = values[seq_right];
        if !iterator.is_leaf()
        {
            iterator.reset_to_level_one(dimension);
            if iterator.seq().is_some()
            {
                LinearHierarchisation::recurse(values, iterator, dimension, left_boundary, right_boundary);
            }
            iterator.reset_to_left_level_zero(dimension);
        }
    }
}

pub struct LinearDehierarchisation<const D: usize, const DIM_OUT: usize>;

impl<const D: usize, const DIM_OUT: usize> LinearDehierarchisation<D, DIM_OUT>
{
    #[inline]
    fn recurse(values: &mut [[f64; DIM_OUT]], iterator: &mut GridIterator<D>, dimension: usize, left_value: [f64; DIM_OUT], right_value: [f64; DIM_OUT])
    {
        if let Some(seq) = iterator.seq()
        {
            let mut mid_value = values[seq];
            for i in 0..DIM_OUT
            {
                mid_value[i] += (left_value[i] + right_value[i]) * 0.5;
            }
            values[seq] = mid_value;
            if !iterator.is_leaf()
            {
                iterator.left_child(dimension);
                if iterator.seq().is_some()
                {
                    Self::recurse(values, iterator, dimension, left_value, mid_value);
                }
                iterator.step_right(dimension);
                if iterator.seq().is_some()
                {
                    Self::recurse(values, iterator, dimension, mid_value, right_value);
                }
                iterator.up(dimension);
            }
        }
    }
}
impl<const D: usize, const DIM_OUT: usize> SweepFunction<D, [f64; D], [f64; DIM_OUT]> for LinearDehierarchisation<D, DIM_OUT>
{
    #[inline]
    fn execute_in_place(&mut self, values: &mut [[f64; DIM_OUT]], iterator: &mut GridIterator<D>, _storage: &SparseGridStorage<D>, dimension: usize) {
        Self::recurse(values, iterator, dimension, [0.0; DIM_OUT], [0.0; DIM_OUT]);
    }
}

pub struct LinearBoundaryDehierarchisation<const D: usize, const DIM_OUT: usize>(LinearDehierarchisation<D, DIM_OUT>);

impl<const D: usize, const DIM_OUT: usize> SweepFunction<D, [f64; D], [f64; DIM_OUT]> for LinearBoundaryDehierarchisation<D, DIM_OUT>
{
    #[inline]
    fn execute_in_place(&mut self, values: &mut [[f64; DIM_OUT]], iterator: &mut crate::iterators::grid_iterator::GridIterator<D>, _storage: &SparseGridStorage<D>, dimension: usize) {
        // left boundary
        iterator.reset_to_left_level_zero(dimension);
        let seq_left = iterator.seq().unwrap();
        let left_boundary = values[seq_left];
        iterator.reset_to_right_level_zero(dimension);
        let seq_right = iterator.seq().unwrap();
        let right_boundary = values[seq_right];
        if !iterator.is_leaf()
        {
            iterator.reset_to_level_one(dimension);
            if iterator.seq().is_some()
            {
               LinearDehierarchisation::recurse(values, iterator, dimension, left_boundary, right_boundary);
            }
            iterator.reset_to_left_level_zero(dimension);
        }
    }
}

#[derive(Clone, Copy)]
pub struct LinearHierarchisationOperation<const D: usize, const DIM_OUT: usize>;

impl<const D: usize, const DIM_OUT: usize> HierarchisationOperation<D, DIM_OUT> for LinearHierarchisationOperation<D, DIM_OUT>
{        
    #[inline]
    fn hierarchize(&self, node_values: &mut [[f64; DIM_OUT]], storage: &SparseGridStorage<D>) {
        use crate::algorithms::sweep;
        let mut func = LinearHierarchisation;
        for d in 0..D
        {
            sweep::sweep_1d_in_place(&mut func, storage, node_values, d);
        }
    }
    #[inline]
    fn dehierarchize(&self, alpha: &mut [[f64; DIM_OUT]], storage: &SparseGridStorage<D>) {
        use crate::algorithms::sweep;
        let mut func = LinearDehierarchisation;
        for d in 0..D
        {
            sweep::sweep_1d_in_place(&mut func, storage, alpha, d);
        }
    }
}

#[derive(Clone, Copy)]
pub struct LinearBoundaryHierarchisationOperation<const D: usize, const DIM_OUT: usize>;

impl<const D: usize, const DIM_OUT: usize> HierarchisationOperation<D, DIM_OUT> for LinearBoundaryHierarchisationOperation<D, DIM_OUT>
{    
    #[inline]
    fn hierarchize(&self, node_values: &mut [[f64; DIM_OUT]], storage: &SparseGridStorage<D>) {
        use crate::algorithms::sweep;
        let mut func = LinearBoundaryHierarchisation(LinearHierarchisation);
        for d in 0..D
        {
            sweep::sweep_1d_boundary_in_place(&mut func, storage, node_values, d);
        }
    }
    #[inline]
    fn dehierarchize(&self, alpha: &mut [[f64; DIM_OUT]], storage: &SparseGridStorage<D>) {
        use crate::algorithms::sweep;
        let mut func = LinearDehierarchisation;
        for d in 0..D
        {
            sweep::sweep_1d_in_place(&mut func, storage, alpha, d);
        }
    }
}
