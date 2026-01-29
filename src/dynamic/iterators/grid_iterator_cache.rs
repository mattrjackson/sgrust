use crate::dynamic::storage::{GridPoint, SparseGridData};

use super::dynamic_grid_iterator::GridIteratorT;



///
/// This iterator uses array access to retrieve neighbors more quickly (~10x speedup relative to hash queries)
/// at the expense of pretty significant memory overhead. Eventually the data behind building the iterator
/// will be an optional format
/// 
pub(crate) struct AdjacencyGridIterator<'a>
{
    pub(crate) seq: usize,    
    storage: &'a SparseGridData,
    index: &'a [u32]
}


impl<'a> AdjacencyGridIterator<'a>
{
    /// Compute the index of the left boundary node (e.g. the left zero-level node) for the given dimension.
    #[inline]
    pub fn compute_lzero(&self, dim: usize) -> Option<u32> {      
        let index = self.seq;
        let offset = self.offset(dim);  
        Some(self.storage.adjacency_data.left_zero[offset + index])  
    }

    /// Compute the index of the right boundary node (e.g. the right zero-level node) for the given dimension.
    #[inline]
    pub fn compute_rzero(&self, dim: usize) -> Option<u32> {
        let index = self.seq;
        let offset = self.offset(dim);  
        Some(self.storage.adjacency_data.right_zero[offset + index])  
    }
    #[inline]
    pub(crate) fn new(storage: &'a SparseGridData) -> Self
    {
        Self { seq: 0, storage, index: &storage.index}
    }
    #[inline(always)]
    fn offset(&self, dim: usize) -> usize
    {
        dim * self.storage.len()
    }
    

    
    #[inline]
    #[allow(unused)]
    pub(crate) fn has_left_leaf(&self, dim: usize) -> bool
    {
        self.storage.adjacency_data[self.offset(dim) + self.seq].has_left_child()
    }
    #[inline]
    #[allow(unused)]
    pub(crate) fn has_right_leaf(&self, dim: usize) -> bool
    {
        self.storage.adjacency_data[self.offset(dim) + self.seq].has_right_child()
    } 
  
}

impl GridIteratorT for AdjacencyGridIterator<'_>
{
    #[inline]
    fn reset_to_level_zero(&mut self) -> bool
    {   
        self.seq = self.storage.adjacency_data.zero_index;
        self.index = &self.storage.index[self.storage.num_inputs*self.storage.adjacency_data.zero_index..];
        true
    }
    #[inline]
    fn reset_to_left_level_zero(&mut self, dim: usize) -> bool
    {
        if let Some(index) = self.compute_lzero(dim)
        {
            self.seq = index as usize;
            self.index = &self.storage.index[self.storage.num_inputs*index as usize..];
            true
        }
        else
        {
            false
        }        
    }
    #[inline]
    fn reset_to_right_level_zero(&mut self, dim: usize) -> bool
    {
        if let Some(index) = self.compute_rzero(dim)
        {           
            self.seq = index as usize;      
            self.index = &self.storage.index[self.storage.num_inputs*index as usize..];      
            true
        }
        else
        {
            false
        }     
    } 

    #[inline]
    fn reset_to_level_one(&mut self, dim: usize) -> bool
    {      
        let index = self.storage.adjacency_data[self.offset(dim) + self.seq].level_one();          
        if index != u32::MAX {
            self.seq = index as usize;
            self.index = &self.storage.index[self.storage.num_inputs * index as usize..];
            true
        } else {
            false
        }
    }
    #[inline]
    fn left_child(&mut self, dim: usize) -> bool
    {
        let adj = &self.storage.adjacency_data[self.offset(dim) + self.seq];
        if adj.has_left_child()
        {
            let index = (self.seq as i64 + adj.down_left()) as usize;
            self.seq = index;
            self.index = &self.storage.index[self.storage.num_inputs*index..];
            true
        }
        else
        {
            false
        }
    }
    #[inline]
    fn right_child(&mut self, dim: usize) -> bool
    {
        let adj = &self.storage.adjacency_data[self.offset(dim) + self.seq];
        if adj.has_right_child()
        {
            let index = (self.seq as i64 + adj.down_right()) as usize;
            self.seq = index;
            self.index = &self.storage.index[self.storage.num_inputs*index..];
            true
        }
        else
        {
            false
        }
    }
   
    #[inline]
    fn is_leaf(&self) -> bool
    {
        self.storage.flags[self.seq].is_leaf()
    }
    
    fn index(&self) ->  Option<usize> {
        Some(self.seq)
    }
    
    fn up(&mut self, dim: usize) -> bool
    {
        let offset = self.offset(dim);
        if self.storage.adjacency_data[offset + self.seq].has_parent()
        {
            let index = (self.seq as i64 + self.storage.adjacency_data[offset + self.seq].up()) as usize;
            self.seq = index;
            self.index = &self.storage.index[self.storage.num_inputs*index..];
            true
        }
        else
        {
            false
        }
    }
    
    fn is_inner_point(&self) -> bool {
        self.storage.flags[self.seq].is_inner()
    }
    
    fn point(&self) -> &GridPoint {
        panic!("Adjacency Iterator doesn't return a GridPoint for Dynamic")
        //&self.storage[self.seq]
    }
    #[inline]
    fn point_index(&self, dim: usize) -> u32
    {
        self.index[dim]
    }
}