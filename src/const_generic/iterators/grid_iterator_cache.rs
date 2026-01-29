
use crate::const_generic::storage::{GridPoint, SparseGridData};

use super::grid_iterator::GridIteratorT;



///
/// This iterator uses array access to retrieve neighbors more quickly (~10x speedup relative to hash queries)
/// at the expense of pretty significant memory overhead. Eventually the data behind building the iterator
/// will be an optional format
/// 
pub(crate) struct AdjacencyGridIterator<'a, const D: usize>
{
    pub(crate) seq: usize,    
    storage: &'a SparseGridData<D>,
    is_leaf: bool,
}


impl<'a, const D: usize> AdjacencyGridIterator<'a, D>
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
    pub(crate) fn new(storage: &'a SparseGridData<D>) -> Self
    {
        Self { seq: 0, storage, is_leaf: storage[0].is_leaf()}
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

impl<const D: usize> GridIteratorT<D> for AdjacencyGridIterator<'_, D>
{
    #[inline]
    fn reset_to_level_zero(&mut self) -> bool
    {   
        self.seq = self.storage.adjacency_data.zero_index;
        self.is_leaf = self.storage[self.storage.adjacency_data.zero_index].is_leaf();
        true
    }
    #[inline]
    fn reset_to_left_level_zero(&mut self, dim: usize) -> bool
    {
        if let Some(index) = self.compute_lzero(dim)
        {
            self.seq = index as usize;
            self.is_leaf = self.storage[index as usize].is_leaf();
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
            self.is_leaf = self.storage[index as usize].is_leaf();
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
        if index == u32::MAX
        {
            return false;
        }
        self.seq = index as usize;
        self.is_leaf = self.storage[index as usize].is_leaf();
        true
    }
    #[inline]
    fn left_child(&mut self, dim: usize) -> bool
    {
        let adj = &self.storage.adjacency_data[self.offset(dim) + self.seq];
        if adj.has_left_child()
        {
            let index = (self.seq as i64 + adj.down_left()) as usize;
            self.seq = index;
            self.is_leaf = self.storage[index].is_leaf();                
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
            self.is_leaf = self.storage[index].is_leaf();                
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
        self.is_leaf
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
            self.is_leaf = self.storage[index].is_leaf();
            true
        }
        else
        {
            false
        }
    }
    
    fn is_inner_point(&self) -> bool {
        self.storage[self.seq].is_inner_point()
    }
    
    fn point(&self) -> &GridPoint<D> {
        &self.storage[self.seq]
    }
}