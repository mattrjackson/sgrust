use crate::dynamic::storage::{GridPoint, SparseGridData};

use super::dynamic_grid_iterator::GridIteratorT;



///
/// This iterator uses array access to retrieve neighbors more quickly (~10x speedup relative to hash queries)
/// at the expense of pretty significant memory overhead. Eventually the data behind building the iterator
/// will be an optional format
/// 
pub(crate) struct AdjacencyGridIterator<'a>
{
    pub(crate) seq: Option<usize>,    
    storage: &'a SparseGridData,
    index: &'a [u32]
}


impl<'a> AdjacencyGridIterator<'a>
{

    /// This computes the left level one (e.g. index=1 and level=1 for the given dimension).
    #[inline]
    pub fn compute_level_one(&self, dim: usize) -> Option<u32>
    {
        let mut index = self.seq.unwrap();
        let offset = self.offset(dim);
        let level = self.storage.level(index, dim);
        // First we determine which way we iterate up or down the level hierarchy...
        let mut data = &self.storage.adjacency_data[offset+index];
        if level == 0
        {
            if data.inner.has_left_child() 
            {            
                index = self.compute_left_child(dim) as usize;              
            }
        }
        else
        {
            while data.inner.has_parent()
            {            
                index = (index as i64 + data.inner.up()) as usize;
                if self.storage.level(index, dim) == 1
                {
                    break;
                }
                data = &self.storage.adjacency_data[offset + index];
            }
        }
        // If this isn't level one, we have no level one node defined, so return None.
        if  self.storage.level(index, dim) != 1
        {
            return None;
        }
        // now we retrieve the active index and check its value. Again our goal here
        // is to retrieve the node with index 1.
        let active_index = self.storage.index(index, dim) as usize;
        
        if active_index == 0 {
            // Need to step to the right until we find index = 1.
            let mut right_index = index;
            data = &self.storage.adjacency_data[offset + right_index];
            while data.inner.has_right()
            {
                right_index = (right_index as i64 + data.inner.right()) as usize;
                if right_index == 1 
                {
                    // We've found the node with index 1.
                    return Some(right_index as u32)
                }
                data = &self.storage.adjacency_data[offset + right_index];
            }
        }
        else if active_index == 1
        {
            return Some(index as u32);
        }
        else
        {
            let mut left_index = index;
            data = &self.storage.adjacency_data[offset + left_index];
            while data.inner.has_left()
            {
                left_index = (left_index as i64 + data.inner.left()) as usize;
                if left_index == 1
                {
                    // We've found the node with index 1.
                    return Some(left_index as u32)
                }
                data = &self.storage.adjacency_data[offset + left_index];
            }
        }           
        // Otherwise, we didn't find the node with index 1, so return None. 
        None
    }

    /// Compute the index of the left boundary node (e.g. the left zero-level node) for the given dimension.
    #[inline]
    pub fn compute_lzero(&self, dim: usize) -> Option<u32> {      
        let index = self.seq.unwrap();
        let offset = self.offset(dim);  
        Some(self.storage.adjacency_data[offset + index].left_zero)  
    }
    
    #[inline]
    fn compute_left_child(&self, dim: usize) -> u32
    {
        let offset = self.offset(dim);
        let original_node = &self.storage.adjacency_data[offset + self.seq.unwrap()];
        let node_index = self.storage.index(self.seq.unwrap(), dim);
        if original_node.inner.has_child()
        {
            let index = self.seq.unwrap() as i64 + original_node.inner.down();
            let dim_index = self.storage.index(index as usize, dim);
            if dim_index == (2*node_index - 1)
            {
                index as u32
            }
            else
            {               
                (index + self.storage.adjacency_data[offset + index as usize].inner.left()) as u32
            }
        }
        else
        {
            // We check that the left child exists before this is called, so this should never panic.
            panic!("No left child found for node");
        }
    }
    #[inline]
    fn compute_right_child(&self, dim: usize) -> u32
    {
        let offset = self.offset(dim);
        let node_index = self.seq.unwrap();
        let original_node = &self.storage.adjacency_data[offset + node_index];
        let node_index = self.storage.index(node_index, dim);
        if original_node.inner.has_child()
        {
            let index = self.seq.unwrap() as i64 + original_node.inner.down();
            let dim_index = self.storage.index(index as usize, dim);
            if dim_index == 2*node_index + 1
            {
                index as u32
            }
            else
            {   
                // We check that the right child exists before this is called, so this should never panic.           
                (index  + self.storage.adjacency_data[offset + index as usize].inner.right()) as u32
            }
        }
        else
        {
            panic!("No right child found for node");
        }
    }

    /// Compute the index of the right boundary node (e.g. the right zero-level node) for the given dimension.
    #[inline]
    pub fn compute_rzero(&self, dim: usize) -> Option<u32> {
        let mut index = self.seq.unwrap();
        let offset = self.offset(dim);
        while self.storage.adjacency_data[offset + index].inner.has_right()
        {
            index = (index as i64 + self.storage.adjacency_data[offset + index].inner.right()) as usize            
        }        
        if self.storage.index(index, dim) != 1
        {
            None
        }
        else {
            Some(index as u32)    
        }
    }
    pub(crate) fn new(storage: &'a SparseGridData) -> Self
    {
        Self { seq: Some(0), storage, index: &storage.index}
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
        self.storage.adjacency_data[self.offset(dim) + self.seq.unwrap()].inner.has_left_child()
    }
    #[inline]
    #[allow(unused)]
    pub(crate) fn has_right_leaf(&self, dim: usize) -> bool
    {
        self.storage.adjacency_data[self.offset(dim) + self.seq.unwrap()].inner.has_right_child()
    } 
  
}

impl GridIteratorT for AdjacencyGridIterator<'_>
{
    #[inline]
    fn reset_to_level_zero(&mut self) -> bool
    {   
        self.seq = Some(self.storage.adjacency_data.zero_index);
        self.index = &self.storage.index[self.storage.num_inputs*self.storage.adjacency_data.zero_index..];
        true
    }
    #[inline]
    fn reset_to_left_level_zero(&mut self, dim: usize) -> bool
    {
        if let Some(index) = self.compute_lzero(dim)
        {
            self.seq = Some(index as usize);
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
            self.seq = Some(index as usize);      
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
        match self.compute_level_one(dim)
        {
            Some(index) => 
            {
                self.seq = Some(index as usize);
                self.index = &self.storage.index[self.storage.num_inputs*index as usize..];
                true
            },
            None => false,
        }
    }
    #[inline]
    fn left_child(&mut self, dim: usize) -> bool
    {
        match self.storage.adjacency_data[self.offset(dim) + self.seq.unwrap()].inner.has_left_child()
        {
            true => 
            {                
                let index = self.compute_left_child(dim) as usize;
                self.seq = Some(index);
                self.index = &self.storage.index[self.storage.num_inputs*index..];
                true
            },
            false => 
            {                
                false
            },
        }        
    }
    #[inline]
    fn right_child(&mut self, dim: usize) -> bool
    {
        match self.storage.adjacency_data[self.offset(dim) + self.seq.unwrap()].inner.has_right_child()
        {
            true => 
            {                
                let index = self.compute_right_child(dim) as usize;
                self.seq = Some(index);
                self.index = &self.storage.index[self.storage.num_inputs*index..];
                true
            },
            false => 
            {             
                false
            },
        }         
    }
   
    #[inline]
    fn is_leaf(&self) -> bool
    {
        self.storage.flags[self.seq.unwrap()].is_leaf()
    }
    
    fn index(&self) ->  Option<usize> {
        self.seq
    }
    
    fn up(&mut self, dim: usize) -> bool
    {
        let offset = self.offset(dim);
        if self.storage.adjacency_data[offset + self.seq.unwrap()].inner.has_parent()
        {
            let index = (self.seq.unwrap() as i64 + self.storage.adjacency_data[offset + self.seq.unwrap()].inner.up()) as usize;
            self.seq = Some(index);
            self.index = &self.storage.index[self.storage.num_inputs*index..];
            true
        }
        else
        {
            false
        }
    }
    
    fn is_inner_point(&self) -> bool {
        self.storage.flags[self.seq.unwrap()].is_inner()
    }
    
    fn point(&self) -> &GridPoint {
        panic!("Adjacency Iterator doesn't return a GridPoint for Dynamic")
        //&self.storage[self.seq.unwrap()]
    }
    #[inline]
    fn point_index(&self, dim: usize) -> u32
    {
        self.index[dim]
    }
}