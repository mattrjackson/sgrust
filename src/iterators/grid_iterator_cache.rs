
use crate::storage::linear_grid::{GridPoint, SparseGridData};

use super::grid_iterator::GridIteratorT;



///
/// This iterator uses array access to retrieve neighbors more quickly (~10x speedup relative to hash queries)
/// at the expense of pretty significant memory overhead. Eventually the data behind building the iterator
/// will be an optional format
/// 
pub(crate) struct AdjacencyGridIterator<'a, const D: usize>
{
    pub(crate) index: Option<usize>,    
    storage: &'a SparseGridData<D>,
}


impl<'a, const D: usize> AdjacencyGridIterator<'a, D>
{

    /// This computes the left level one (e.g. index=1 and level=1 for the given dimension).
    #[inline]
    pub fn compute_level_one(&self, dim: usize) -> Option<u32>
    {
        let mut index = self.index.unwrap();
        let offset = self.offset(dim);
        let level = self.storage[index].level[dim];
        // First we determine which way we iterate up or down the level hierarchy...
        if level == 0
        {
            if self.storage.adjacency_data[offset + index].has_left_child() 
            {            
                index = self.compute_left_child(dim) as usize;              
            }
        }
        else
        {
            while self.storage.adjacency_data[offset + index].has_parent()
            {            
                index = (index as i64 + self.storage.adjacency_data[offset + index].up()) as usize;
                if self.storage[index].level[dim] == 1
                {
                    break;
                }
            }
        }
        // Now get the node for current level. 
        let node= &self.storage[index];
        // If this isn't level one, we have no level one node defined, so return None.
        if node.level[dim] != 1
        {
            return None;
        }
        // now we retrieve the active index and check its value. Again our goal here
        // is to retrieve the node with index 1.
        let active_index = self.storage[index].index[dim] as usize;
        match active_index.cmp(&1)
        {
            std::cmp::Ordering::Less => {
                // Need to step to the right until we find index = 1.
                let mut right_index = index;
                while self.storage.adjacency_data[self.offset(dim) + right_index].has_right()
                {
                    right_index = (right_index as i64 + self.storage.adjacency_data[self.offset(dim) + right_index].right()) as usize;
                    if right_index == 1 
                    {
                        // We've found the node with index 1.
                        return Some(right_index as u32)
                    }
                }
            },
            std::cmp::Ordering::Equal => 
            {
                return Some(index as u32);
            },
            std::cmp::Ordering::Greater => 
            {
                let mut left_index = index;
                while self.storage.adjacency_data[self.offset(dim) + left_index].has_left()
                {
                    left_index = (left_index as i64 + self.storage.adjacency_data[self.offset(dim) + left_index].left()) as usize;
                    if left_index == 1
                    {
                        // We've found the node with index 1.
                        return Some(left_index as u32)
                    }
                }
            },
        }           
        // Otherwise, we didn't find the node with index 1, so return None. 
        None
    }

    /// Compute the index of the left boundary node (e.g. the left zero-level node) for the given dimension.
    #[inline]
    pub fn compute_lzero(&self, dim: usize) -> Option<u32> {        
        let mut index = self.index.unwrap();
        let offset = self.offset(dim);
        while self.storage.adjacency_data[offset + index].has_left()
        {
            index = (index as i64 + self.storage.adjacency_data[offset + index].left()) as usize
        }
        if self.storage[index].index[dim] != 0
        {
            None
        }
        else {
            Some(index as u32)    
        }        
    }
    #[inline]
    fn compute_left_child(&self, dim: usize) -> u32
    {
        let original_node = &self.storage.adjacency_data[self.offset(dim) + self.index.unwrap()];
        let node_index = self.storage[self.index.unwrap()].index[dim];
        if original_node.has_child()
        {
            let index = self.index.unwrap() as i64 + original_node.down();
            let dim_index = self.storage[index as usize].index[dim];
            if dim_index == (2*node_index - 1)
            {
                index as u32
            }
            else
            {               
                (index + self.storage.adjacency_data[self.offset(dim) + index as usize].left()) as u32
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
        let node_index = self.index.unwrap();
        let original_node = &self.storage.adjacency_data[self.offset(dim) + node_index];
        let node_index = self.storage[node_index].index[dim];
        if original_node.has_child()
        {
            let index = self.index.unwrap() as i64 + original_node.down();
            let dim_index = self.storage[index as usize].index[dim];
            if dim_index == 2*node_index + 1
            {
                index as u32
            }
            else
            {   
                // We check that the right child exists before this is called, so this should never panic.           
                (index  + self.storage.adjacency_data[ self.offset(dim) + index as usize].right()) as u32
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
        let mut index = self.index.unwrap();
        let offset = self.offset(dim);
        while self.storage.adjacency_data[ offset + index].has_right()
        {
            index = (index as i64 + self.storage.adjacency_data[offset + index].right()) as usize            
        }        
        if self.storage[index].index[dim] != 1
        {
            None
        }
        else {
            Some(index as u32)    
        }
    }
    pub(crate) fn new(storage: &'a SparseGridData<D>) -> Self
    {
        Self { index: Some(0), storage}
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
        self.storage.adjacency_data[self.offset(dim) + self.index.unwrap()].has_left_child()
    }
    #[inline]
    #[allow(unused)]
    pub(crate) fn has_right_leaf(&self, dim: usize) -> bool
    {
        self.storage.adjacency_data[self.offset(dim) + self.index.unwrap()].has_right_child()
    } 
  
}

impl<const D: usize> GridIteratorT<D> for AdjacencyGridIterator<'_, D>
{
    #[inline]
    fn reset_to_level_zero(&mut self) -> bool
    {   
        self.index = Some(self.storage.adjacency_data.zero_index);
        true
    }
    #[inline]
    fn reset_to_left_level_zero(&mut self, dim: usize) -> bool
    {
        if let Some(index) = self.compute_lzero(dim)
        {
            self.index = Some(index as usize);
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
            self.index = Some(index as usize);            
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
                self.index = Some(index as usize);
                true
            },
            None => false,
        }
    }
    #[inline]
    fn left_child(&mut self, dim: usize) -> bool
    {
        match self.storage.adjacency_data[self.offset(dim) + self.index.unwrap()].has_left_child()
        {
            true => 
            {                
                self.index = Some(self.compute_left_child(dim) as usize);
                
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
        match self.storage.adjacency_data[self.offset(dim) + self.index.unwrap()].has_right_child()
        {
            true => 
            {                
                self.index = Some(self.compute_right_child(dim) as usize);
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
        self.storage[self.index.unwrap()].is_leaf()
    }
    
    fn index(&self) ->  Option<usize> {
        self.index
    }
    
    fn up(&mut self, dim: usize) -> bool
    {
        if self.storage.adjacency_data[self.offset(dim) + self.index.unwrap()].has_parent()
        {
            self.index = Some((self.index.unwrap() as i64 + self.storage.adjacency_data[self.offset(dim) + self.index.unwrap()].up()) as usize);
            true
        }
        else
        {
            false
        }
    }
    
    fn is_inner_point(&self) -> bool {
        self.storage[self.index.unwrap()].is_inner_point()
    }
    
    fn point(&self) -> &GridPoint<D> {
        &self.storage[self.index.unwrap()]
    }
}