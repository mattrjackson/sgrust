use std::collections::HashMap;

use rustc_hash::FxBuildHasher;
use serde::{Deserialize, Serialize};

use crate::storage::linear_grid::{GridPoint, SparseGridStorage};

use super::grid_iterator::GridIterator;

#[derive(Serialize, Deserialize)]
#[derive(Copy, Clone, Default)]
pub(crate) struct NodeProperties
{
    pub left: Option<u32>,
    pub right: Option<u32>,
    pub left_child: Option<u32>,
    pub right_child: Option<u32>,
    pub up: Option<u32>,
    pub is_inner: bool,
    pub is_leaf: bool,
}

impl NodeProperties
{
    pub fn is_complete(&self) -> bool
    {
        self.left.is_some() && self.right.is_some()
        && self.left_child.is_some()&& self.right_child.is_some() 
        && self.up.is_some()
    }
}
#[derive(Serialize, Deserialize, Clone)]
pub(crate) struct GridIteratorData<const D: usize>
{
    array: Vec<NodeProperties>,
    zero_level_node: usize,
    storage_len: usize,
}

impl<const D: usize> GridIteratorData<D>
{
    pub(crate) fn new(storage: & SparseGridStorage<D> ) -> Self
    {
        let mut iterator = GridIterator::new(storage);
        let mut array = vec![NodeProperties::default(); D*storage.len()];
        
        for dim in 0..D
        {
            for i in 0..storage.len()
            {
                update_indices(&mut iterator, &mut array, &storage.map, i, dim);     
            }
        }
        iterator.reset_to_level_zero();
        let zero_idx = iterator.index();
        let zero_index = *storage.map.get(zero_idx).unwrap_or(&0);
        GridIteratorData { zero_level_node: zero_index,  array, storage_len: storage.len()}
    }
}

///
/// This iterator uses array access to retrieve neighbors more quickly (~10x speedup relative to hash queries)
/// at the expense of pretty significant memory overhead. Eventually the data behind building the iterator
/// will be an optional format
/// 
pub(crate) struct GridIteratorWithCache<'b, const D: usize>
{
    pub(crate) index: usize,    
    data: &'b GridIteratorData<D>
}

fn update_indices<const D: usize>(iterator: &mut GridIterator<D>, array:&mut [NodeProperties],
    map: &HashMap<GridPoint<D>, usize, FxBuildHasher>, seq: usize, dim: usize)
{
    let offset =  dim * map.len();
    let active_index = offset + seq;
    if !array[active_index].is_complete()
    {
        let node_index= iterator.storage[seq]; 
        iterator.set_index(node_index);
        iterator.step_left(dim);
        if iterator.valid_seq()
        {
            let left_index = map[iterator.index()];
            // assign left node            
            array[active_index].left = Some(left_index as u32);
            // assign right node for left of current node...
            array[offset + left_index].right = Some(seq as u32);
        }
        iterator.set_index(node_index);
        iterator.step_right(dim);     
        if iterator.valid_seq()
        {
            // assign right node
            let right_index = map[iterator.index()];
            array[active_index].right = Some(right_index as u32);
            // assign left node for right of current node...
            array[offset + right_index].left = Some(seq as u32);            
        }        
        iterator.set_index(node_index);
        iterator.up(dim);       
        if iterator.valid_seq()  // parent exists
        {
            let parent_index = map[iterator.index()];
            // assign parent
            array[active_index].up = Some(parent_index as u32);            
        }
        iterator.set_index(node_index);
        // get left child
        iterator.left_child(dim);
        if iterator.valid_seq()
        {
            // assign left child
            let lc_index = map[iterator.index()];            
            array[active_index].left_child = Some(lc_index as u32);                        
        }
        
        iterator.set_index(node_index);
        // get right child        
        iterator.right_child(dim);
        if iterator.valid_seq()
        {
            // assign right child
            let rc_index = map[iterator.index()];      
            array[active_index].right_child = Some(rc_index as u32);            
        } 
        iterator.reset_to_left_level_zero(dim);
        // Handle left level zero
        if let Some(&index) = map.get(iterator.index()) {             
            // we know what the correct index is...
            let lzero = index as u32;
            // first let's find the leftmost node linked in our data structure.
            let mut left_index = seq as u32;
            while let Some(index) = array[offset + left_index as usize].left
            {
                left_index = index;
            }
            // If the leftmost node isn't the boundary, we need to  update our data structure
            // such that if it has boundaries, we set its left boundary neighbor.
            if left_index != lzero
            {
                array[offset + left_index as usize].left = Some(lzero);
            }
        } 
        iterator.reset_to_right_level_zero(dim);
        if let Some(&index) = map.get(iterator.index())
        {
            // first let's find the rightmost node linked in our data structure.
            let rzero = index as u32;
            let mut right_index = seq as u32;
            while let Some(index) = array[offset + right_index as usize].right
            {
                right_index = index;
            }
            // If the rightmost node isn't the boundary, we need to update our data structure
            // such that if it has boundaries, we set its right boundary neighbor.
            if right_index != rzero
            {
                array[offset + right_index as usize].right = Some(rzero);
            }
        }       
    }

    
}

impl<'b, const D: usize> GridIteratorWithCache<'b,  D>
{

    /// This computes the left level one (e.g. index=1 and level=1 for the given dimension).
    #[inline]
    pub fn compute_level_one(&self, dim: usize, storage: &SparseGridStorage<D>) -> Option<u32>
    {
        let mut index = self.index;
        let level = storage[index].level[dim];
        // First we determine which way we iterate up or down the level hierarchy...
        if level == 0
        {
            while let Some(parent_index) = self.data.array[self.offset(dim) + index].left_child 
            {            
                index = parent_index as usize;
                if storage[index].level[dim] == 1
                {
                    break;
                }
            }
        }
        else
        {
            while let Some(parent_index) = self.data.array[self.offset(dim) + index].up 
            {            
                index = parent_index as usize;
                if storage[index].level[dim] == 1
                {
                    break;
                }
            }
        }
        // Now get the node for current level. 
        let node= &storage.list()[index];
        // If this isn't level one, we have no level one node defined, so return None.
        if node.level[dim] != 1
        {
            return None;
        }
        // now we retrieve the active index and check its value. Again our goal here
        // is to retrieve the node with index 1.
        let active_index = storage.list()[index].index[dim] as usize;
        match active_index.cmp(&1)
        {
            std::cmp::Ordering::Less => {
                // Need to step to the right until we find index = 1.
                let mut right_index = index;
                while let Some(index) = self.data.array[self.offset(dim) +  right_index].right
                {
                    right_index = index as usize;
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
                while let Some(index) = self.data.array[self.offset(dim) +  left_index].right
                {
                    left_index = index as usize;
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
    pub fn compute_lzero(&self, dim: usize, storage: &SparseGridStorage<D>) -> Option<u32> {        
        let mut index = self.index;
        while let Some(left_node) = self.data.array[self.offset(dim) + index].left
        {
            index = left_node as usize
        }
        if storage[index].index[dim] != 0
        {
            None
        }
        else {
            Some(index as u32)    
        }
        
    }

    /// Compute the index of the right boundary node (e.g. the right zero-level node) for the given dimension.
    #[inline]
    pub fn compute_rzero(&self, dim: usize, storage: &SparseGridStorage<D>) -> Option<u32> {
        let mut index = self.index;
        while let Some(right_node) = self.data.array[self.offset(dim) + index].right
        {
            index = right_node as usize
        }        
        if storage[index].index[dim] != 1
        {
            None
        }
        else {
            Some(index as u32)    
        }
    }
    pub(crate) fn new(data: &'b GridIteratorData<D>) -> Self
    {
        Self { index: 0, data }
    }
    #[inline(always)]
    fn offset(&self, dim: usize) -> usize
    {
        dim * self.data.storage_len
    }
    #[inline]
    pub(crate) fn reset_to_level_zero(&mut self)
    {   
        self.index = self.data.zero_level_node;
    }
    #[inline]
    pub(crate) fn reset_to_left_level_zero(&mut self, dim: usize, storage: &SparseGridStorage<D>) -> bool
    {
        if let Some(index) = self.compute_lzero(dim, storage)
        {
            self.index = index as usize;
            true
        }
        else
        {
            false
        }        
    }
    #[inline]
    pub(crate) fn reset_to_right_level_zero(&mut self, dim: usize, storage: &SparseGridStorage<D>) -> bool
    {
        if let Some(index) = self.compute_rzero(dim, storage)
        {           
            self.index = index as usize;            
            true
        }
        else
        {
            false
        }     
    } 

    #[inline]
    #[allow(unused)]
    pub(crate) fn is_valid(&self) -> bool
    {
        self.index < self.data.storage_len
    }

    #[inline]
    pub(crate) fn reset_to_level_one(&mut self, dim: usize, storage: &SparseGridStorage<D>) -> bool
    {
        match self.compute_level_one(dim, storage)
        {
            Some(index) => 
            {
                self.index = index as usize;
                true
            },
            None => false,
        }
    }
    #[inline]
    pub(crate) fn left_child(&mut self, dim: usize) -> bool
    {
        match self.data.array[self.offset(dim) + self.index].left_child
        {
            Some(index) => 
            {
                self.index = index as usize;
                true
            },
            None => 
            {                
                false
            },
        }        
    }
    #[inline]
    pub(crate) fn right_child(&mut self, dim: usize) -> bool
    {
        match self.data.array[self.offset(dim) + self.index].right_child
        {
            Some(index) => 
            {
                self.index = index as usize;
                true
            },
            None => 
            {
                false
            }
        }        
    }
   
    #[inline]
    pub(crate) fn is_leaf(&self, storage: &SparseGridStorage<D>) -> bool
    {
        storage[self.index].is_leaf()
    }
    #[inline]
    #[allow(unused)]
    pub(crate) fn is_left_leaf(&self, dim: usize) -> bool
    {
        self.data.array[self.offset(dim) + self.index].left_child.is_some()
    }
    #[inline]
    #[allow(unused)]
    pub(crate) fn is_right_leaf(&self, dim: usize) -> bool
    {
        self.data.array[self.offset(dim) + self.index].right_child.is_some()
    } 
    #[inline]
    #[allow(unused)]
    pub(crate) fn get_grid_depth(&mut self, dim: usize) -> usize
    {
        let mut depth = 1;
        let orig_index = self.index;
        loop
        {            
            if self.is_left_leaf(dim) || self.is_right_leaf(dim)
            {
                depth += 1;
            }
            else
            {                
               break;
            }
        }
        self.index = orig_index;
        depth
    }
}