use std::collections::HashMap;

use rustc_hash::FxBuildHasher;
use serde::{Deserialize, Serialize};

use crate::storage::linear_grid::{GridPoint, SparseGridStorage};

use super::grid_iterator::GridIterator;

#[derive(Serialize, Deserialize)]
#[derive(Copy, Clone, Default)]
pub(crate) struct NodeProperties
{
    /// Left neighbor
    pub(crate) left: u32,
    /// Right neighbor
    pub(crate) right: u32,
    /// Child (could be either left or right child)
    pub(crate) down: u32,
    /// Parent node
    pub(crate) up: u32,
    /// Flags used to compactly store state of presence of neighbors (e.g. left/right/up/down/left_child/right_child).
    flags: u8,
}

impl NodeProperties
{
    const HAS_LEFT: usize = 0;
    const HAS_RIGHT: usize = 1;
    const HAS_LEFT_CHILD: usize = 2;
    const HAS_RIGHT_CHILD: usize = 3;
    const HAS_PARENT: usize = 4;
    const HAS_CHILD: usize = 5;
    #[inline]
    fn set_flag<const INDEX: usize>(&mut self, value: bool) {
        if value {
            self.flags |= 1 << INDEX;
        } else {
            self.flags &= !(1 << INDEX);
        }
    }

    #[inline]
    fn get_flag<const INDEX: usize>(&self) -> bool {
        (self.flags & (1 << INDEX)) != 0
    }

    #[inline]
    pub fn has_left(&self) -> bool {
        self.get_flag::<{ Self::HAS_LEFT }>()
    }

    #[inline]
    pub fn set_has_left(&mut self, value: bool) {
        self.set_flag::<{ Self::HAS_LEFT }>(value);
    }

    #[inline]
    pub fn has_right(&self) -> bool {
        self.get_flag::<{ Self::HAS_RIGHT }>()
    }

    #[inline]
    pub fn set_has_right(&mut self, value: bool) {
        self.set_flag::<{ Self::HAS_RIGHT }>(value);
    }

    #[inline]
    pub fn has_left_child(&self) -> bool {
        self.get_flag::<{ Self::HAS_LEFT_CHILD }>()
    }

    #[inline]
    pub fn set_has_left_child(&mut self, value: bool) {
        self.set_flag::<{ Self::HAS_LEFT_CHILD }>(value);
    }

    #[inline]
    pub fn has_right_child(&self) -> bool {
        self.get_flag::<{ Self::HAS_RIGHT_CHILD }>()
    }

    #[inline]
    pub fn set_has_right_child(&mut self, value: bool) {
        self.set_flag::<{ Self::HAS_RIGHT_CHILD }>(value);
    }

    #[inline]
    pub fn has_parent(&self) -> bool {
        self.get_flag::<{ Self::HAS_PARENT }>()
    }

    #[inline]
    pub fn set_has_parent(&mut self, value: bool) {
        self.set_flag::<{ Self::HAS_PARENT }>(value);
    }

    #[inline]
    pub fn has_child(&self) -> bool {
        self.get_flag::<{ Self::HAS_CHILD }>()
    }

    #[inline]
    pub fn set_has_child(&mut self, value: bool) {
        self.set_flag::<{ Self::HAS_CHILD }>(value);
    }

    #[inline]
    pub fn is_complete(&self) -> bool {
        self.has_left() && self.has_right() && self.has_left_child() && self.has_right_child() && self.has_parent()
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
            array[active_index].left = left_index as u32;
            array[active_index].set_has_left(true);
            // assign right node for left of current node...
            array[offset + left_index].right = seq as u32;
            array[offset + left_index].set_has_right(true);
        }
        iterator.set_index(node_index);
        iterator.step_right(dim);     
        if iterator.valid_seq()
        {
            // assign right node
            let right_index = map[iterator.index()];
            array[active_index].right = right_index as u32;
            array[active_index].set_has_right(true);
            // assign left node for right of current node...
            array[offset + right_index].left = seq as u32;   
            array[offset + right_index].set_has_left(true);
        }        
        iterator.set_index(node_index);
        iterator.up(dim);       
        if iterator.valid_seq()  // parent exists
        {
            let parent_index = map[iterator.index()];
            // assign parent
            array[active_index].up = parent_index as u32;        
            array[active_index].set_has_parent(true);
        }
        iterator.set_index(node_index);
        // get left child
        iterator.left_child(dim);
        if iterator.valid_seq()
        {
            // assign left child
            let lc_index = map[iterator.index()];            
                                 
            array[active_index].set_has_left_child(true);
            array[active_index].down = lc_index as u32;
            array[active_index].set_has_child(true);
        }
        
        iterator.set_index(node_index);
        // get right child        
        iterator.right_child(dim);
        if iterator.valid_seq()
        {
            // assign right child
            let rc_index = map[iterator.index()];      
            array[active_index].set_has_right_child(true);
            // this potentially overwrites the down index if left child exists,
            // but that's ok. We just need one of them, or to know neither exist.
            array[active_index].down = rc_index as u32;   
            array[active_index].set_has_child(true);
        } 
        iterator.reset_to_left_level_zero(dim);
        // Handle left level zero
        if let Some(&index) = map.get(iterator.index()) {             
            // we know what the correct index is...
            let lzero = index as u32;
            // first let's find the leftmost node linked in our data structure.
            let mut left_index = seq as u32;
            while array[offset + left_index as usize].has_left()
            {
                left_index = array[offset + left_index as usize].left;
            }
            // If the leftmost node isn't the boundary, we need to  update our data structure
            // such that if it has boundaries, we set its left boundary neighbor.
            if left_index != lzero
            {
                array[offset + left_index as usize].left = lzero;
                array[offset + left_index as usize].set_has_left(true);
            }
        } 
        iterator.reset_to_right_level_zero(dim);
        if let Some(&index) = map.get(iterator.index())
        {
            // first let's find the rightmost node linked in our data structure.
            let rzero = index as u32;
            let mut right_index = seq as u32;
            while array[offset + right_index as usize].has_right()
            {
                right_index =  array[offset + right_index as usize].right;
            }
            // If the rightmost node isn't the boundary, we need to update our data structure
            // such that if it has boundaries, we set its right boundary neighbor.
            if right_index != rzero
            {
                array[offset + right_index as usize].right = rzero;
                array[offset + right_index as usize].set_has_right(true);
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
            if self.data.array[self.offset(dim) + index].has_left_child() 
            {            
                index = self.compute_left_child(dim, storage).unwrap() as usize;              
            }
        }
        else
        {
            while self.data.array[self.offset(dim) + index].has_parent()
            {            
                index = self.data.array[self.offset(dim) + index].up as usize;
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
                while self.data.array[self.offset(dim) +  right_index].has_right()
                {
                    right_index = self.data.array[self.offset(dim) +  right_index].right as usize;
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
                while self.data.array[self.offset(dim) +  left_index].has_left()
                {
                    left_index = self.data.array[self.offset(dim) +  left_index].left as usize;
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
        let offset = self.offset(dim);
        while self.data.array[offset + index].has_left()
        {
            index = self.data.array[offset + index].left as usize
        }
        if storage[index].index[dim] != 0
        {
            None
        }
        else {
            Some(index as u32)    
        }        
    }
    #[inline]
    pub fn compute_left_child(&self, dim: usize, storage: &SparseGridStorage<D>) -> Option<u32>
    {
        let original_node = &self.data.array[self.offset(dim) + self.index];
        let node_index = storage[self.index].index[dim];
        if original_node.has_child()
        {
            let index = original_node.down;
            let dim_index = storage[original_node.down as usize].index[dim];
            if dim_index == (2*node_index - 1)
            {
                Some(index)
            }
            else
            {               
                Some(self.data.array[self.offset(dim) + index as usize].left)
            }
        }
        else
        {
            None
        }
    }
    #[inline]
    fn compute_right_child(&self, dim: usize, storage: &SparseGridStorage<D>) -> Option<u32>
    {
        let original_node = &self.data.array[self.offset(dim) + self.index];
        let node_index = storage[self.index].index[dim];
        if original_node.has_child()
        {
            let index = original_node.down;
            let dim_index = storage[original_node.down as usize].index[dim];
            if dim_index == 2*node_index + 1
            {
                Some(index)
            }
            else
            {              
                Some(self.data.array[self.offset(dim) + index as usize].right)
            }
        }
        else
        {
            None
        }
    }

    /// Compute the index of the right boundary node (e.g. the right zero-level node) for the given dimension.
    #[inline]
    pub fn compute_rzero(&self, dim: usize, storage: &SparseGridStorage<D>) -> Option<u32> {
        let mut index = self.index;
        while self.data.array[self.offset(dim) + index].has_right()
        {
            index = self.data.array[self.offset(dim) + index].right as usize
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
    pub(crate) fn left_child(&mut self, dim: usize, storage: &SparseGridStorage<D>) -> bool
    {
        match self.data.array[self.offset(dim) + self.index].has_left_child()
        {
            true => 
            {                
                self.index = self.compute_left_child(dim, storage).unwrap() as usize;
                
                true
            },
            false => 
            {                
                assert_eq!(self.compute_left_child(dim, storage), None);
                false
            },
        }        
    }
    #[inline]
    pub(crate) fn right_child(&mut self, dim: usize, storage: &SparseGridStorage<D>) -> bool
    {
        match self.data.array[self.offset(dim) + self.index].has_right_child()
        {
            true => 
            {                
                self.index = self.compute_right_child(dim, storage).unwrap() as usize;
                true
            },
            false => 
            {        
                assert_eq!(self.compute_right_child(dim, storage), None);        
                false
            },
        }         
    }
   
    #[inline]
    pub(crate) fn is_leaf(&self, storage: &SparseGridStorage<D>) -> bool
    {
        storage[self.index].is_leaf()
    }
    #[inline]
    #[allow(unused)]
    pub(crate) fn has_left_leaf(&self, dim: usize) -> bool
    {
        self.data.array[self.offset(dim) + self.index].has_left_child()
    }
    #[inline]
    #[allow(unused)]
    pub(crate) fn has_right_leaf(&self, dim: usize) -> bool
    {
        self.data.array[self.offset(dim) + self.index].has_right_child()
    } 
    #[inline]
    #[allow(unused)]
    pub(crate) fn get_grid_depth(&mut self, dim: usize) -> usize
    {
        let mut depth = 1;
        let orig_index = self.index;
        loop
        {            
            if self.has_left_leaf(dim) || self.has_right_leaf(dim)
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