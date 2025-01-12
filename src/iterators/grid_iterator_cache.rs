use std::collections::HashMap;

use rustc_hash::FxBuildHasher;
use serde::{Deserialize, Serialize};

use crate::storage::linear_grid::{GridPoint, SparseGridData, SparseGridStorage};

use super::grid_iterator::{GridIterator, GridIteratorT};

#[derive(Serialize, Deserialize)]
#[derive(Copy, Clone, Default)]
pub(crate) struct NodeProperties {
    data: u128, // Stores all fields: left, right, down, up, and flags
}

impl NodeProperties {
    /// Number of bits for each field
    const FIELD_BITS: u128 = 30;
    const FIELD_BITS_LIMIT: u128 = Self::FIELD_BITS - 1;
    /// Mask for extracting the 30 bits (including sign bit)
    const FIELD_MASK: u128 = (1 << Self::FIELD_BITS) - 1;
    /// Offset for each field
    const LEFT_OFFSET: u128 = 0;
    const RIGHT_OFFSET: u128 = Self::FIELD_BITS;
    const DOWN_OFFSET: u128 = 2 * Self::FIELD_BITS;
    const UP_OFFSET: u128 = 3 * Self::FIELD_BITS;
    const FLAG_OFFSET: u128 = 120; // Flags occupy the highest 8 bits

    fn pack_signed(value: i32) -> u32 {
        (value as u32) & Self::FIELD_MASK as u32
    }
    #[inline]
    fn unpack_signed(value: u128) -> i32 {
        let result = (value & Self::FIELD_MASK) as i32;
        // If the sign bit (Self::FIELD_BITS_LIMITth bit) is set, extend the sign
        if result & (1 << Self::FIELD_BITS_LIMIT) != 0 {
            result | !((1 << Self::FIELD_BITS) - 1) // Extend the sign to a full i32
        } else {
            result
        }
    }
    #[inline]
    pub fn left(&self) -> i32 {
        Self::unpack_signed(self.data >> Self::LEFT_OFFSET)
    }

    pub fn set_left(&mut self, left: i32) {
        assert!((-(1 << Self::FIELD_BITS_LIMIT)..(1 << Self::FIELD_BITS_LIMIT)).contains(&left), "left out of range");
        self.data = (self.data & !(Self::FIELD_MASK << Self::LEFT_OFFSET))
            | ((Self::pack_signed(left) as u128) << Self::LEFT_OFFSET);
    }
    #[inline]
    pub fn right(&self) -> i32 {
        Self::unpack_signed(self.data >> Self::RIGHT_OFFSET)
    }

    pub fn set_right(&mut self, right: i32) {
        assert!((-(1 << Self::FIELD_BITS_LIMIT)..(1 << Self::FIELD_BITS_LIMIT)).contains(&right), "right out of range");
        self.data = (self.data & !(Self::FIELD_MASK << Self::RIGHT_OFFSET))
            | ((Self::pack_signed(right) as u128) << Self::RIGHT_OFFSET);
    }
    #[inline]
    pub fn down(&self) -> i32 {
        Self::unpack_signed(self.data >> Self::DOWN_OFFSET)
    }

    pub fn set_down(&mut self, down: i32) {
        assert!((-(1 << Self::FIELD_BITS_LIMIT)..(1 << Self::FIELD_BITS_LIMIT)).contains(&down), "down out of range");
        self.data = (self.data & !(Self::FIELD_MASK << Self::DOWN_OFFSET))
            | ((Self::pack_signed(down) as u128) << Self::DOWN_OFFSET);
    }
    #[inline]
    pub fn up(&self) -> i32 {
        Self::unpack_signed(self.data >> Self::UP_OFFSET)
    }

    pub fn set_up(&mut self, up: i32) {
        assert!((-(1 << Self::FIELD_BITS_LIMIT)..(1 << Self::FIELD_BITS_LIMIT)).contains(&up), "up out of range");
        self.data = (self.data & !(Self::FIELD_MASK << Self::UP_OFFSET))
            | ((Self::pack_signed(up) as u128) << Self::UP_OFFSET);
    }
    #[inline(always)]
    pub fn flags(&self) -> u8 {
        ((self.data >> Self::FLAG_OFFSET) & 0xFF) as u8
    }

    pub fn set_flags(&mut self, flags: u8) {
        self.data = (self.data & !(0xFF << Self::FLAG_OFFSET)) | ((flags as u128) << Self::FLAG_OFFSET);
    }
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
        let mut flags = self.flags();
        if value {
            flags |= 1 << INDEX;
        } else {
            flags &= !(1 << INDEX);
        }
        self.set_flags(flags);
    }

    #[inline]
    fn get_flag<const INDEX: usize>(&self) -> bool {
        (self.flags() & (1 << INDEX)) != 0
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
#[derive(Default,Serialize, Deserialize, Clone)]
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
        let zero_idx = iterator.node();
        let zero_index = *storage.map.get(zero_idx).unwrap_or(&0);
        GridIteratorData { zero_level_node: zero_index,  array, storage_len: storage.len()}
    }
}

///
/// This iterator uses array access to retrieve neighbors more quickly (~10x speedup relative to hash queries)
/// at the expense of pretty significant memory overhead. Eventually the data behind building the iterator
/// will be an optional format
/// 
pub(crate) struct GridIteratorWithCache<'a, 'b, const D: usize>
{
    pub(crate) index: usize,    
    data: &'b GridIteratorData<D>,
    storage: &'a SparseGridData<D>,
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
            let left_index = map[iterator.node()];
            // assign left node            
            array[active_index].set_left((left_index as i64 - seq as i64).try_into().unwrap());
            array[active_index].set_has_left(true);
            // assign right node for left of current node...
            array[offset + left_index].set_right((seq as i64 - left_index as i64).try_into().unwrap());
            array[offset + left_index].set_has_right(true);
        }
        iterator.set_index(node_index);
        iterator.step_right(dim);     
        if iterator.valid_seq()
        {
            // assign right node
            let right_index = map[iterator.node()];
            array[active_index].set_right((right_index as i64 - seq as i64).try_into().unwrap());
            array[active_index].set_has_right(true);
            // assign left node for right of current node...
            array[offset + right_index].set_left((seq as i64 - right_index as i64).try_into().unwrap());
            array[offset + right_index].set_has_left(true);
        }        
        iterator.set_index(node_index);
        iterator.up(dim);       
        if iterator.valid_seq()  // parent exists
        {
            let parent_index = map[iterator.node()];
            // assign parent
            array[active_index].set_up((parent_index as i64 - seq as i64).try_into().unwrap());      
            array[active_index].set_has_parent(true);
        }
        iterator.set_index(node_index);
        // get left child
        iterator.left_child(dim);
        if iterator.valid_seq()
        {
            // assign left child
            let lc_index = map[iterator.node()];            
                                 
            array[active_index].set_has_left_child(true);
            array[active_index].set_down((lc_index as i64 - seq as i64).try_into().unwrap());
            array[active_index].set_has_child(true);
        }
        
        iterator.set_index(node_index);
        // get right child        
        iterator.right_child(dim);
        if iterator.valid_seq()
        {
            // assign right child
            let rc_index = map[iterator.node()];      
            array[active_index].set_has_right_child(true);
            // this potentially overwrites the down index if left child exists,
            // but that's ok. We just need one of them, or to know neither exist.
            array[active_index].set_down((rc_index as i64 - seq as i64).try_into().unwrap());
            array[active_index].set_has_child(true);
        } 
        iterator.reset_to_left_level_zero(dim);
        // Handle left level zero
        if let Some(&index) = map.get(iterator.node()) {             
            // we know what the correct index is...
            let lzero = index as u32;
            // first let's find the leftmost node linked in our data structure.
            let mut left_index = seq as u32;
            while array[offset + left_index as usize].has_left()
            {
                left_index = (left_index as i64 + array[offset + left_index as usize].left() as i64) as u32;
            }
            // If the leftmost node isn't the boundary, we need to update our data structure
            // such that if it has boundaries, we set its left boundary neighbor.
            if left_index != lzero
            {
                array[offset + left_index as usize].set_left((lzero as i64 - left_index as i64).try_into().unwrap());
                array[offset + left_index as usize].set_has_left(true);
            }
        } 
        iterator.reset_to_right_level_zero(dim);
        if let Some(&index) = map.get(iterator.node())
        {
            // first let's find the rightmost node linked in our data structure.
            let rzero = index as u32;
            let mut right_index = seq as u32;
            while array[offset + right_index as usize].has_right()
            {
                right_index = (right_index as i64 + array[offset + right_index as usize].right() as i64) as u32;
            }
            // If the rightmost node isn't the boundary, we need to update our data structure
            // such that if it has boundaries, we set its right boundary neighbor.
            if right_index != rzero
            {
                array[offset + right_index as usize].set_right((rzero as i64 - right_index as i64).try_into().unwrap());
                array[offset + right_index as usize].set_has_right(true);
            }
        }       
    }

    
}

impl<'a, 'b, const D: usize> GridIteratorWithCache<'a, 'b,  D>
{

    /// This computes the left level one (e.g. index=1 and level=1 for the given dimension).
    #[inline]
    pub fn compute_level_one(&self, dim: usize) -> Option<u32>
    {
        let mut index = self.index;
        let level = self.storage[index].level[dim];
        // First we determine which way we iterate up or down the level hierarchy...
        if level == 0
        {
            if self.data.array[self.offset(dim) + index].has_left_child() 
            {            
                index = self.compute_left_child(dim) as usize;              
            }
        }
        else
        {
            while self.data.array[self.offset(dim) + index].has_parent()
            {            
                index = (index as i64 + self.data.array[self.offset(dim) + index].up() as i64) as usize;
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
                while self.data.array[self.offset(dim) + right_index].has_right()
                {
                    right_index = (right_index as i64 + self.data.array[self.offset(dim) + right_index].right() as i64) as usize;
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
                while self.data.array[self.offset(dim) + left_index].has_left()
                {
                    left_index = (left_index as i64 + self.data.array[self.offset(dim) + left_index].left() as i64) as usize;
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
        let mut index = self.index;
        let offset = self.offset(dim);
        while self.data.array[offset + index].has_left()
        {
            index = (index as i64 + self.data.array[offset + index].left() as i64) as usize
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
        let original_node = &self.data.array[self.offset(dim) + self.index];
        let node_index = self.storage[self.index].index[dim];
        if original_node.has_child()
        {
            let index = self.index as i64 + original_node.down() as i64;
            let dim_index = self.storage[index as usize].index[dim];
            if dim_index == (2*node_index - 1)
            {
                index as u32
            }
            else
            {               
                (index + self.data.array[self.offset(dim) + index as usize].left() as i64) as u32
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
        let original_node = &self.data.array[self.offset(dim) + self.index];
        let node_index = self.storage[self.index].index[dim];
        if original_node.has_child()
        {
            let index = self.index as i64 + original_node.down() as i64;
            let dim_index = self.storage[index as usize].index[dim];
            if dim_index == 2*node_index + 1
            {
                index as u32
            }
            else
            {   
                // We check that the right child exists before this is called, so this should never panic.           
                (index  + self.data.array[ self.offset(dim) + index as usize].right() as i64) as u32
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
        let mut index = self.index;
        let offset = self.offset(dim);
        while self.data.array[ offset + index].has_right()
        {
            index = (index as i64 + self.data.array[offset + index].right() as i64) as usize            
        }        
        if self.storage[index].index[dim] != 1
        {
            None
        }
        else {
            Some(index as u32)    
        }
    }
    pub(crate) fn new(data: &'b GridIteratorData<D>, storage: &'a SparseGridData<D>) -> Self
    {
        Self { index: 0, data, storage}
    }
    #[inline(always)]
    fn offset(&self, dim: usize) -> usize
    {
        dim * self.data.storage_len
    }
    

    #[inline]
    #[allow(unused)]
    pub(crate) fn is_valid(&self) -> bool
    {
        self.index < self.data.storage_len
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
  
}

impl<const D: usize> GridIteratorT<D> for GridIteratorWithCache<'_, '_, D>
{
    #[inline]
    fn reset_to_level_zero(&mut self) -> bool
    {   
        self.index = self.data.zero_level_node;
        true
    }
    #[inline]
    fn reset_to_left_level_zero(&mut self, dim: usize) -> bool
    {
        if let Some(index) = self.compute_lzero(dim)
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
    fn reset_to_right_level_zero(&mut self, dim: usize) -> bool
    {
        if let Some(index) = self.compute_rzero(dim)
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
    fn reset_to_level_one(&mut self, dim: usize) -> bool
    {
        match self.compute_level_one(dim)
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
    fn left_child(&mut self, dim: usize) -> bool
    {
        match self.data.array[self.offset(dim) + self.index].has_left_child()
        {
            true => 
            {                
                self.index = self.compute_left_child(dim) as usize;
                
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
        match self.data.array[self.offset(dim) + self.index].has_right_child()
        {
            true => 
            {                
                self.index = self.compute_right_child(dim) as usize;
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
        self.storage[self.index].is_leaf()
    }
    
    fn index(&self) -> usize {
        self.index
    }
    
    fn up(&mut self, dim: usize) -> bool
    {
        if self.data.array[self.offset(dim) + self.index].has_parent()
        {
            self.index = (self.index as i64 + self.data.array[self.offset(dim) + self.index].up() as i64) as usize;
            true
        }
        else
        {
            false
        }
    }
    
    fn is_inner_point(&self) -> bool {
        self.storage[self.index].is_inner_point()
    }
    
    fn node(&self) -> &GridPoint<D> {
        &self.storage[self.index]
    }
}