use std::{hash::{Hash, Hasher}, ops::{Index, IndexMut}};
use rustc_hash::{FxHashMap, FxHasher};
use serde::{Deserialize, Serialize};
use serde_with::serde_as;
use bitfield_struct::bitfield;
use crate::iterators::grid_iterator::{GridIteratorT, HashMapGridIterator};

use super::adjacency_data::{NodeAdjacency, NodeAdjacencyData};

#[bitfield(u8, new=false)]
#[derive(Serialize, Deserialize)]
pub struct GridPointFlags
{
    pub is_leaf: bool,
    pub is_inner: bool,
    #[bits(6)]
    pub _empty: u8
}
impl GridPointFlags
{
    pub fn new<const D: usize>(level: &[u8; D], is_leaf: bool) -> Self
    {
        let mut r = Self::default();        
        r.set_is_leaf(is_leaf);
        r.set_is_inner(!level.contains(&0));
        r
    }
    /// update `is_inner` flag...
    pub fn update_is_inner<const D: usize>(&mut self, level: &[u8; D]) 
    {
        self.set_is_inner(!level.contains(&0));
    }
}
#[serde_as]
#[derive(Copy, Clone, Serialize, Deserialize, Debug)]
pub struct GridPoint<const D: usize>
{
    #[serde_as(as = "[_; D]")]
    pub level: [u8; D],
    #[serde_as(as = "[_; D]")]
    pub index: [u32; D],
    pub(crate) flags: GridPointFlags,
}
impl<const D: usize> Hash for GridPoint<D>
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.level.hash(state);
        self.index.hash(state);
    }
}
impl<const D: usize> Default for GridPoint<D> 
{
    fn default() -> Self {
        Self { level: [1; D], index: [1; D], flags: GridPointFlags(0) }
    }
}
impl<const D: usize> PartialOrd for GridPoint<D>
{    
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(std::cmp::Ord::cmp(self, other))
    }
}
impl<const D:usize> Ord for GridPoint<D>{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.index.cmp(&other.index).then(self.level.cmp(&other.level))
    }
}

impl<const D: usize> PartialEq for GridPoint<D>
{
    fn eq(&self, other: &Self) -> bool {
        self.level == other.level && self.index == other.index
    }
}
impl<const D: usize> Eq for GridPoint<D>{}

impl<const D: usize> GridPoint<D>
{
    pub fn new (level: [u8; D], index: [u32; D], is_leaf: bool) -> Self
    {
        let flags= GridPointFlags::new(&level, is_leaf);
        Self { level, index, flags }
    }   
    pub fn is_leaf(&self) -> bool
    {
        self.flags.is_leaf()
    }
    pub fn set_is_leaf(&mut self, is_leaf: bool)
    {
        self.flags.set_is_leaf(is_leaf);
    }

    pub fn set_is_inner(&mut self, is_inner: bool)
    {        
        self.flags.set_is_inner(is_inner);
    }
    /// 
    /// This is an inner point if no indices are zero...
    /// 
    pub fn is_inner_point(&self) -> bool
    {        
        self.flags.is_inner()
    }   
    pub fn level_sum(&self) -> u8
    {
        self.level.iter().sum()
    }
    #[inline]
    pub fn level_max(&self) -> u8
    {
        *self.level.iter().max().unwrap()
    }
    pub fn level_min(&self) -> u8
    {
        *self.level.iter().min().unwrap()
    }
    
    pub fn left_child(&self, dim: usize) -> GridPoint<D>
    {
        let mut r = *self;
        if r.index[dim] == 0
        {
            r.index[dim] = u32::MAX;
            return r;
        }
        r.index[dim] = 2*self.index[dim] - 1;
        r.level[dim] += 1;
        r
    }
    pub fn right_child(&self, dim: usize) -> GridPoint<D>
    {
        let mut r = *self;
        r.index[dim] = 2*self.index[dim] + 1;
        r.level[dim] += 1;
        r
    }
    ///
    /// returns an index with the top level in direction dim
    /// 
    pub fn root(&self, dim: usize) -> GridPoint<D>
    {
        let mut r = *self;
        r.index[dim] = 1;
        r.level[dim] = 1;
        r
    }

    ///
    /// This only works for grids without boundaries.    
    /// 
    pub fn parent(&self, dim: usize) -> GridPoint<D>
    {
        let mut r = *self;
        if self.level[dim] == 0
        {
            r.index[dim] = u32::MAX;
            return r;
        }
        r.index[dim] = (self.index[dim] >> 1) | 1;
        r.level[dim] -= 1;
        r
    }

    pub fn unit_coordinate(&self) -> [f64; D]
    {
        let mut coor = [0.0; D];
        #[allow(clippy::needless_range_loop)]
        for d in 0..D
        {
            coor[d] = self.index[d] as f64 / (1 << self.level[d]) as f64;
        }
        coor
    }

    pub const fn zero_index() -> Self
    {
        Self{ level: [0; D], index: [0; D], flags: GridPointFlags(0) }
    }

}

impl<const D: usize> From<GridPoint<D>> for u64
{
    fn from(val: GridPoint<D>) -> Self {
        let hasher = &mut FxHasher::default();
        val.hash(hasher);
        hasher.finish()
    }
}


#[serde_as]
#[derive(Copy, Clone, Serialize, Deserialize)]
pub struct BoundingBox<const D: usize>
{
    #[serde_as(as = "[_; D]")]
    pub lower: [f64; D],
    #[serde_as(as = "[_; D]")]
    pub upper: [f64; D],
}

impl<const D: usize> Default for BoundingBox<D>
{
    #[inline]
    fn default() -> Self {
        Self { lower: [0.0; D], upper: [1.0; D] }
    }
}
impl<const D: usize> BoundingBox<D>
{
    #[inline]
    pub fn new(lower: [f64; D], upper: [f64; D]) -> Self
    {
        Self { lower, upper }
    }
    #[inline]
    pub fn width(&self, dim: usize) -> f64
    {
        self.upper[dim] - self.lower[dim]
    }

    ///
    /// Volume of hypercube (width(dim1)*...*width(dim_n))
    /// 
    #[inline]
    pub fn volume(&self) -> f64
    {
        let mut volume = 1.0;
        for d in 0..D
        {
            volume *= self.width(d);
        }
        volume
    }
    #[inline]
    pub fn to_unit_coordinate(&self, point: &[f64; D]) -> [f64; D]
    {
        let mut r = [0.0;D];
        for i in 0..D
        {
            r[i] = (point[i] - self.lower[i])/(self.upper[i] - self.lower[i]);
        }
        r
    }
    #[inline]
    pub fn to_real_coordinate(&self, point: &[f64; D]) -> [f64; D]
    {
        let mut r = [0.0;D];
        for i in 0..D
        {
            r[i] = self.lower[i] + (self.upper[i] - self.lower[i]) * point[i];
        }
        r
    }
    #[inline]
    pub fn contains(&self, point: &[f64; D]) -> bool
    {
        #[allow(clippy::needless_range_loop)]
        for d in 0..D
        {
            if self.lower[d] > point[d] || self.upper[d] < point[d]
            {
                return false;
            }
        }
        true
    }
}



impl<const D: usize> SparseGridData<D>
{
    #[inline(always)]
    pub fn len(&self) -> usize
    {
        self.nodes.len()
    }
    #[inline(always)]
    pub fn is_empty(&self) -> bool
    {
        self.nodes.is_empty()
    }
    #[inline(always)]
    pub fn has_boundary(&self) -> bool
    {
        self.has_boundary
    }

    #[inline]
    pub fn contains(&self, point: &GridPoint<D>) -> bool
    {
        let hasher = &mut FxHasher::default();
        point.hash(hasher);
        self.map.contains_key(&hasher.finish())
    }
    #[inline]
    pub fn get_mut(&mut self, point: &GridPoint<D>) -> Option<&mut GridPoint<D>>
    {
        let hasher = &mut FxHasher::default();
        point.hash(hasher);
        match self.map.get(&hasher.finish())
        {
            Some(&seq) => Some(&mut self.nodes[seq]),
            None => None
        }
    }

    #[inline]
    pub fn get(&mut self, point: &GridPoint<D>) -> Option<&GridPoint<D>>
    {
        let hasher = &mut FxHasher::default();
        point.hash(hasher);
        match self.map.get(&hasher.finish())
        {
            Some(&seq) => Some(&self.nodes[seq]),
            None => None
        }
    }

    #[inline(always)]
    pub fn nodes(&self) -> &Vec<GridPoint<D>>
    {
        &self.nodes
    }
    
    pub fn nodes_mut(&mut self) -> &mut Vec<GridPoint<D>>
    {
        &mut self.nodes
    }
    
    #[inline]
    pub fn bounding_box(&self) -> &BoundingBox<D>
    {
        &self.bounding_box
    }
    #[inline]
    pub fn bounding_box_mut(&mut self) -> &mut BoundingBox<D>
    {
        &mut self.bounding_box
    }
    #[inline]
    pub fn index_of(&self, point: &GridPoint<D>) -> Option<usize>
    {
        let hasher = &mut FxHasher::default();
        point.hash(hasher);
        self.map.get(&hasher.finish()).copied()
    }

    pub fn remove(&mut self, points_to_keep: &indexmap::IndexSet<usize>)
    {
        let mut new_list = Vec::with_capacity(points_to_keep.len());
        for &i in points_to_keep.iter()
        {            
            new_list.push(self.nodes[i]);
        }
        self.nodes = new_list;
        // TODO: This doesn't work for some reason...
        //self.map.retain(|_, v| points_to_keep.contains(v));     
        self.generate_map();
        self.update_leaves();
    }

    fn update_leaves(&mut self)
    {
        #[allow(clippy::needless_range_loop)]
        for i in 0..self.len()
        {
            let point = &mut self.nodes[i];
            let mut is_leaf = true;
            
            for dim in 0..D
            {                
                if point.level[dim] > 0
                {   
                    // Check if this point has any children. If not it is a leaf.
                    if self.map.contains_key(&point.left_child(dim).into()) ||
                        self.map.contains_key(&point.right_child(dim).into()) 
                    {
                        is_leaf = false;
                        break; 
                    }                
                }
                else
                {                    
                    // don't remove boundary nodes that are used by other nodes...
                    if self.map.contains_key(&point.root(dim).into()) 
                    {
                        is_leaf = false;
                        break;
                    }
                }
            }                
            point.flags.set_is_leaf(is_leaf);            
        }        
    }

    ///
    /// Inserts a point
    /// 
    #[inline]
    pub fn insert_point(&mut self, mut point: GridPoint<D>) -> usize
    {
        // make sure our is_inner flag is up-to-date...
        point.flags.update_is_inner(&point.level);
        self.nodes.push(point);        
        self.map.insert(point.into(), self.len() - 1);
        self.len() - 1
    }

    #[inline]
    pub fn update(&mut self, mut point: GridPoint<D>, index: usize)
    {
        // make sure our is_inner flag is up-to-date...
        point.flags.update_is_inner(&point.level);
        let key: u64 = point.into();
        self.map.insert(key, index);
        self.nodes[index] = point;
    }
    ///
    /// Return the real coordinates for each node...
    /// 
    pub fn points(&self) -> PointIterator<'_, D> {
        PointIterator::new(&self.nodes, self.bounding_box)
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct SparseGridData<const D: usize>
{
    pub(crate) bounding_box: BoundingBox<D>,
    pub(crate) nodes: Vec<GridPoint<D>>,
    pub(crate) adjacency_data: NodeAdjacencyData,    
    pub(crate) has_boundary: bool,
    #[serde(skip_serializing, skip_deserializing)]
    pub(crate) map: FxHashMap<u64, usize>
}

impl<const D: usize> Index<usize> for SparseGridData<D>
{
    type Output = GridPoint<D>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.nodes[index]
    }
}
impl<const D: usize> IndexMut<usize> for SparseGridData<D>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.nodes[index]
    }
}

impl<const D: usize> Default for SparseGridData<D>
{
    fn default() -> Self {
        Self { nodes: Default::default(), bounding_box: Default::default(), adjacency_data: NodeAdjacencyData::default(), has_boundary: false, map: FxHashMap::default() }
    }
}

impl<const D: usize> SparseGridData<D>
{
    pub fn generate_map(&mut self)
    {
        
        let mut map = FxHashMap::default();
        for (i, node) in self.nodes.iter().enumerate()
        {
            let mut hasher = FxHasher::default();
            node.hash(&mut hasher);            
            map.insert(hasher.finish(), i);
        }
        self.map = map;
    }

    pub fn generate_adjacency_data(&mut self)
    {
        let mut array = vec![NodeAdjacency::default(); D*self.len()];
        if self.map.len() != self.len()
        {
            self.generate_map();
        }
        for dim in 0..D
        {
            for i in 0..self.len()
            {
                self.generate_adjacency_data_for_index(&mut array, i, dim);     
            }
        }
        self.adjacency_data.zero_index = self.index_of(&GridPoint::zero_index()).unwrap_or(usize::MAX);
        self.adjacency_data.data = array;
    }

    fn generate_adjacency_data_for_index(&mut self, array:&mut [NodeAdjacency], seq: usize, dim: usize)
    {

        let mut iterator = HashMapGridIterator::new(self);
        let offset =  dim * self.len();
        let active_index = offset + seq;
        if !array[active_index].inner.is_complete()
        {
            let node_index= iterator.storage[seq]; 
            iterator.set_index(node_index);
            iterator.step_left(dim);
            if let Some(left_index) = iterator.index()
            {                
                // assign left node            
                array[active_index].inner.set_left(left_index as i64 - seq as i64);
                array[active_index].inner.set_has_left(true);
                // assign right node for left of current node...
                array[offset + left_index].inner.set_right(seq as i64 - left_index as i64);
                array[offset + left_index].inner.set_has_right(true);
            }
            iterator.set_index(node_index);
            iterator.step_right(dim);     
            if let Some(right_index) = iterator.index()
            {
                // assign right node
                array[active_index].inner.set_right(right_index as i64 - seq as i64);
                array[active_index].inner.set_has_right(true);
                // assign left node for right of current node...
                array[offset + right_index].inner.set_left(seq as i64 - right_index as i64);
                array[offset + right_index].inner.set_has_left(true);
            }        
            iterator.set_index(node_index);
            iterator.up(dim);       
            if let Some(parent_index) = iterator.index()  // parent exists
            {
                // assign parent
                array[active_index].inner.set_up(parent_index as i64 - seq as i64);      
                array[active_index].inner.set_has_parent(true);
            }
            iterator.set_index(node_index);
            // get left child
            iterator.left_child(dim);
            if let Some(lc_index) = iterator.index()
            {
                // assign left child
                array[active_index].inner.set_has_left_child(true);
                array[active_index].inner.set_down(lc_index as i64 - seq as i64);
                array[active_index].inner.set_has_child(true);
            }
            
            iterator.set_index(node_index);
            // get right child        
            iterator.right_child(dim);
            if let Some(rc_index) = iterator.index()
            {
                // assign right child
                array[active_index].inner.set_has_right_child(true);
                // this potentially overwrites the down index if left child exists,
                // but that's ok. We just need one of them, or to know neither exist.
                array[active_index].inner.set_down(rc_index as i64 - seq as i64);
                array[active_index].inner.set_has_child(true);
            } 
            iterator.reset_to_left_level_zero(dim);            
            // Handle left level zero
            if let Some(lzero) = iterator.index() {                             
                // we know what the correct index is...
                // first let's find the leftmost node linked in our data structure.
                let mut left_index = seq as u32;
                while array[offset + left_index as usize].inner.has_left()
                {
                    left_index = (left_index as i64 + array[offset + left_index as usize].inner.left()) as u32;
                }
                // If the leftmost node isn't the boundary, we need to update our data structure
                // such that if it has boundaries, we set its left boundary neighbor.
                if left_index != lzero as u32
                {
                    array[offset + left_index as usize].inner.set_left(lzero as i64 - left_index as i64);
                    array[offset + left_index as usize].inner.set_has_left(true);
                }
                array[active_index].left_zero = lzero as u32;
            } 
            iterator.reset_to_right_level_zero(dim);
            if let Some(rzero) = iterator.index()
            {
                // first let's find the rightmost node linked in our data structure.
                let mut right_index = seq as u32;
                while array[offset + right_index as usize].inner.has_right()
                {
                    right_index = (right_index as i64 + array[offset + right_index as usize].inner.right()) as u32;
                }
                // If the rightmost node isn't the boundary, we need to update our data structure
                // such that if it has boundaries, we set its right boundary neighbor.
                if right_index != rzero as u32
                {
                    array[offset + right_index as usize].inner.set_right(rzero as i64 - right_index as i64);
                    array[offset + right_index as usize].inner.set_has_right(true);
                }
            }       
        }
    }
}

pub struct PointIterator<'a, const D: usize>
{
    data: &'a Vec<GridPoint<D>>,
    bounding_box: BoundingBox<D>,
    index: usize
}

impl<'a, const D: usize>  PointIterator<'a, D>
{
    pub fn new( data: &'a Vec<GridPoint<D>>, bounding_box: BoundingBox<D>) -> Self
    {
        Self { data, bounding_box, index: 0 }
    }
}

impl<const D: usize> Iterator for PointIterator<'_, D>
{
    type Item=[f64; D];

    fn next(&mut self) -> Option<Self::Item> {

        if self.index < self.data.len()
        {
            let mut point = self.data[self.index].unit_coordinate();
            point = self.bounding_box.to_real_coordinate(&point);
            self.index += 1;
            Some(point)
        }
        else
        {
            None
        }
    }
}
