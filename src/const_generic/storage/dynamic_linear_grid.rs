use std::{hash::{Hash, Hasher}, ops::{Index, IndexMut}};
use rustc_hash::{FxHashMap, FxHasher};
use serde::{Deserialize, Serialize};
use serde_with::serde_as;
use bitfield_struct::bitfield;
use crate::iterators::{dynamic_grid_iterator::{DynamicHashMapGridIterator, GridIteratorT}};

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
    pub fn new(level: &[u8], is_leaf: bool) -> Self
    {
        let mut r = Self::default();        
        r.set_is_leaf(is_leaf);
        r.set_is_inner(!level.contains(&0));
        r
    }
    /// update `is_inner` flag...
    pub fn update_is_inner(&mut self, level: &[u8]) 
    {
        self.set_is_inner(!level.contains(&0));
    }
}
#[serde_as]
#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct GridPoint
{    
    pub level: Vec<u8>,    
    pub index: Vec<u32>,
    pub(crate) flags: GridPointFlags,
}
impl Hash for GridPoint
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.level.hash(state);
        self.index.hash(state);
    }
}
impl Default for GridPoint 
{
    fn default() -> Self {
        Self { level: vec![], index: vec![], flags: GridPointFlags(0) }
    }
}
impl PartialOrd for GridPoint
{    
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(std::cmp::Ord::cmp(self, other))
    }
}
impl Ord for GridPoint{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.index.cmp(&other.index).then(self.level.cmp(&other.level))
    }
}

impl PartialEq for GridPoint
{
    fn eq(&self, other: &Self) -> bool {
        self.level == other.level && self.index == other.index
    }
}
impl Eq for GridPoint{}

impl GridPoint
{
    pub fn new (level: &[u8], index: &[u32], is_leaf: bool) -> Self
    {
        let flags= GridPointFlags::new(&level, is_leaf);
        Self { level: level.to_vec(), index: index.to_vec(), flags }
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
    
    pub fn left_child(&self, dim: usize) -> GridPoint
    {
        let mut r = self.clone();
        if r.index[dim] == 0
        {
            r.index[dim] = u32::MAX;
            return r;
        }
        r.index[dim] = 2*self.index[dim] - 1;
        r.level[dim] += 1;
        r
    }
    pub fn right_child(&self, dim: usize) -> GridPoint
    {
        let mut r = self.clone();
        r.index[dim] = 2*self.index[dim] + 1;
        r.level[dim] += 1;
        r
    }
    ///
    /// returns an index with the top level in direction dim
    /// 
    pub fn root(&self, dim: usize) -> GridPoint
    {
        let mut r = self.clone();
        r.index[dim] = 1;
        r.level[dim] = 1;
        r
    }

    ///
    /// This only works for grids without boundaries.    
    /// 
    pub fn parent(&self, dim: usize) -> GridPoint
    {
        let mut r = self.clone();
        if self.level[dim] == 0
        {
            r.index[dim] = u32::MAX;
            return r;
        }
        r.index[dim] = (self.index[dim] >> 1) | 1;
        r.level[dim] -= 1;
        r
    }

    pub fn unit_coordinate(&self) -> Vec<f64>
    {
        let mut coor = vec![0.0; self.index.len()];
        #[allow(clippy::needless_range_loop)]
        for d in 0..self.index.len()
        {
            coor[d] = self.index[d] as f64 / (1 << self.level[d]) as f64;
        }
        coor
    }

    pub fn zero_index(num_inputs: usize) -> Self
    {
        Self{ level: vec![0; num_inputs], index: vec![0; num_inputs], flags: GridPointFlags(0) }
    }

}

impl From<GridPoint> for u64
{
    fn from(val: GridPoint) -> Self {
        let hasher = &mut FxHasher::default();
        val.hash(hasher);
        hasher.finish()
    }
}


#[derive(Clone, Serialize, Deserialize)]
pub struct BoundingBox
{
    pub lower: Vec<f64>,    
    pub upper: Vec<f64>
}

impl Default for BoundingBox
{
    #[inline]
    fn default() -> Self {
        Self { lower: vec![], upper: vec![] }
    }
}
impl BoundingBox
{
    #[inline]
    pub fn new(lower: &[f64], upper: &[f64]) -> Self
    {
        Self { lower: lower.to_vec(), upper: upper.to_vec() }
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
        for d in 0..self.lower.len()
        {
            volume *= self.width(d);
        }
        volume
    }
    #[inline]
    pub fn to_unit_coordinate(&self, point: &[f64]) -> Vec<f64>
    {
        let mut r = vec![0.0; point.len()];
        for i in 0..r.len()
        {
            r[i] = (point[i] - self.lower[i])/(self.upper[i] - self.lower[i]);
        }
        r
    }
    #[inline]
    pub fn to_real_coordinate(&self, point: &[f64]) -> Vec<f64>
    {
        let mut r = vec![0.0; point.len()];
        for i in 0..point.len()
        {
            r[i] = self.lower[i] + (self.upper[i] - self.lower[i]) * point[i];
        }
        r
    }
    #[inline]
    pub fn contains(&self, point: &[f64]) -> bool
    {
        #[allow(clippy::needless_range_loop)]
        for d in 0..point.len()
        {
            if self.lower[d] > point[d] || self.upper[d] < point[d]
            {
                return false;
            }
        }
        true
    }
}



impl SparseGridData
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
    pub fn contains(&self, point: &GridPoint) -> bool
    {
        let hasher = &mut FxHasher::default();
        point.hash(hasher);
        self.map.contains_key(&hasher.finish())
    }
    #[inline]
    pub fn get_mut(&mut self, point: &GridPoint) -> Option<&mut GridPoint>
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
    pub fn get(&mut self, point: &GridPoint) -> Option<&GridPoint>
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
    pub fn nodes(&self) -> &Vec<GridPoint>
    {
        &self.nodes
    }
    
    pub fn nodes_mut(&mut self) -> &mut Vec<GridPoint>
    {
        &mut self.nodes
    }
    
    #[inline]
    pub fn bounding_box(&self) -> &BoundingBox
    {
        &self.bounding_box
    }
    #[inline]
    pub fn bounding_box_mut(&mut self) -> &mut BoundingBox
    {
        &mut self.bounding_box
    }
    #[inline]
    pub fn index_of(&self, point: &GridPoint) -> Option<usize>
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
            new_list.push(self.nodes[i].clone());
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
            
            for dim in 0..self.num_inputs
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
    pub fn insert_point(&mut self, mut point: GridPoint) -> usize
    {
        // make sure our is_inner flag is up-to-date...
        point.flags.update_is_inner(&point.level);
        self.nodes.push(point.clone());        
        self.map.insert(point.into(), self.len() - 1);
        self.len() - 1
    }

    #[inline]
    pub fn update(&mut self, mut point: GridPoint, index: usize)
    {
        // make sure our is_inner flag is up-to-date...
        point.flags.update_is_inner(&point.level);
        let key: u64 = point.clone().into();
        self.map.insert(key, index);
        self.nodes[index] = point;
    }
    ///
    /// Return the real coordinates for each node...
    /// 
    pub fn points(&self) -> PointIterator<'_> {
        PointIterator::new(&self.nodes, self.bounding_box.clone())
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct SparseGridData
{
    pub(crate) bounding_box: BoundingBox,
    pub(crate) nodes: Vec<GridPoint>,
    pub(crate) adjacency_data: NodeAdjacencyData,    
    pub(crate) has_boundary: bool,
    #[serde(skip_serializing, skip_deserializing)]
    pub(crate) map: FxHashMap<u64, usize>,
    pub(crate) num_inputs: usize,
    pub(crate) num_outputs: usize,
}

impl Index<usize> for SparseGridData
{
    type Output = GridPoint;

    fn index(&self, index: usize) -> &Self::Output {
        &self.nodes[index]
    }
}
impl IndexMut<usize> for SparseGridData
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.nodes[index]
    }
}

impl Default for SparseGridData
{
    fn default() -> Self {
        Self { nodes: Default::default(), bounding_box: Default::default(), adjacency_data: NodeAdjacencyData::default(), has_boundary: false, map: FxHashMap::default(), num_inputs: 0, num_outputs: 0 }
    }
}

impl SparseGridData
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
        let mut array = vec![NodeAdjacency::default(); self.num_inputs*self.len()];
        if self.map.len() != self.len()
        {
            self.generate_map();
        }
        for dim in 0..self.num_inputs
        {
            for i in 0..self.len()
            {
                self.generate_adjacency_data_for_index(&mut array, i, dim);     
            }
        }
        self.adjacency_data.zero_index = self.index_of(&GridPoint::zero_index(self.num_inputs)).unwrap_or(usize::MAX);
        self.adjacency_data.data = array;
    }

    fn generate_adjacency_data_for_index(&mut self, array:&mut [NodeAdjacency], seq: usize, dim: usize)
    {

        let mut iterator = DynamicHashMapGridIterator::new(self);
        let offset =  dim * self.len();
        let active_index = offset + seq;
        if !array[active_index].inner.is_complete()
        {
            let node_index= &iterator.storage[seq]; 
            iterator.set_index(node_index.clone());
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
            iterator.set_index(node_index.clone());
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
            iterator.set_index(node_index.clone());
            iterator.up(dim);       
            if let Some(parent_index) = iterator.index()  // parent exists
            {
                // assign parent
                array[active_index].inner.set_up(parent_index as i64 - seq as i64);      
                array[active_index].inner.set_has_parent(true);
            }
            iterator.set_index(node_index.clone());
            // get left child
            iterator.left_child(dim);
            if let Some(lc_index) = iterator.index()
            {
                // assign left child
                array[active_index].inner.set_has_left_child(true);
                array[active_index].inner.set_down(lc_index as i64 - seq as i64);
                array[active_index].inner.set_has_child(true);
            }
            
            iterator.set_index(node_index.clone());
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

pub struct PointIterator<'a>
{
    data: &'a Vec<GridPoint>,
    bounding_box: BoundingBox,
    index: usize
}

impl<'a>  PointIterator<'a>
{
    pub fn new( data: &'a Vec<GridPoint>, bounding_box: BoundingBox) -> Self
    {
        Self { data, bounding_box, index: 0 }
    }
}

impl Iterator for PointIterator<'_>
{
    type Item=Vec<f64>;

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
