use std::hash::{Hash, Hasher};
use rustc_hash::{FxHashMap, FxHasher};
use serde::{Deserialize, Serialize};
use serde_with::serde_as;
use bitfield_struct::bitfield;
use crate::{dynamic::iterators::dynamic_grid_iterator::{DynamicHashMapGridIterator, GridIteratorT}};

use crate::adjacency_data::{NodeAdjacency, NodeAdjacencyData};

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
impl From<&GridPoint> for u64
{
    fn from(val: &GridPoint) -> Self {
        let hasher = &mut FxHasher::default();
        val.hash(hasher);
        hasher.finish()
    }
}

pub struct GridPointRef<'a> {
    pub(crate) index: &'a [u32],
    pub(crate) level: &'a [u8],
    pub(crate) flags: &'a GridPointFlags
}
impl GridPointRef<'_>
{
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
}

impl PartialEq for GridPointRef<'_>
{
    fn eq(&self, other: &Self) -> bool {
        self.level == other.level && self.index == other.index
    }
}
impl Eq for GridPointRef<'_>{}

impl PartialOrd for GridPointRef<'_>
{    
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(std::cmp::Ord::cmp(self, other))
    }
}

impl Ord for GridPointRef<'_>{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.index.cmp(&other.index).then(self.level.cmp(&other.level))
    }
}

impl<'a> From<(&'a [u32], &'a [u8], &'a GridPointFlags)> for GridPoint
{
    fn from((index, level, flags): (&[u32], &[u8], &GridPointFlags)) -> Self {
        Self { index: index.to_owned(), level: level.to_owned(), flags: *flags }
    }
}


impl<'a> From<(&'a [u32], &'a [u8], &'a GridPointFlags)> for GridPointRef<'a>
{
    fn from((index, level, flags): (&'a [u32], &'a [u8], &'a GridPointFlags)) -> Self {
        Self { index, level, flags }
    }
}

impl<'a> std::hash::Hash for GridPointRef<'a>
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.level.hash(state);
        self.index.hash(state);
    }
}
impl From<GridPointRef<'_>> for u64
{
    fn from(val: GridPointRef<'_>) -> Self {
        let hasher = &mut FxHasher::default();
        val.hash(hasher);
        hasher.finish()
    }
}

impl From<GridPointRef<'_>> for GridPoint
{
    fn from(value: GridPointRef<'_>) -> Self {
        GridPoint { level: value.level.to_owned(), index: value.index.to_owned(), flags: *value.flags }
    }
}

pub struct GridPointMutRef<'a> {
    #[allow(unused)]
    pub(crate) index: &'a mut [u32],
    #[allow(unused)]
    pub(crate) level: &'a mut [u8],
    pub(crate) flags: &'a mut GridPointFlags
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
    pub fn with_dim(num_inputs: usize) -> Self
    {
        Self { lower: vec![0.0; num_inputs], upper: vec![1.0; num_inputs] }
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
    pub fn to_unit_coordinate(&self, point: &[f64]) -> [f64; 128]
    {
        let mut r = [0.0; 128];
        for i in 0..point.len()
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
    pub fn to_real_coordinate_in_place(&self, point: &mut [f64])
    {        
        for i in 0..point.len()
        {
            point[i] = self.lower[i] + (self.upper[i] - self.lower[i]) * point[i];
        }
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


#[derive(Serialize, Deserialize, Clone)]
pub struct SparseGridData
{
    pub bounding_box: BoundingBox,
    pub(crate) index: Vec<u32>,
    pub(crate) level: Vec<u8>,
    pub(crate) flags: Vec<GridPointFlags>,
    pub(crate) num_inputs: usize,
    pub(crate) map: FxHashMap<u64, usize>,
    pub(crate) adjacency_data: NodeAdjacencyData,
    pub(crate) has_boundary: bool,
    pub(crate) num_outputs: usize,
}

impl SparseGridData
{
    pub fn new(num_inputs: usize, num_outputs: usize) -> Self
    {
        Self { bounding_box: BoundingBox::with_dim(num_inputs), index: Vec::new(), level: Vec::new(), flags: Vec::new(), num_inputs, num_outputs, map: FxHashMap::default(), adjacency_data: NodeAdjacencyData::default(), has_boundary: false }
    }
    #[inline]
    pub fn num_inputs(&self) -> usize
    {
        self.num_inputs
    }
    #[inline]
    pub fn point(&self,seq: usize) -> GridPoint
    {
        GridPoint::from((&self.index[seq*self.num_inputs..(seq+1)*self.num_inputs], &self.level[seq*self.num_inputs..(seq+1)*self.num_inputs], &self.flags[seq]))
    }

    #[inline]
    pub fn is_empty(&self) -> bool
    {
        self.index.is_empty()
    }

    #[inline]
    pub fn len(&self) -> usize
    {
        self.index.len() / self.num_inputs
    }

    #[inline(always)]
    pub fn has_boundary(&self) -> bool
    {
        self.has_boundary
    }
    
    #[inline]
    pub fn index(&self, seq: usize, dim: usize) -> u32
    {
        self.index[self.num_inputs*seq + dim]
    }    

    #[inline]
    pub fn set_index(&mut self, seq: usize, dim: usize, value: u32)
    {
        self.index[self.num_inputs*seq + dim] = value;
    }

    #[inline(always)]
    pub fn level(&self, seq: usize, dim: usize) -> u8
    {
        self.level[self.num_inputs*seq + dim]
    }

    #[inline]
    pub fn get_mut(&mut self, point: &GridPoint) -> Option<GridPointMutRef<'_>>
    {
        if let Some(&seq) = self.map.get(&point.into())
        {
            Some(GridPointMutRef { index: &mut self.index[seq*self.num_inputs..(seq+1)*self.num_inputs], level: &mut self.level[seq*self.num_inputs..(seq+1)*self.num_inputs], flags: &mut self.flags[seq] })
        }
        else {
            None
        }
    }


    #[inline]
    pub fn set_level(&mut self, seq: usize, dim: usize, value: u8)
    {
        self.level[self.num_inputs*seq + dim] = value;
    }

    #[inline]
    pub fn is_leaf(&self, seq: usize) -> bool
    {
        self.flags[seq].is_leaf()
    }

    #[inline]
    pub fn set_is_leaf(&mut self, seq: usize, value: bool)
    {
        self.flags[seq].set_is_leaf(value);
    }

    #[inline]
    pub fn is_inner_point(&self, seq: usize) -> bool
    {        
        self.level_min(seq) > 0
    }   

    #[inline]
    pub fn level_sum(&self, seq: usize) -> u32
    {
        self.level[seq*self.num_inputs..(seq+1)*self.num_inputs].iter().map(|&i| i as u32).sum()
    }

    #[inline]
    pub fn level_max(&self, seq: usize) -> u8
    {
        *self.level[seq*self.num_inputs..(seq+1)*self.num_inputs].iter().max().unwrap()
    }

    #[inline]
    pub fn level_min(&self, seq: usize) -> u8
    {
        *self.level[seq*self.num_inputs..(seq+1)*self.num_inputs].iter().min().unwrap()
    }
    
    pub fn left_child(&self, seq: usize, dim: usize) -> GridPoint
    {
        let mut r = GridPoint::from((&self.index[seq*self.num_inputs..(seq+1)*self.num_inputs],
             &self.level[seq*self.num_inputs..(seq+1)*self.num_inputs], &self.flags[seq]));
        if r.index[dim] == 0
        {
            r.index[dim] = u32::MAX;
            return r;
        }
        r.index[dim] = 2*self.index[seq*self.num_inputs + dim] - 1;
        r.level[dim] += 1;
        r
    }
    pub fn right_child(&self, seq: usize, dim: usize) -> GridPoint
    {
        let mut r = GridPoint::from((&self.index[seq*self.num_inputs..(seq+1)*self.num_inputs], &self.level[seq*self.num_inputs..(seq+1)*self.num_inputs], &self.flags[seq]));
        r.index[dim] = 2*self.index[seq*self.num_inputs + dim] + 1;
        r.level[dim] += 1;
        r
    }
    ///
    /// returns an index with the top level in direction dim
    /// 
    pub fn root(&self, seq: usize, dim: usize) -> GridPoint
    {
        let mut r = GridPoint::from((&self.index[seq*self.num_inputs..(seq+1)*self.num_inputs], &self.level[seq*self.num_inputs..(seq+1)*self.num_inputs], &self.flags[seq]));
        r.index[dim] = 1;
        r.level[dim] = 1;
        r
    }

    pub fn insert_point(&mut self, mut point: GridPoint)
    {
         // make sure our is_inner flag is up-to-date...
        point.flags.update_is_inner(&point.level);
        let key: u64 = (&point).into();
        self.flags.push(point.flags);
        self.index.extend(point.index);
        self.level.extend(point.level);        
        self.map.insert(key, self.flags.len() - 1);
    }
   
    #[inline]
    pub fn update(&mut self, mut point: GridPoint, index: usize)
    {
        // make sure our is_inner flag is up-to-date...
        point.flags.update_is_inner(&point.level);
        let key: u64 = (&point).into();
        self.map.insert(key, index);
        self.index.chunks_exact_mut(self.num_inputs).nth(index).unwrap().copy_from_slice(&point.index);
        self.level.chunks_exact_mut(self.num_inputs).nth(index).unwrap().copy_from_slice(&point.level);   
        self.flags[index] = point.flags;
    }
    ///
    /// Return the nodes in the grid...
    /// 
    pub fn nodes(&self) -> NodeIterator<'_> {
        NodeIterator::new(&self)
    }

    ///
    /// Return the real coordinates for each node...
    /// 
    pub fn points(&self) -> PointIterator<'_>
    {
        PointIterator::new(&self)
    }

    pub fn generate_map(&mut self)
    {        
        let mut map = FxHashMap::default();
        for (i, node) in self.nodes().enumerate()
        {                
            map.insert(node.into(), i);
        }
        self.map = map;
    }
    #[inline]
    pub fn map(&self) -> &FxHashMap<u64, usize>
    {
        &self.map
    }
    #[inline]
    pub fn reset_map(&mut self)
    {
        self.map.clear();
    }
    #[inline]
    pub fn map_initialized(&self) -> bool
    {
        self.len() == self.map.len()
    }
    #[inline]
    pub fn contains(&self, point: &GridPoint) -> bool
    {
        self.map.contains_key(&point.into())
    }
    #[inline]
    pub fn get_index(&self, index: &GridPoint) -> Option<usize>
    {
        self.map.get(&index.into()).copied()
    }

    pub fn remove(&mut self, points_to_keep: &indexmap::IndexSet<usize>)
    {
        let mut indices = Vec::with_capacity(points_to_keep.len()*self.num_inputs);
        let mut levels = Vec::with_capacity(points_to_keep.len()*self.num_inputs);
        let mut flags =  Vec::with_capacity(points_to_keep.len());
        for &i in points_to_keep
        {
            indices.extend(&self.index[i*self.num_inputs..(i+1)*self.num_inputs]);
            levels.extend(&self.level[i*self.num_inputs..(i+1)*self.num_inputs]);
            flags.push(self.flags[i]);            
        }
        self.index = indices;
        self.level = levels;
        self.flags = flags;
        self.generate_map();
        self.update_leaves();
    }

    fn update_leaves(&mut self)
    {
        #[allow(clippy::needless_range_loop)]
        for i in 0..self.len()
        {
            //let point = &mut list[i];
            let mut is_leaf = true;
            let point = self.point(i);
            for dim in 0..self.num_inputs
            {                
                if self.level[i*self.num_inputs + dim] > 0
                {                    
                    // Check if this point has any children. If not it is a leaf.
                    if self.contains(&point.left_child(dim)) ||
                        self.contains(&point.right_child(dim)) 
                    {
                        is_leaf = false;
                        break; 
                    }                
                }
                else
                {                    
                    // don't remove boundary nodes that are used by other nodes...
                    if self.contains(&point.root(dim)) 
                    {
                        is_leaf = false;
                        break;
                    }
                }
            }                
            self.flags[i].set_is_leaf(is_leaf);
        }
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
        let mut iterator = DynamicHashMapGridIterator::new(&self);
        let offset =  dim * self.len();
        let active_index = offset + seq;
        if !array[active_index].inner.is_complete()
        {
            let node_index= self.point(seq); 
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
    #[inline]
    pub fn index_of(&self, point: &GridPoint) -> Option<usize>
    {
        self.map.get(&point.into()).copied()
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

    pub fn unit_coordinate(&self, seq: usize) -> Vec<f64>
    {
        let mut coor = vec![0.0; self.num_inputs];
        #[allow(clippy::needless_range_loop)]
        for d in 0..self.num_inputs
        {
            coor[d] = self.index[seq*self.num_inputs + d] as f64 / (1 << self.level[seq*self.num_inputs + d]) as f64;
        }
        coor
    }

    // #[inline]
    // pub fn points(&self) -> Vec<f64>
    // {
    //     let mut list = Vec::with_capacity(self.len()*self.num_inputs);
    //     let mut coor = vec![0.0; self.num_inputs];
    //     for (index, level, _flag) in self.points()
    //     {
    //         for d in 0..self.num_inputs
    //         {
    //             coor[d] = index[d] as f64 / (1 << level[d]) as f64;
    //         }
    //         self.bounding_box.to_real_coordinate_in_place(&mut coor);
    //         list.extend(coor);
    //     }
    //     list
    // }
}

pub struct NodeIterator<'a> {
    storage: &'a SparseGridData,
    current_seq: usize,
}
impl<'a> NodeIterator<'a>
{
    pub fn new( storage: &'a SparseGridData) -> Self
    {
        Self { storage, current_seq: 0 }
    }
}

impl<'a> Iterator for NodeIterator<'a> {
    type Item = GridPointRef<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_seq < self.storage.len() {
            let start = self.current_seq * self.storage.num_inputs;
            let end = start + self.storage.num_inputs;
            self.current_seq += 1;
            Some((&self.storage.index[start..end], &self.storage.level[start..end], &self.storage.flags[self.current_seq - 1]).into())
        } else {
            None
        }
    }
}

pub struct PointIterator<'a> {
    storage: &'a SparseGridData,
    current_seq: usize,
}
impl<'a> PointIterator<'a>
{
    pub fn new( storage: &'a SparseGridData) -> Self
    {
        Self { storage, current_seq: 0 }
    }
}

impl<'a> Iterator for PointIterator<'a> {
    type Item = Vec<f64>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_seq < self.storage.len() {
            let mut point = self.storage.unit_coordinate(self.current_seq);
            self.storage.bounding_box.to_real_coordinate_in_place(&mut point);
            self.current_seq += 1;
            Some(point)
        } else {
            None
        }
    }
}