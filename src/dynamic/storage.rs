use std::{hash::{Hash, Hasher}, u32};
use nohash_hasher::BuildNoHashHasher;
use rustc_hash::{FxHashMap, FxHasher};
use serde::{Deserialize, Serialize};
use serde_with::serde_as;
use bitfield_struct::bitfield;
use crate::{dynamic::iterators::dynamic_grid_iterator::{DynamicHashMapGridIterator, GridIteratorT}, errors::SGError};

use crate::adjacency_data::{NodeAdjacency, NodeAdjacencyData};
pub type FastU64Map<V> = std::collections::HashMap<u64, V, BuildNoHashHasher<u64>>;

#[bitfield(u8, new=false)]
#[derive(Serialize, Deserialize, PartialEq, Eq)]
pub struct GridPointFlags
{
    pub is_leaf: bool,
    pub is_inner: bool,
    #[bits(6)]
    pub _empty: u8
}

#[cfg(feature = "rkyv")]
impl rkyv::Archive for GridPointFlags {
    type Archived = u8;
    type Resolver = ();

    fn resolve(&self, _resolver: Self::Resolver, out: rkyv::Place<Self::Archived>) {
        out.write(self.0);
    }
}
#[cfg(feature = "rkyv")]
impl<S: rkyv::rancor::Fallible + ?Sized> rkyv::Serialize<S> for GridPointFlags {
    fn serialize(&self, _serializer: &mut S) -> Result<Self::Resolver, S::Error> {
        Ok(())
    }
}
#[cfg(feature = "rkyv")]
impl<D: rkyv::rancor::Fallible + ?Sized> rkyv::Deserialize<GridPointFlags, D> for u8 {
    fn deserialize(&self, _deserializer: &mut D) -> Result<GridPointFlags, D::Error> {
        Ok(GridPointFlags(*self))
    }
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
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
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
        *self.level.iter().max().unwrap_or(&0)
    }
    pub fn level_min(&self) -> u8
    {
        *self.level.iter().min().unwrap_or(&0)
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
        *self.level.iter().max().unwrap_or(&0)
    }
    pub fn level_min(&self) -> u8
    {
        *self.level.iter().min().unwrap_or(&0)
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
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
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
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
pub struct SparseGridData
{
    pub bounding_box: BoundingBox,
    pub(crate) index: Vec<u32>,
    pub(crate) level: Vec<u8>,
    pub(crate) flags: Vec<GridPointFlags>,
    pub(crate) num_inputs: usize,
    pub(crate) map: FastU64Map<u32>,
    pub(crate) adjacency_data: NodeAdjacencyData,
    pub(crate) has_boundary: bool,
    pub(crate) num_outputs: usize,
}

impl SparseGridData
{
    pub fn new(num_inputs: usize, num_outputs: usize) -> Self
    {
        Self { bounding_box: BoundingBox::with_dim(num_inputs), index: Vec::new(), level: Vec::new(), flags: Vec::new(), num_inputs, num_outputs, map: FastU64Map::default(), adjacency_data: NodeAdjacencyData::default(), has_boundary: false }
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
        if let Some(seq) = self.map.get(&point.into())
        {
            let seq = *seq as usize;
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
        *self.level[seq*self.num_inputs..(seq+1)*self.num_inputs].iter().max().unwrap_or(&0)
    }

    #[inline]
    pub fn level_min(&self, seq: usize) -> u8
    {
        *self.level[seq*self.num_inputs..(seq+1)*self.num_inputs].iter().min().unwrap_or(&0)
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
        self.map.insert(key, self.flags.len() as u32 - 1);
    }
   
    #[inline]
    pub fn update(&mut self, mut point: GridPoint, index: usize) -> Result<(), SGError>
    {
        // make sure our is_inner flag is up-to-date...
        point.flags.update_is_inner(&point.level);
        let key: u64 = (&point).into();
        self.map.insert(key, index as u32);
        self.index.chunks_exact_mut(self.num_inputs).nth(index).ok_or_else(||SGError::InvalidIndex)?.copy_from_slice(&point.index);
        self.level.chunks_exact_mut(self.num_inputs).nth(index).ok_or_else(||SGError::InvalidIndex)?.copy_from_slice(&point.level);   
        self.flags[index] = point.flags;
        Ok(())
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
        let mut map = FastU64Map::default();
        for (i, node) in self.nodes().enumerate()
        {                
            map.insert(node.into(), i as u32);
        }
        self.map = map;
    }
    #[inline]
    pub fn map(&self) -> &FastU64Map<u32>
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
        self.map.get(&index.into()).map(|&v| v as usize)
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
        self.generate_adjacency_data();
    }

    fn update_leaves(&mut self)
    {
        #[allow(clippy::needless_range_loop)]
        for i in 0..self.len()
        {
            let point = self.point(i);
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
                    // Boundary node in this dimension - check for level-1 child
                    let mut child = point.clone();
                    child.level[dim] = 1;
                    child.index[dim] = 1;
                    if self.map.contains_key(&child.into())
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
        let total_size = self.num_inputs * self.len();
        let mut array = vec![NodeAdjacency::default(); total_size];
        let mut left_zero = vec![u32::MAX; total_size];
        let mut right_zero = vec![u32::MAX; total_size];
        if self.map.len() != self.len()
        {
            self.generate_map();
        }
        for dim in 0..self.num_inputs
        {
            for i in 0..self.len()
            {
                self.generate_adjacency_data_for_index(&mut array, &mut left_zero, &mut right_zero, i, dim);     
            }
        }
        self.adjacency_data.zero_index = self.index_of(&GridPoint::zero_index(self.num_inputs)).unwrap_or(usize::MAX);
        self.adjacency_data.data = array;
        self.adjacency_data.left_zero = left_zero;
        self.adjacency_data.right_zero = right_zero;
    }

    fn generate_adjacency_data_for_index(&mut self, array:&mut [NodeAdjacency], left_zero: &mut [u32], right_zero: &mut [u32], seq: usize, dim: usize)
    {
        let mut iterator = DynamicHashMapGridIterator::new(&self);
        let offset =  dim * self.len();
        let active_index = offset + seq;
        let node_index= self.point(seq); 
        
        // Only process adjacency links if not already complete
        if !array[active_index].is_complete()
        {            
            iterator.set_index(node_index.clone());
            iterator.up(dim);       
            if let Some(parent_index) = iterator.index()  // parent exists
            {
                // assign parent
                array[active_index].set_up(parent_index as i64 - seq as i64);      
                array[active_index].set_has_parent(true);
            }
            iterator.set_index(node_index.clone());
            // get left child
            iterator.left_child(dim);
            if let Some(lc_index) = iterator.index()
            {
                // assign left child - down_left offset
                array[active_index].set_down_left(lc_index as i64 - seq as i64);
            }
            
            iterator.set_index(node_index.clone());
            // get right child        
            iterator.right_child(dim);
            if let Some(rc_index) = iterator.index()
            {
                // assign right child - down_right offset
                array[active_index].set_down_right(rc_index as i64 - seq as i64);
            }             
        }
        
        // Always populate boundary and level indices (even for complete nodes)
        iterator.set_index(node_index.clone());
        iterator.reset_to_left_level_zero(dim);
        if let Some(lzero) = iterator.index() {
            left_zero[active_index] = lzero as u32;
        }
        
        iterator.set_index(node_index.clone());
        iterator.reset_to_right_level_zero(dim);
        if let Some(rzero) = iterator.index() {
            right_zero[active_index] = rzero as u32;
        }
        iterator.set_index(node_index);
        iterator.reset_to_level_one(dim);
        if let Some(level_one_idx) = iterator.index() {
            array[active_index].set_level_one(level_one_idx as u32);
        }
        else
        {
            array[active_index].set_level_one(u32::MAX);
        }
    }
    #[inline]
    pub fn index_of(&self, point: &GridPoint) -> Option<usize>
    {
        self.map.get(&point.into()).map(|v|*v as usize)
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
    pub storage: &'a SparseGridData,
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