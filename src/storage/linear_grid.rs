use std::{collections::HashSet, hash::{Hash, Hasher}, ops::{Index, IndexMut}};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use serde_with::serde_as;

#[serde_as]
#[derive(Copy, Clone, Eq, Serialize, Deserialize, Debug)]
pub struct GridPoint<const D: usize>
{
    #[serde_as(as = "[_; D]")]
    pub level: [u32; D],
    #[serde_as(as = "[_; D]")]
    pub index: [u32; D],
    pub(crate) is_leaf: bool,
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
        Self { level: [1; D], index: [1; D], is_leaf: false }
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

impl<const D: usize> GridPoint<D>
{
    pub fn new (level: [u32; D], index: [u32; D], is_leaf: bool) -> Self
    {
         Self { level, index, is_leaf}
    }   
    pub fn is_leaf(&self) -> bool
    {
        self.is_leaf
    }
    /// 
    /// This is an inner point if no indices are zero...
    /// 
    pub fn is_inner_point(&self) -> bool
    {        
        self.level_min() > 0
    }   
    pub fn level_sum(&self) -> u32
    {
        self.level.iter().sum()
    }
    pub fn level_max(&self) -> u32
    {
        *self.level.iter().max().unwrap()
    }
    pub fn level_min(&self) -> u32
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

}


impl<const D: usize> PartialEq for GridPoint<D>
{
    fn eq(&self, other: &Self) -> bool {
        self.level == other.level && self.index == other.index
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

///
/// Storage for a linear grid. This differs from `SparseGridData` 
/// in that it contains a map to allow retrieval the storage index
/// for a given grid point.
/// 
#[derive(Default, Clone, Serialize, Deserialize)]
pub struct SparseGridStorage<const D: usize>
{
    pub map: FxHashMap<GridPoint<D>, usize>,
    pub data: SparseGridData<D>
}

impl<const D: usize> SparseGridStorage<D>
{
    pub fn new(data: SparseGridData<D>) -> Self
    {
        let map: FxHashMap<GridPoint<D>, usize> = FxHashMap::from_iter(data.iter().enumerate().map(|(i, item)|(*item, i)));
        Self { data, map }
    }
    #[inline(always)]
    pub fn len(&self) -> usize
    {
        self.data.list.len()
    }

    #[inline(always)]
    pub fn has_boundary(&self) -> bool
    {
        self.data.has_boundary
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    #[inline]
    pub fn contains(&self, point: &GridPoint<D>) -> bool
    {
        self.map.contains_key(point)
    }
    #[inline]
    pub fn get_mut(&mut self, point: &GridPoint<D>) -> Option<&mut GridPoint<D>>
    {
        match self.map.get(point)
        {
            Some(&seq) => Some(&mut self.data.list[seq]),
            None => None
        }
    }

    #[inline(always)]
    pub fn list(&self) -> &Vec<GridPoint<D>>
    {
        &self.data.list
    }
    
    #[inline(always)]
    pub fn list_mut(&mut self) -> &mut Vec<GridPoint<D>>
    {
        &mut self.data.list
    }
    pub fn store(&mut self, index: GridPoint<D>)
    {
        self.map.insert(index, self.data.list.len());
        self.data.list.push(index);     
    }
    #[inline]
    pub fn bounding_box(&self) -> Option<BoundingBox<D>>
    {
        self.data.bounding_box
    }
    #[inline]
    pub fn bounding_box_mut(&mut self) -> &mut Option<BoundingBox<D>>
    {
        &mut self.data.bounding_box
    }
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, GridPoint<D>>
    {
        self.data.list.iter()
    }
    #[inline]
    pub fn sequence_number(&self, index: &GridPoint<D>) -> Option<usize>
    {
        self.map.get(index).copied()
    }

    pub fn remove(&mut self, points: &HashSet<usize>)
    {
        let mut new_list = Vec::with_capacity(self.len() - points.len());
        for i in 0..self.len()
        {
            if points.contains(&i)
            {
                continue;
            }
            new_list.push(self.data.list[i]);
        }
        self.data.list = new_list;
        self.map = FxHashMap::from_iter(self.data.list.iter().enumerate().map(|(i, item)|(*item, i)));
        Self::update_leaves(&mut self.data.list, &self.map);
    }

    fn update_leaves(list: &mut [GridPoint<D>], map: &FxHashMap<GridPoint<D>, usize>)
    {
        #[allow(clippy::needless_range_loop)]
        for i in 0..list.len()
        {
            let point = &mut list[i];
            let mut is_leaf = true;

            for dim in 0..D
            {                
                if point.level[dim] > 0
                {                    
                    // Check if this point has any children. If not it is a leaf.
                    if map.contains_key(&point.left_child(dim)) ||
                        map.contains_key(&point.right_child(dim)) 
                    {
                        is_leaf = false;
                        break; 
                    }                
                }
                else
                {                    
                    // don't remove boundary nodes that are used by other nodes...
                    if map.contains_key(&point.root(dim)) 
                    {
                        is_leaf = false;
                    }
                }
            }                
            point.is_leaf = is_leaf;            
        }
    }

    ///
    /// Inserts a point
    /// 
    pub fn insert_point(&mut self, point: GridPoint<D>) -> usize
    {
        self.data.list.push(point);
        let value = self.data.list.len() - 1;
        self.map.insert(point, value);
        value
    }

    pub fn update(&mut self, point: GridPoint<D>, index: usize)
    {
        self.map.remove(&point);
        self.map.insert(point, index);
        self.data.list[index] = point;
    }

    pub fn points(&self) -> Vec<[f64; D]>
    {
        let mut list = Vec::new();
        for index in &self.data.list
        {
            let mut point = index.unit_coordinate();
            if let Some(bbox) = self.data.bounding_box
            {
                point = bbox.to_real_coordinate(&point);
            }
            list.push(point)
        }
        list
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct SparseGridData<const D: usize>
{
    pub(crate) list: Vec<GridPoint<D>>,
    pub(crate) bounding_box: Option<BoundingBox<D>>,
    pub(crate) has_boundary: bool,
}


impl<const D: usize> Default for SparseGridData<D>
{
    fn default() -> Self {
        Self { list: Default::default(), bounding_box: Default::default(), has_boundary: false }
    }
}

impl<const D: usize,  Idx: std::slice::SliceIndex<[GridPoint<D>]>> Index<Idx> for SparseGridStorage<D> 
{
    type Output=Idx::Output;    
    fn index(&self, index: Idx) -> &Self::Output 
    {
        &self.data.list[index]
    }
}
impl<const D: usize, Idx: std::slice::SliceIndex<[GridPoint<D>]>> IndexMut<Idx> for SparseGridStorage<D>
{
    fn index_mut(&mut self, index: Idx) -> &mut Self::Output {
        &mut self.data.list[index]
    }
}
impl<const D: usize> SparseGridData<D>
{
    pub fn new() -> Self
    {
        Self { ..Default::default() }
    }
    
    pub fn is_empty(&self) -> bool
    {
        self.list.is_empty()
    }
    pub fn len(&self) -> usize
    {
        self.list.len()
    }   

    pub fn iter(&self) -> std::slice::Iter<'_, GridPoint<D>>
    {
        self.list.iter()
    }
}

pub struct PointIterator<'a, const D: usize>
{
    iterator: &'a mut std::slice::Iter<'a, GridPoint<D>>,
    bounding_box: Option<BoundingBox<D>>,
}

impl<'a, const D: usize>  PointIterator<'a, D>
{
    pub fn new( iterator: &'a mut std::slice::Iter<'a, GridPoint<D>>, bounding_box: Option<BoundingBox<D>>) -> Self
    {
        Self { iterator, bounding_box }
    }
}

impl<const D: usize> Iterator for PointIterator<'_, D>
{
    type Item=[f64; D];

    fn next(&mut self) -> Option<Self::Item> {

        if let Some(index) = self.iterator.next()
        {
            let mut point = index.unit_coordinate();
            if let Some(bbox) = self.bounding_box
            {
                point = bbox.to_real_coordinate(&point);
            }
            Some(point)
        }
        else
        {
            None
        }
    }
}