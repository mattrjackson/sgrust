use crate::dynamic::storage::{GridPoint, SparseGridData};

pub trait GridIteratorT
{
    fn point(&self) -> &GridPoint;    
    fn point_index(&self, dim: usize) -> u32;
    fn index(&self) -> Option<usize>;
    fn reset_to_level_zero(&mut self) -> bool;
    fn reset_to_left_level_zero(&mut self, dim: usize) -> bool;
    fn reset_to_right_level_zero(&mut self, dim: usize) -> bool;
    fn reset_to_level_one(&mut self, dim: usize) -> bool;
    fn left_child(&mut self, dim: usize) -> bool;
    fn right_child(&mut self, dim: usize) -> bool;
    fn up(&mut self, dim: usize) -> bool;
    fn is_inner_point(&self) -> bool;
    fn is_leaf(&self) -> bool;

}

pub struct DynamicHashMapGridIterator<'a>
{
    pub(crate) storage: &'a  SparseGridData,
    index: GridPoint,
    seq: Option<usize>,
}

impl<'a> DynamicHashMapGridIterator<'a>
{

    pub(crate) fn new(storage: &'a SparseGridData ) -> Self
    {        
        let point = storage.point(0);
        let seq = storage.index_of(&point);
        Self { storage, index: point, seq}
    }
    
    pub(crate) fn set_index(&mut self, point: GridPoint)
    {
        self.index = point;
        
        self.seq = self.storage.index_of(&self.index);
    }

    
    pub(crate) fn step_left(&mut self, dim: usize)
    {
        let i = self.index.index[dim];     
        if i < 2 
        {
            self.seq = None;
            return;
        }   
        self.index.index[dim] = i - 2;
        
        self.seq = self.storage.index_of(&self.index);
    }
    pub(crate) fn step_right(&mut self, dim: usize)
    {
        let i = self.index.index[dim];
        self.index.index[dim] = i + 2;
        self.seq = self.storage.index_of(&self.index);
    }

    #[allow(unused)]
    pub(crate) fn is_left_leaf(&self, dim: usize) -> bool
    {
        let i = self.index.index[dim];
        let l = self.index.level[dim];
        let mut index = self.index.clone();
        index.index[dim] = 2 * i - 1;
        index.level[dim] = l + 1;        
        self.storage.contains(&index)
    }

    #[allow(unused)]
    pub(crate) fn is_right_leaf(&self, dim: usize) -> bool
    {
        let i = self.index.index[dim];
        let l = self.index.level[dim];
        let mut index = self.index.clone();
        index.index[dim] = 2 * i + 1;
        index.level[dim] = l + 1;        
        self.storage.contains(&index)
    }

    #[inline(always)]
    pub(crate) fn seq(&self) -> Option<usize>
    {
        self.seq
    }

    #[allow(unused)]
    pub(crate) fn get_grid_depth(&mut self, dim: usize) -> usize
    {
        let mut depth = 1;
        let orig_index = self.index.index[dim];
        let orig_level = self.index.level[dim];
        loop
        {
            
            if self.is_left_leaf(dim) || self.is_right_leaf(dim)
            {
                depth += 1;
            }
            else
            {                
                let cur_index = self.index.index[dim];                
                let mut index_found = false;
                for i in (cur_index+2..(1<< depth)).step_by(2)
                {
                    self.index.index[dim] = i;
                    if self.seq.is_some()
                    {
                        if self.is_left_leaf(dim)
                        {
                            depth += 1;
                            self.left_child(dim);
                            index_found = true;
                            break;
                        }
                        else if self.is_right_leaf(dim)
                        {
                            depth += 1;
                            self.right_child(dim);
                            index_found = true;
                            break;
                        }
                    }
                }
                if !index_found
                {
                    break;
                }
            }
        }
        self.index.index[dim] = orig_index;
        self.index.level[dim] = orig_level;
        depth
    }
}

impl GridIteratorT for DynamicHashMapGridIterator<'_>
{
    #[inline(always)]
    fn point(&self) -> &GridPoint
    {
        &self.index
    }
    
    fn index(&self) ->  Option<usize>
    {
        self.seq
    }

    fn reset_to_level_zero(&mut self) -> bool
    {        
        self.index.index.fill(0);
        self.index.level.fill(0);
        self.seq = self.storage.index_of(&self.index);
        self.seq.is_some()
    }
    fn reset_to_left_level_zero(&mut self, dim: usize) -> bool
    {
        self.index.index[dim] = 0;
        self.index.level[dim] = 0;
        
        self.seq = self.storage.index_of(&self.index);        
        self.seq.is_some()
    }
    fn reset_to_right_level_zero(&mut self, dim: usize) -> bool
    {
        self.index.level[dim] = 0;
        self.index.index[dim] = 1;       
        
        self.seq = self.storage.index_of(&self.index);  
        self.seq.is_some()      
    } 
    fn reset_to_level_one(&mut self, dim: usize) -> bool
    {
        self.index.level[dim] = 1;
        self.index.index[dim] = 1;
        
        self.seq = self.storage.index_of(&self.index);
        self.seq.is_some()
    }
    fn left_child(&mut self, dim: usize) -> bool
    {
        let i = self.index.index[dim];
        if i == 0
        {
            self.seq = None;
            return false;
        }
        let l = self.index.level[dim];
        self.index.level[dim] = l + 1;
        self.index.index[dim] = (2 * i) - 1;
        
        self.seq = self.storage.index_of(&self.index);
        self.seq.is_some()
    }
    fn right_child(&mut self, dim: usize) -> bool
    {
        let i = self.index.index[dim];
        let l = self.index.level[dim];
        self.index.level[dim] = l + 1;
        self.index.index[dim] = 2 * i + 1;
        
        self.seq = self.storage.index_of(&self.index);
        self.seq.is_some()
    }
    fn up(&mut self, dim: usize) -> bool
    {
        let mut i = self.index.index[dim];
        let l = self.index.level[dim];
        if l == 0
        {
            self.seq = None;
            return false;
        }
        i /= 2;
        i += if i % 2 == 0 {1} else {0};       
        self.index.level[dim] = l - 1;
        self.index.index[dim] = i;
        
        self.seq = self.storage.index_of(&self.index);
        self.seq.is_some()
    }

    
    #[allow(unused)]
    fn is_inner_point(&self) -> bool
    {
        self.index.is_inner_point()
    }

    fn is_leaf(&self) -> bool
    {
        if let Some(seq) = self.seq
        {
            self.storage.flags[seq].is_leaf()
        }
        else
        {
            true
        }
    }

    fn point_index(&self, dim: usize) -> u32
    {
        self.index.index[dim]
    }
}