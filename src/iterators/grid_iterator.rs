use crate::storage::linear_grid::{GridPoint, SparseGridStorage};



pub(crate) struct GridIterator<'a, const D: usize>
{
    pub(crate) storage: &'a SparseGridStorage<D>,
    index: GridPoint<D>,
    seq: Option<usize>,
}

impl<'a, const D: usize> GridIterator<'a, D>
{

    pub(crate) fn new(storage: &'a SparseGridStorage<D> ) -> Self
    {
        let index = GridPoint::default();
        Self { storage, index, seq: storage.sequence_number(&index) }
    }
    #[inline(always)]
    pub(crate) fn index(&self) -> &GridPoint<D>
    {
        &self.index
    }
    #[inline(always)]
    pub fn valid_seq(&self) -> bool
    {
        self.seq.is_some()
    }
    pub(crate) fn set_index(&mut self, point: GridPoint<D>)
    {
        self.index = point;
        
        self.seq = self.storage.sequence_number(&self.index);
    }

    pub(crate) fn reset_to_level_zero(&mut self)
    {        
        self.index.index = [0; D];      
        self.index.level = [0; D];  
        
        self.seq = self.storage.sequence_number(&self.index);
    }
    pub(crate) fn reset_to_left_level_zero(&mut self, dim: usize)
    {
        self.index.index[dim] = 0;
        self.index.level[dim] = 0;
        
        self.seq = self.storage.sequence_number(&self.index);        
    }
    pub(crate) fn reset_to_right_level_zero(&mut self, dim: usize)
    {
        self.index.level[dim] = 0;
        self.index.index[dim] = 1;       
        
        self.seq = self.storage.sequence_number(&self.index);        
    } 
    pub(crate) fn reset_to_level_one(&mut self, dim: usize)
    {
        self.index.level[dim] = 1;
        self.index.index[dim] = 1;
        
        self.seq = self.storage.sequence_number(&self.index);
    }
    pub(crate) fn left_child(&mut self, dim: usize)
    {
        let i = self.index.index[dim];
        if i == 0
        {
            self.seq = None;
            return;
        }
        let l = self.index.level[dim];
        self.index.level[dim] = l + 1;
        self.index.index[dim] = (2 * i) - 1;
        
        self.seq = self.storage.sequence_number(&self.index);
    }
    pub(crate) fn right_child(&mut self, dim: usize)
    {
        let i = self.index.index[dim];
        let l = self.index.level[dim];
        self.index.level[dim] = l + 1;
        self.index.index[dim] = 2 * i + 1;
        
        self.seq = self.storage.sequence_number(&self.index);
    }
    pub(crate) fn up(&mut self, dim: usize)
    {
        let mut i = self.index.index[dim];
        let l = self.index.level[dim];
        if l == 0
        {
            self.seq = None;
            return;
        }
        i /= 2;
        i += if i % 2 == 0 {1} else {0};       
        self.index.level[dim] = l - 1;
        self.index.index[dim] = i;
        
        self.seq = self.storage.sequence_number(&self.index);
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
        
        self.seq = self.storage.sequence_number(&self.index);
    }
    pub(crate) fn step_right(&mut self, dim: usize)
    {
        let i = self.index.index[dim];        
        self.index.index[dim] = i + 2;
        
        self.seq = self.storage.sequence_number(&self.index);
    }

    #[allow(unused)]
    pub(crate) fn is_inner_point(&self) -> bool
    {
        self.index.is_inner_point()
    }

    pub(crate) fn hint(&self) -> bool
    {
        if let Some(seq) = self.seq
        {
            self.storage[seq].is_leaf()
        }
        else
        {
            true
        }
    }

    #[allow(unused)]
    pub(crate) fn hint_left(&self, dim: usize) -> bool
    {
        let i = self.index.index[dim];
        let l = self.index.level[dim];
        let mut index = self.index;
        index.index[dim] = 2 * i - 1;
        index.level[dim] = l + 1;        
        self.storage.contains(&index)
    }

    #[allow(unused)]
    pub(crate) fn hint_right(&self, dim: usize) -> bool
    {
        let i = self.index.index[dim];
        let l = self.index.level[dim];
        let mut index = self.index;
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
            
            if self.hint_left(dim) || self.hint_right(dim)
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
                        if self.hint_left(dim)
                        {
                            depth += 1;
                            self.left_child(dim);
                            index_found = true;
                            break;
                        }
                        else if self.hint_right(dim)
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