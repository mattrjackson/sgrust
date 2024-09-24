use crate::{basis::{base::Basis, global::GlobalBasis}, grids::combination_grid::TensorSelectionStrategy};

use super::step_iterator::GridType;

pub struct AdvancedStepIterator
{
    ndim: usize,
    max_levels: Vec<u32>,
    level_bound: u32,
    exactness_bound: u32,
    exactness_sum: u32,
    index_sum: u32,
    index_head: Vec<u32>,
    first: bool,    
    basis: Vec<GlobalBasis>,
    selection_strategy: TensorSelectionStrategy,
}
impl AdvancedStepIterator
{
    pub fn new(max_levels: &[u32], grid_type: GridType, basis: Vec<GlobalBasis>, mut selection_strategy: TensorSelectionStrategy, exactness_bound: u32) -> Self
    {
        let ndim = max_levels.len();
        let level_bound = if grid_type == GridType::Sparse {  *max_levels.iter().max().unwrap() } else { max_levels.iter().sum() };        
        // For a full grid, there is no point in using an exactness test
        selection_strategy = if grid_type == GridType::Full { TensorSelectionStrategy::Level } else { selection_strategy };
        
        Self { level_bound, max_levels: max_levels.to_owned(), index_sum: 0, ndim,  index_head: vec![0; ndim], first: true,
            basis, selection_strategy, exactness_bound, exactness_sum: 0 }
    } 
    #[inline]
    fn in_bounds(&self) ->bool
    {
        match self.selection_strategy
        {
            TensorSelectionStrategy::Level => self.level_bound > self.index_sum,
            TensorSelectionStrategy::InterpolationExactness | TensorSelectionStrategy::QuadratureExactness => self.exactness_bound >= self.exactness_sum,
        }
    }
    
    fn sum_exactness(&self) -> u32
    {
        let mut sum = 0;
        
        for i in 0..self.ndim
        {
            if self.index_head[i] > 0
            {
                sum += match self.selection_strategy
                {
                    TensorSelectionStrategy::Level => {0},
                    TensorSelectionStrategy::QuadratureExactness => self.basis[i].quadrature_exactness(self.index_head[i] - 1) + 1,
                    TensorSelectionStrategy::InterpolationExactness => self.basis[i].interpolation_exactness(self.index_head[i] - 1) + 1,
                }
            }
        }
        sum
    }

}
impl Iterator for AdvancedStepIterator
{
    type Item = Vec<u32>;

    fn next(&mut self) -> Option<Self::Item> {
        loop
        {
            if self.first
            {
                self.first = false;
                return Some(self.index_head.clone());
            }
            self.index_head[self.ndim - 1] += 1;
            self.exactness_sum = self.sum_exactness();
            self.index_head[self.ndim - 1] -= 1;
            if self.in_bounds() && self.max_levels[self.ndim-1] > self.index_head[self.ndim-1]
            {
                self.index_sum += 1;
                self.index_head[self.ndim - 1] += 1;
                self.exactness_sum = self.sum_exactness();                
            }        
            else
            {
                let mut dim = self.ndim - 1;
                while self.index_head[dim] == 0 && dim != usize::MAX 
                {
                    dim -= 1;
                }
                match dim
                {
                    0 => 
                    {
                        self.index_head[0] = 0;
                        self.index_sum = 0;         
                        self.exactness_sum = self.sum_exactness(); 
                        return None;
                    }
                    _ =>
                    {
                        self.index_sum -= self.index_head[dim] - 1;                        
                        self.index_head[dim] = 0;
                        self.index_head[dim-1] += 1;
                        self.exactness_sum = self.sum_exactness();
                    }
                }
            }
            let mut break_loop = true;
            for i in 0..self.ndim
            {
                if self.index_head[i] > self.max_levels[i]
                {
                    break_loop = false;
                    break;
                }
            }
            if break_loop { break; }
        }
        Some(self.index_head.clone())
    
    }
}