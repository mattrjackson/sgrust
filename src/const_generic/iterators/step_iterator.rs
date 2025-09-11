#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum GridType
{
    Sparse,
    Full,
}
pub struct StepIterator
{
    ndim: usize,
    max_levels: Vec<u32>,
    level_bound: u32,
    index_sum: u32,
    index_head: Vec<u32>,
    first: bool,    
}
impl StepIterator
{
    pub fn new(max_levels: &[u32], grid_type: GridType) -> Self
    {
        let ndim = max_levels.len();
        let level_bound = if grid_type == GridType::Sparse {  *max_levels.iter().max().unwrap_or(&0) } else { max_levels.iter().sum() };

        Self { level_bound, max_levels: max_levels.to_owned(), index_sum: 0, ndim,  index_head: vec![0; ndim], first: true }
    } 
    pub fn last_dimension_count(&self) -> u32
    {
        (self.level_bound - self.index_sum + 1).min(self.max_levels[self.ndim-1]+1)
    }
}
impl Iterator for StepIterator
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
            if self.level_bound > self.index_sum && self.max_levels[self.ndim-1] > self.index_head[self.ndim-1]
            {
                self.index_sum += 1;
                self.index_head[self.ndim - 1] += 1;
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
                        return None;
                    }
                    _ =>
                    {
                        self.index_sum -= self.index_head[dim] - 1;
                        self.index_head[dim] = 0;
                        self.index_head[dim-1] += 1;
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