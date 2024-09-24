#[derive(Clone)]
pub struct BoundedSumIterator<const D: usize>
{
    max_levels: [usize; D],
    level_bound: usize,
    index_sum: usize,
    index_head: [usize; D],
    first: bool

}
impl<const D: usize> Default for BoundedSumIterator<D>
{
    fn default() -> Self {
        Self { max_levels: [0; D], level_bound: Default::default(), index_sum: Default::default(), index_head: [0; D], first: Default::default() }
    }
}
impl<const D: usize> BoundedSumIterator<D>
{
    pub fn new( max_levels: [usize; D]) -> Self
    {
        Self { level_bound: *max_levels.iter().max().unwrap(), max_levels, index_sum: 0, index_head: [0; D], first: true }
    } 
    pub fn last_dimension_count(&self) -> usize
    {
        (self.level_bound - self.index_sum + 1).min(self.max_levels[D-1]+1)
    }

    pub fn num_levels(&self) -> usize
    {
        self.level_bound
    }
    pub fn max_levels(&self) -> [usize; D]
    {
        self.max_levels
    }
}
impl<const D: usize> Iterator for BoundedSumIterator<D>
{
    type Item = ([usize; D], usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.first
        {
            self.first = false;
            return Some((self.index_head, self.last_dimension_count()));
        }
        if self.level_bound > self.index_sum && self.max_levels[D-2] > self.index_head[D-2]
        {
            self.index_sum += 1;
            self.index_head[D - 2] += 1;
        }
        else
        {
            let mut dim = D - 2;
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
                },                
                _ =>
                {
                    self.index_sum -= self.index_head[dim] - 1;
                    self.index_head[dim] = 0;
                    self.index_head[dim-1] += 1;
                }
            }
        }
        Some((self.index_head, self.last_dimension_count()))
    }
}

//TODO: can probably remove this.
pub struct StepIterator<const D: usize>
{
    max_levels: [usize; D],
    level_bound: usize,
    index_sum: usize,
    index_head: [usize; D],
    first: bool,    
}

impl<const D: usize> StepIterator<D>
{
    #[allow(unused)]
    pub fn new( max_levels: [usize; D]) -> Self
    {
        Self { level_bound: *max_levels.iter().max().unwrap(), max_levels, index_sum: 0, index_head: [0; D], first: true }
    } 
    #[allow(unused)]
    pub fn last_dimension_count(&self) -> usize
    {
        (self.level_bound - self.index_sum + 1).min(self.max_levels[D-1]+1)
    }
}
impl<const D: usize> Iterator for StepIterator<D>
{
    type Item = [usize; D];

    fn next(&mut self) -> Option<Self::Item> {
        if self.first
        {
            self.first = false;
            return Some(self.index_head);
        }
        if self.level_bound > self.index_sum && self.max_levels[D-2] > self.index_head[D-2]
        {
            self.index_sum += 1;
            self.index_head[D - 1] += 1;
        }
        else
        {
            let mut dim = D - 1;
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
                },
                _ =>                
                {
                    self.index_sum -= self.index_head[dim] - 1;
                    self.index_head[dim] = 0;
                    self.index_head[dim-1] += 1;
                }
            }
        }
        Some(self.index_head)
    }
}

pub struct PointIterator<const D: usize>([f64; D], bool);

impl<const D: usize> Iterator for PointIterator<D>
{
    type Item = ([usize; D], usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.1
        {
            self.1 = false;
            Some(([0; D], self.0.len()))
        }
        else
        {
            None
        }        

    }
}

#[test]
fn test_iterator()
{
    let ndim = 3;
    let iterator = BoundedSumIterator::new([2,2,3]);
    for (item, last_dim_count) in iterator
    {        
        for i in 0..last_dim_count
        {
            #[allow(clippy::needless_range_loop)]
            for d in 0..ndim-1
            {
                print!("{:?},", item[d]);
            }
            println!("{i}");
        }
    }
}

#[test]
fn test_step_iterator()
{
    let ndim = 2;
    let iterator = StepIterator::new([2,2]);
    for item in iterator
    {  
        #[allow(clippy::needless_range_loop)]
        for d in 0..ndim
        {
            print!("{:?},", item[d]);
        }        
        println!();
    }
}