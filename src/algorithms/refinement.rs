use crate::storage::linear_grid::{GridPoint, SparseGridStorage};

///
/// This trait defines operations used for refinement or coarsening. These
/// two operations are never done simulataneously, but provide a common
/// interface to allow user-specified constraints to control either operation.
/// 
pub trait RefinementFunctor<const D: usize, const DIM_OUT: usize>
{
    ///
    /// Return criteria for determining refinement threshold
    /// 
    fn eval(&self, storage: &SparseGridStorage<D>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]], seq: usize) -> f64;

    ///
    /// Return threshold for refinement or coarsening
    /// 
    fn threshold(&self) -> f64;

    ///
    /// Returns the maximum number of points to be refined. If
    /// set to none there is no limit to the maximum number of points.
    /// 
    fn max_num_refined(&self) -> Option<usize>
    {
        None
    }

    ///
    /// Returns the maximum number of points that may be removed
    /// 
    fn max_num_removed(&self) -> Option<usize>
    {
        None
    }
}

pub trait SparseGridRefinement<const D: usize, const DIM_OUT: usize>
{
    ///
    /// Refine grid based on criteria computed using functor
    /// 
    fn refine(&self, storage: &mut SparseGridStorage<D>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]], functor: &dyn RefinementFunctor<D, DIM_OUT>) -> Vec<usize>;

    ///
    /// Returns the number of grid points that can be refined.
    ///      
    fn get_num_refinable_points(&self, storage: &SparseGridStorage<D>) -> usize;

    ///
    /// Refine a grid point along a single direction
    /// 
    fn refine_1d(&self, storage: &mut SparseGridStorage<D>, point: GridPoint<D>, dim: usize);

    fn create_point(&self, storage: &mut SparseGridStorage<D>, mut point: GridPoint<D>)
    {
        if !storage.contains(&point)
        {
            point.is_leaf = false;
            self.create_point(storage, point);
        }
        else if let Some(point) = storage.get_mut(&point)
        {
            point.is_leaf = false;
        }
    }
}

fn iterate_refinable_points<const D: usize, Op: FnMut((usize,  &GridPoint<D>))>(storage: &SparseGridStorage<D>, operation: &mut Op)
{
    for (seq, point) in storage.iter().enumerate()
    {
        let parent = point;
        let mut point = *point;
        for d in 0..D
        {               
            let index = point.index[d];
            let level = point.level[d];
            if level == 0
            {
                point.level[d] = 1;
                point.index[d] = 1;
                if !storage.contains(&point)
                {
                    operation((seq, parent));                    
                    break;
                }
            }
            else
            {
                // check if left child exists
                point.index[d] = 2 * index - 1;
                point.level[d] = level + 1;
                // Child doesn't exist. we can refine this node.
                if !storage.contains(&point)
                {
                    operation((seq, parent));
                    break;
                }

                // check if right child exists
                point.index[d] = 2 * index + 1;

                // Child doesn't exist. we can refine this node.
                if !storage.contains(&point)
                {
                    operation((seq, parent));
                    break;
                }                
                point.index[d] = index;
                point.level[d] = level;
            }
            point.level[d] = level;
            point.index[d] = index;
        }
    }
}



///
/// Base refinement struct. Boolean determines whether or not boundaries are enabled...
/// 
pub struct BaseRefinement<const D: usize, const DIM_OUT: usize>(pub bool);

impl<const D: usize, const DIM_OUT: usize> BaseRefinement<D, DIM_OUT>
{
    fn create_point(&self, storage: &mut SparseGridStorage<D>, point: GridPoint<D>)
    {
        for dim in 0..D
        {
            if !self.0  // no boundaries
            { 
                self.create_point_1d(point, storage, dim);
            }
            else // has boundaries
            {
                self.create_point_1d_with_boundary(point, storage, dim);
            }
        }
        storage.insert_point(point);

        // deal with boundaries
        if self.0
        {
            self.create_gridpoint_level_zero_consistency(storage, point);
        }
    }
    fn create_gridpoint_level_zero_consistency(&self, storage: &mut SparseGridStorage<D>, mut point: GridPoint<D>)
    {
        if D == 1 // only needed for D > 1
        {
            return;
        }
        for dim in 0..D
        {        
            let level = point.level[dim];               
            let index = point.index[dim];
            if level == 0
            {                
                point.level[dim] = 0;
                // loop through left, right boundaries
                for i in 0..2
                {
                    point.index[dim] = i;
                    if storage.contains(&point)
                    {
                        let leaf_l = point.is_leaf;
                        // check the boundary not being evaluated
                        point.index[dim] = if i == 0 { 1 } else { 0 };
                        if !storage.contains(&point)
                        {
                            let leaf_r = point.is_leaf;
                            point.is_leaf = leaf_l;
                            self.create_point(storage, point);
                            point.is_leaf = leaf_r;
                        }
                        else if let Some(point) = storage.get_mut(&point)
                        {
                            point.is_leaf = leaf_l;   
                        }
                    }
                }               
                point.level[dim] = level;
                point.index[dim] = index;   
            }                        
        }
    }
    fn create_gridpoint_internal(&self, storage: &mut SparseGridStorage<D>, mut point: GridPoint<D>)
    {

        if let Some(point) = storage.get_mut(&point)
        {
            point.is_leaf = false;
        }
        else
        {
            point.is_leaf = false;
            self.create_point(storage, point);
        }
    }
    fn create_point_1d_with_boundary(&self, mut point: GridPoint<D>, storage: &mut SparseGridStorage<D>, dim: usize)
    {
        let level = point.level[dim];
        let index = point.index[dim];

        // stuff for boundaries...
        if level == 1 && D > 1
        {
            // check if we need some additional points on the boundaries,
            // only needed on a N dim grid
            
            // test if there are boundaries in every dimension for this grid point
            // left boundary
            point.index[dim] = 0;
            point.level[dim] = 0;
            self.create_gridpoint_internal(storage, point);
    
            // right boundary
            point.level[dim] = 0;
            point.index[dim] = 1;
            self.create_gridpoint_internal(storage, point);
    
            // restore values
            point.level[dim] = level;
            point.index[dim] = index;
            
        }        
        self.create_point_1d(point, storage, dim);        
    }
    fn create_point_1d(&self, mut point: GridPoint<D>, storage: &mut SparseGridStorage<D>, dim: usize)
    {
        let level = point.level[dim];
        let index = point.index[dim];
        if level > 1
        {
            if((index + 1) / 2) % 2 == 1
            {
                point.index[dim] = (index + 1) / 2;
                point.level[dim] = level - 1;
            }
            else 
            {
                point.index[dim] = (index - 1) / 2;
                point.level[dim] = level - 1;
            }  
            self.create_gridpoint_internal(storage, point);
            // reset level and index back to original values
            point.level[dim] = level;
            point.index[dim] = index;
        }
    }

    fn refine_gridpoint(&self, storage: &mut SparseGridStorage<D>, index: usize) 
    {
        let point = storage[index];
        storage[index].is_leaf = false;
        for d in 0..D
        {
            self.refine_1d(storage, point, d);
        }
    }
}

impl<const D: usize, const DIM_OUT: usize> SparseGridRefinement<D, DIM_OUT> for BaseRefinement<D, DIM_OUT>
{
    fn refine(&self, storage: &mut SparseGridStorage<D>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]], functor: &dyn RefinementFunctor<D, DIM_OUT>) -> Vec<usize> {
        let mut refinable_nodes: Vec<(usize, f64)> = Vec::new();
        let original_number = storage.len();
        iterate_refinable_points(storage,&mut |(seq, _point)|
        {            
            refinable_nodes.push((seq, functor.eval(storage, alpha, values, seq)));
        });
        refinable_nodes.sort_by(|a, b|a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)));
        for (seq, value) in refinable_nodes
        {
            if value >= functor.threshold()
            {
                self.refine_gridpoint(storage,seq);
            }
        }
        (original_number..storage.len()).collect()
    }

    fn get_num_refinable_points(&self, storage: &SparseGridStorage<D>) -> usize
    {        
        let mut count = 0;
        iterate_refinable_points(storage, &mut |_point| { count += 1; });        
        count
    }

    fn refine_1d(&self, storage: &mut SparseGridStorage<D>, mut point: GridPoint<D>, dim: usize) 
    {
        let index = point.index[dim];
        let level = point.level[dim];
        if level == 0
        {
            point.level[dim] = 1;
            point.index[dim] = 1;
            if !storage.contains(&point)
            {
                point.is_leaf = true;
                self.create_point(storage, point);
            }
        }
        else
        {
            point.level[dim] = level + 1;
            point.index[dim] = 2 * index - 1;
            if !storage.contains(&point)
            {
                point.is_leaf = true;
                self.create_point(storage, point);
            }
            point.level[dim] = level + 1;
            point.index[dim] = index * 2 + 1;
            if !storage.contains(&point)
            {
                point.is_leaf = true;
                self.create_point(storage, point);
            }
        }
    }
}
