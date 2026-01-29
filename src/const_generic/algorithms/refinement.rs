use crate::const_generic::storage::{GridPoint, PointIterator, SparseGridData};

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefinementMode
{
    #[default]
    Isotropic,
    Anisotropic,
}
#[derive(Default, Debug, Clone)]
pub struct RefinementOptions
{
    pub threshold: f64,    
    pub refinement_mode: RefinementMode,
    pub level_limits: Option<Vec<u8>>, 
}

impl RefinementOptions
{
    pub fn new(threshold: f64) -> Self
    {
        Self { threshold, ..Default::default() }
    }
}

///
/// This trait defines operations used for refinement or coarsening. These
/// two operations are never done simulataneously, but provide a common
/// interface to allow user-specified constraints to control either operation.
///
pub trait RefinementFunctor<const D: usize, const DIM_OUT: usize> : Send + Sync
{
    ///
    /// Return criteria for determining refinement threshold
    /// `alpha` represents the surplus coefficients for each point
    /// `values` represents the values at each point
    /// returns the error estimate at each node. A common choice is
    /// to just use the absolute value of the surplus.
    ///
    fn eval(&self, points: PointIterator<D>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]]) -> Vec<f64>;

    ///
    /// Return per-dimension error estimates for anisotropic refinement
    /// `alpha` represents the surplus coefficients for each point
    /// `values` represents the values at each point
    /// returns a vector where each element is a Vec<f64> of length D
    /// containing the error estimate for each dimension at that node.
    /// Default implementation returns uniform errors (same as eval() for all dimensions)
    ///
    fn eval_per_dimension(&self, points: PointIterator<D>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]]) -> Vec<Vec<f64>>
    {
        let error_estimate = self.eval(points, alpha, values);
        error_estimate.iter().map(|&err| vec![err; D]).collect()
    }

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


fn iterate_refinable_points<const D: usize, Op: FnMut((usize, &GridPoint<D>))>(storage: &SparseGridData<D>, operation: &mut Op, level_limits: Option<Vec<u8>>)
{
    let level_limits = if let Some(level_limits) = level_limits.as_ref()
    {
        if level_limits.len() != D
        {
            panic!("level limits must have length equal to D");
        }
        std::array::from_fn(|i| level_limits[i])
    }
    else
    {
        [u8::MAX; D]
    };    
    for (seq, point) in storage.nodes().iter().enumerate()
    {
        let parent = point;
        let mut point = *point;
        for d in 0..D
        {               
            if point.level[d] >= level_limits[d]  
            {
                continue; // skip this point, it is too deep
            }
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
    fn create_point(&self, storage: &mut SparseGridData<D>, point: GridPoint<D>)
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
    fn create_gridpoint_level_zero_consistency(&self, storage: &mut SparseGridData<D>, mut point: GridPoint<D>)
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
                        let leaf_l = point.is_leaf();
                        // check the boundary not being evaluated
                        point.index[dim] = if i == 0 { 1 } else { 0 };
                        if !storage.contains(&point)
                        {
                            let leaf_r = point.is_leaf();
                            point.set_is_leaf(leaf_l);
                            self.create_point(storage, point);
                            point.set_is_leaf(leaf_r);
                        }
                        else if let Some(point) = storage.get_mut(&point)
                        {
                            point.set_is_leaf(leaf_l);   
                        }
                    }
                }               
                point.level[dim] = level;
                point.index[dim] = index;   
            }                        
        }
    }
    fn create_gridpoint_internal(&self, storage: &mut SparseGridData<D>, mut point: GridPoint<D>)
    {

        if let Some(point) = storage.get_mut(&point)
        {
            point.set_is_leaf(false);
        }
        else
        {
            point.set_is_leaf(false);
            self.create_point(storage, point);
        }
    }
    fn create_point_1d_with_boundary(&self, mut point: GridPoint<D>, storage: &mut SparseGridData<D>, dim: usize)
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
    fn create_point_1d(&self, mut point: GridPoint<D>, storage: &mut SparseGridData<D>, dim: usize)
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
        }
    }

    fn refine_gridpoint(&self, storage: &mut SparseGridData<D>, index: usize, level_limits: Option<Vec<u8>>, dim_errors: Option<&[f64]>, threshold: f64)
    {
        let level_limits = level_limits.unwrap_or(vec![u8::MAX; D]);
        let point = storage[index];

        storage[index].set_is_leaf(false);
        for dim in 0..D
        {
            if point.level[dim] == level_limits[dim]
            {
                continue; // skip this dimension, it is too deep
            }

            // For anisotropic refinement, only refine dimensions with high error
            if let Some(errors) = dim_errors
            {
                if errors[dim] <= threshold
                {
                    continue; // skip this dimension, error is too low
                }
            }

            self.refine_1d(storage, point, dim);
        }
    }
}

impl<const D: usize, const DIM_OUT: usize> BaseRefinement<D, DIM_OUT>
{
    ///
    /// Refine grid based on criteria computed using functor
    /// 
    pub fn refine(&self, storage: &mut SparseGridData<D>, alpha: &[[f64; DIM_OUT]], values: &[[f64; DIM_OUT]], functor: &dyn RefinementFunctor<D, DIM_OUT>, options: RefinementOptions) -> Vec<usize> {
        let mut refinable_nodes = Vec::new();
        let original_number = storage.len();

        // For anisotropic refinement, use per-dimension errors
        let use_anisotropic = options.refinement_mode == RefinementMode::Anisotropic;

        if use_anisotropic
        {
            let dim_errors = functor.eval_per_dimension(storage.points(), alpha, values);
            iterate_refinable_points(storage, &mut |(seq, _point)|
            {
                // Check if any dimension exceeds threshold
                if dim_errors[seq].iter().any(|&err| err > options.threshold)
                {
                    refinable_nodes.push(seq);
                }
            }, options.level_limits.clone());

            for seq in refinable_nodes
            {
                self.refine_gridpoint(storage, seq, options.level_limits.clone(), Some(&dim_errors[seq]), options.threshold);
            }
        }
        else // Isotropic refinement
        {
            let error_estimate = functor.eval(storage.points(), alpha, values);
            iterate_refinable_points(storage, &mut |(seq, _point)|
            {
                if error_estimate[seq] > options.threshold
                {
                    refinable_nodes.push(seq);
                }
            }, options.level_limits.clone());

            for seq in refinable_nodes
            {
                self.refine_gridpoint(storage, seq, options.level_limits.clone(), None, options.threshold);
            }
        }
        (original_number..storage.len()).collect()
    }

    ///
    /// Returns the number of grid points that can be refined.
    ///      
    pub fn get_num_refinable_points(&self, storage: &SparseGridData<D>, level_limits: Option<Vec<u8>>) -> usize
    {        
        let mut count = 0;
        iterate_refinable_points(storage, &mut |_point| { count += 1; }, level_limits);        
        count
    }

    ///
    /// Refine a grid point along a single direction
    ///
    pub fn refine_1d(&self, storage: &mut SparseGridData<D>, mut point: GridPoint<D>, dim: usize) 
    {
        let index = point.index[dim];
        let level = point.level[dim];
        if level == 0
        {
            point.level[dim] = 1;
            point.index[dim] = 1;
            if !storage.contains(&point)
            {
                point.set_is_leaf(true);
                self.create_point(storage, point);
            }
        }
        else
        {
            point.level[dim] = level + 1;
            point.index[dim] = 2 * index - 1;
            if !storage.contains(&point)
            {
                point.set_is_leaf(true);
                self.create_point(storage, point);
            }
            point.level[dim] = level + 1;
            point.index[dim] = index * 2 + 1;
            if !storage.contains(&point)
            {
                point.set_is_leaf(true);
                self.create_point(storage, point);
            }
        }
    }
}

