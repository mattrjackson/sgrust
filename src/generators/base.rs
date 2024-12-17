use crate::storage::linear_grid::{GridPoint, SparseGridStorage};

pub trait Generator<const D: usize> : Default
{
    ///
    /// Generates a regular sparse grid of level levels, without boundaries.
    /// For details about T, See pages 8-9 of Griebel and Knapek's "Optimized 
    /// Tensor-Product Approximation Spaces".
    /// 
    #[allow(non_snake_case)]
    fn regular(&self, storage: &mut SparseGridStorage<D>, levels: [usize; D], T :Option<f64>);
    ///
    /// Generates a regular sparse grid of level levels, without boundaries
    /// where dimensions are splitted into a groups with only certain number
    /// of dimensions completely connected in a clique.
    /// For details about T, See pages 8-9 of Griebel and Knapek's "Optimized 
    /// Tensor-Product Approximation Spaces".
    /// 
    #[allow(non_snake_case)]
    fn cliques(&self, storage: &mut SparseGridStorage<D>, levels: [usize; D], clique_size: usize, T :Option<f64>);
    ///
    /// Generates a full grid of 2^@level tensors, without boundaries
    /// 
    fn full(&self, storage: &mut SparseGridStorage<D>, level: usize);

    ///
    /// Generates a full grid of level @level, with boundary grid points.
    /// 
    fn full_with_boundaries(&self, storage: &mut SparseGridStorage<D>, level: usize);

    /// Generates a regular sparse grid of level levels, with boundaries.
    /// For details about T, See pages 8-9 of Griebel and Knapek's "Optimized 
    /// Tensor-Product Approximation Spaces".
    #[allow(non_snake_case)]
    fn regular_with_boundaries(&self, storage: &mut SparseGridStorage<D>, levels: [usize; D], boundary_level: Option<usize>, T :Option<f64>);

}

///
/// Generate a regular sparse grid iteratively without grid points on the boundary.
/// 
#[allow(non_snake_case)]
fn regular_generator_iterative<const D: usize>(storage: &mut SparseGridStorage<D>, levels: [usize; D], T :Option<f64>)
{
    let mut point = GridPoint::new([1; D], [1;D], false);
    let t = T.unwrap_or(0.0); // default to zero (sparse grid)
    let n = levels[0] as u32;
    for l in 1..=n
    {
        for i in (1..(1 << l)).step_by(2)
        {
            let is_leaf = l == n;            
            point.level[0] = l;
            point.index[0] = i;
            point.is_leaf = is_leaf;
            storage.insert_point(point);
        }
    }
    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    #[allow(clippy::needless_range_loop)]
    for d in 1..D
    {
        let ngrids = storage.len();
        for g in 0..ngrids
        {
            let mut first = true;
            let mut point: GridPoint<D> = storage[g];
            let level_sum = point.level_sum() - 1;
            let level_max = point.level_max();
            let mut l = 1;
            // TODO: This is the one change I've made from the SG++ implementation - allow non-uniform levels. TBD if this works or I need to revert this section.
            let n = levels[d] as u32;
            while (l + level_sum) as f64 - (t * l.max(level_max) as f64) <= (n + D as u32 - 1) as f64 - (t * n as f64)  && l.max(level_max) <= n
            {
                for i in (1..(1 << l)).step_by(2)
                {
                    let is_leaf = l + level_sum == n + D as u32 - 1;
                    point.level[d] = l;
                    point.index[d] = i;
                    point.is_leaf = is_leaf;
                    if !first
                    {
                      storage.insert_point(point);                      
                    }
                    else
                    {
                        storage.update(point, g);
                        first = false;
                    }
                }
                l += 1;
            }
        }
    }
}

#[allow(non_snake_case)]
pub fn regular<const D: usize>(storage: &mut SparseGridStorage<D>, levels: [usize; D], T :Option<f64>) 
{
    regular_generator_iterative(storage, levels, T)
}

#[allow(non_snake_case)]
pub fn cliques<const D: usize>(storage: &mut SparseGridStorage<D>, levels: [usize; D], clique_size: usize, T: Option<f64>) {
    let mut point = GridPoint::new([1; D], [1;D], false);
    let t = T.unwrap_or(0.0); // default to zero (sparse grid)
    let n = levels[0] as u32;
    for l in 1..=n
    {
        for i in (1..(1 << l)).step_by(2)
        {
            let is_leaf = l == n;            
            point.level[0] = l;
            point.index[0] = i;
            point.is_leaf = is_leaf;
            storage.insert_point(point);
        }
    }
    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    #[allow(clippy::needless_range_loop)]
    for d in 1..D
    {
        let ngrids = storage.len();
        let clique_num = d / clique_size;
        for g in 0..ngrids
        {
            let mut first = true;
            let mut point: GridPoint<D> = storage[g];
            let level_sum = point.level_sum() - 1;
            
            let mut l = 1;
            // TODO: This is the one change I've made from the SG++ implementation - allow non-uniform levels. TBD if this works or I need to revert this section.
            let n = levels[d] as u32;
            let mut dt = 0;
            let mut skip = false;
            while dt < clique_size * clique_num && dt < D 
            {
                if point.level[d] > 1
                {
                    skip = true;
                    break;
                }
                dt += 1;
            }
            if skip
            {
                continue;
            }
            let level_max = point.level_max();
            while (l + level_sum) as f64 - (t * l.max(level_max) as f64) <= (n + D as u32 - 1) as f64 - (t * n as f64)  && l.max(level_max) <= n
            {
                for i in (1..(1 << l)).step_by(2)
                {
                    let is_leaf = l + level_sum == n + D as u32 - 1;
                    point.level[d] = l;
                    point.index[d] = i;
                    point.is_leaf = is_leaf;
                    if !first
                    {
                      storage.insert_point(point);                      
                    }
                    else
                    {
                        storage.update(point, g);
                        first = false;
                    }
                }
                l += 1;
            }
        }
    }
}

fn full_iterative<const D: usize>(storage: &mut SparseGridStorage<D>, level: usize)
{
    let mut point = GridPoint::new([1; D], [1;D], false);   
    let n = level as u32;
    for l in 1..=n
    {
        for i in (1..(1 << l)).step_by(2)
        {
            let is_leaf = l == n;            
            point.level[0] = l;
            point.index[0] = i;
            point.is_leaf = is_leaf;
            storage.insert_point(point);
        }
    }
    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    for d in 1..D
    {
        let ngrids = storage.len();
        for g in 0..ngrids
        {
            let mut first = true;
            let mut point: GridPoint<D> = storage[g];                   
            for l in 1..=n            
            {
                for i in (1..(1 << l)).step_by(2)
                {
                    let is_leaf = point.level_sum() == n * D as u32;
                    point.level[d] = l;
                    point.index[d] = i;
                    point.is_leaf = is_leaf;
                    if !first
                    {
                      storage.insert_point(point);                      
                    }
                    else
                    {
                        storage.update(point, g);
                        first = false;
                    }
                }
            }
        }
    }
}

pub fn full<const D: usize>(storage: &mut SparseGridStorage<D>, level: usize) {
    full_iterative(storage, level);
}

pub fn anisotropic_full<const D: usize>(_storage: &mut SparseGridStorage<D>, _level: &[usize; D]) {
    todo!()
}

pub fn full_with_boundaries<const D: usize>(storage: &mut SparseGridStorage<D>, level: usize) {
    full_with_boundaries_iter(storage, level);
}

#[allow(non_snake_case)]
pub fn regular_with_boundaries<const D: usize>(storage: &mut SparseGridStorage<D>, levels: [usize; D], boundary_level: Option<usize>, T :Option<f64>) {
    let boundary_level = boundary_level.unwrap_or(1);
    if boundary_level > 0
    {
        regular_with_boundaries_iter(storage, levels, Some(boundary_level), T);
    }   
    else
    {
        todo!("Need to implement recursive generator");
    } 
}

#[allow(non_snake_case)]
fn regular_with_boundaries_iter<const D:usize>(storage: &mut SparseGridStorage<D>, levels: [usize; D], boundary_level: Option<usize>, T :Option<f64>)
{
    let boundary_level = boundary_level.unwrap_or(1) as u32;
    let mut point = GridPoint::new([1; D], [1;D], false);
    let t = T.unwrap_or(0.0); // default to zero (sparse grid)
    let n = levels[0] as u32;

    point.level[0] = 0;
    point.index[0] = 0;
    storage.insert_point(point);    
    point.index[0] = 1;
    storage.insert_point(point);
    
    for l in 1..=n
    {
        for i in (1..(1 << l)).step_by(2)
        {
            let is_leaf = l == n;            
            point.level[0] = l;
            point.index[0] = i;
            point.is_leaf = is_leaf;
            storage.insert_point(point);
        }
    }
    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    for d in 1..D as u32
    {
        let ngrids = storage.len();
        let cur_dim = d+1;
        for g in 0..ngrids
        {
            let mut point: GridPoint<D> = storage[g];
            let mut level_sum = 0;
            let mut num_zero_levels = 0;
            
            for j in 0..d as usize
            {  
                let lvl = point.level[j];
                if lvl == 0
                {
                    num_zero_levels += 1;
                }
                level_sum += lvl;
            }

            // generate boundary basis functions,
            // but only if levelSum <=
            // n + curDim - boundaryLevel - (numberOfZeroLevels + 1)
            // (the +1 comes from the fact that the newly generated functions
            // will have an additional zero in the d-th dimension)
            let mut first_point = true;
            if (level_sum + boundary_level + num_zero_levels < n + cur_dim) || (num_zero_levels == cur_dim - 1)
            {
                point.level[d as usize] = 0;
                point.index[d as usize] = 0;
                point.is_leaf = false;
                storage.update(point, g);
                point.index[d as usize] = 1;
                storage.insert_point(point);
                first_point = false;
            }
            // choose upper bound of level sum according whether
            // the new basis function is an interior or a boundary function
            // (the loop below skips l = 0 as the boundary points
            // have been inserted a few lines above)
            let mut upper_bound = if num_zero_levels > 0
            {
                // check if upperBound would be negative
                // (we're working with unsigned integers here)
                if n + cur_dim  < boundary_level + num_zero_levels 
                {
                  continue;
                } else {
                  // upper bound for boundary basis functions
                  (n + cur_dim - num_zero_levels - boundary_level) as f64
                }
            } 
            else 
            {
                // upper bound for interior basis functions
                (n + cur_dim - 1) as f64
            };
            upper_bound -= t * n as f64;
            let level_max = point.level_max();
            let mut l = 1;
            let d = d as usize;
            // TODO: This is the one change I've made from the SG++ implementation - allow non-uniform levels. TBD if this works or I need to revert this section.
            let n = levels[d] as u32;
            while (l + level_sum) as f64 - (t * l.max(level_max) as f64) <= upper_bound && l.max(level_max) <= n
            {
                for i in (1..(1 << l)).step_by(2)
                {
                    let is_leaf = if l+level_sum == n + D as u32 - 1 { num_zero_levels == 0 } else { false};
                    point.level[d] = l;
                    point.index[d] = i;
                    point.is_leaf = is_leaf;
                    if !first_point
                    {
                      storage.insert_point(point);                      
                    }
                    else
                    {
                        storage.update(point, g);
                        first_point = false;
                    }
                }
                l += 1;
            }
        }
    }
}

fn full_with_boundaries_iter<const D:usize>(storage: &mut SparseGridStorage<D>, level: usize)
{
    let mut point = GridPoint::new([1; D], [1;D], false);
    let n = level as u32;

    // point.level[0] = 0;
    // point.index[0] = 0;
    // storage.insert_point(point);    
    // point.index[0] = 1;
    // storage.insert_point(point);
    
    for l in 1..=n
    {
        if l == 1
        {
            point.level[0] = 0;
            point.index[0] = 0;
            point.is_leaf = false;
            storage.insert_point(point);

            point.index[0] = 1;
            storage.insert_point(point);
        }
        for i in (1..(1 << l)).step_by(2)
        {
            let is_leaf = l == n;            
            point.level[0] = l;
            point.index[0] = i;
            point.is_leaf = is_leaf;
            storage.insert_point(point);
        }
    }
    
    // Generate grid points in all other dimensions:
    // loop dim times over intermediate grid, take all grid points and
    // modify them in current dimension d
    for d in 1..D as u32
    {
        let ngrids = storage.len();
        for g in 0..ngrids
        {
            let mut point: GridPoint<D> = storage[g];
      
            
            for l in 1..=n
            {
                if l == 1
                {
                    // handle level 1
                    point.level[d as usize] = 0;
                    point.index[d as usize] = 0;
                    point.is_leaf = false;
                    storage.update(point, g);
                    point.index[d as usize] = 1;
                    storage.insert_point(point);
                }
                point.level[d as usize] = l;
                for i in (1..(1 <<l)).step_by(2)
                {
                    point.index[d as usize] = i;
                    point.is_leaf = point.level_sum() == n * D as u32;
                    storage.insert_point(point);
                }
            }   
        }
    }
}

#[test]
fn test_regular()
{
    let mut storage = SparseGridStorage::<2>::default();   
    regular(&mut storage, [3,3], Some(0.0));
    assert_eq!(storage.len(), 17);
}
#[test]
fn test_truncated_boundaries_1d()
{
    let mut storage = SparseGridStorage::<1>::default();
    regular_with_boundaries(&mut storage, [2], Some(1), None);
    assert_eq!(storage.len(), 5);
}
#[test]
fn test_truncated_boundaries_2d()
{
    let mut storage = SparseGridStorage::<2>::default();
    regular_with_boundaries(&mut storage, [2,2], Some(1), None);
    assert_eq!(storage.len(), 21);
    let mut storage2 = SparseGridStorage::<2>::default();
    regular_with_boundaries(&mut storage2, [3,3], Some(1), None);
    assert_eq!(storage2.len(), 49);   
    assert!(storage2.contains(&GridPoint::new([1,1], [1,1], false)));
    assert!(storage2.contains(&GridPoint::new([1,2], [1,1], false)));
    assert!(storage2.contains(&GridPoint::new([2,2], [3,1], false)));
    assert!(!storage2.contains(&GridPoint::new([3,2], [5,1], false)));
    assert!(storage2.contains(&GridPoint::new([3,1], [5,1], false)));
    assert!(storage2.contains(&GridPoint::new([3,0], [5,0], false)));
    assert!(storage2.contains(&GridPoint::new([0,0], [0,0], false)));
}