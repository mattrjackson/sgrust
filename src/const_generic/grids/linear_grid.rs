use serde_with::serde_as;
use crate::const_generic::algorithms::refinement::{RefinementFunctor, RefinementOptions};
use crate::errors::SGError;
use crate::const_generic::storage::PointIterator;
use crate::const_generic::algorithms::hierarchisation::{LinearBoundaryHierarchisationOperation, LinearHierarchisationOperation};
use crate::const_generic::storage::{BoundingBox, SparseGridData};
use crate::const_generic::generators::*;
use serde::{Serialize, Deserialize};

use super::sparse_grid::SparseGridBase;

#[derive(Default)]
pub struct LinearGridGenerator<const D: usize>;

impl<const D: usize> Generator<D> for LinearGridGenerator<D>
{
    #[allow(non_snake_case)]
    fn regular(&self, storage: &mut SparseGridData<D>, levels: [usize; D], T :Option<f64>) {
        regular(storage, levels, T);
    }

    #[allow(non_snake_case)]
    fn cliques(&self, storage: &mut SparseGridData<D>, levels: [usize; D], clique_size: usize, T :Option<f64>) {
        cliques(storage, levels, clique_size, T);
    }

    fn full(&self, storage: &mut SparseGridData<D>, level: usize) {
       full(storage, level);
    }
    fn full_with_boundaries(&self, storage: &mut SparseGridData<D>, level: usize) {
        full_with_boundaries(storage, level);
    }

    #[allow(non_snake_case)]
    fn regular_with_boundaries(&self, storage: &mut SparseGridData<D>, levels: [usize; D], boundary_level: Option<usize>, T :Option<f64>) {
        regular_with_boundaries(storage, levels, boundary_level, T);
    }
}

#[serde_as]
#[derive(Default, Serialize, Deserialize, Clone)]
pub struct LinearGrid<const D: usize, const DIM_OUT: usize>(pub(crate) SparseGridBase<D, DIM_OUT>);

impl<const D: usize, const DIM_OUT: usize> LinearGrid<D, DIM_OUT>
{
    pub fn new() -> Self {
        Self(SparseGridBase::new())
    }

    pub fn base(&self) -> &SparseGridBase<D, DIM_OUT> {
        &self.0
    }

    pub fn base_mut(&mut self) -> &mut SparseGridBase<D, DIM_OUT> {
        &mut self.0
    }    
    
    pub fn hierarchize(&mut self) {
        if !self.0.has_boundary()
        {
            let op = LinearHierarchisationOperation;  
            self.0.hierarchize(&op);  
        }
        else
        {
            let op = LinearBoundaryHierarchisationOperation;
            self.0.hierarchize(&op);  
        }
    }
    pub fn refine<F: RefinementFunctor<D, DIM_OUT>, EF: Fn(&[f64;D])->[f64; DIM_OUT]>(&mut self, functor: &F, eval_fun: &EF, options: RefinementOptions, max_iterations: usize)
    {
        if !self.0.has_boundary()
        {
            let op = LinearHierarchisationOperation;  
            self.0.refine(functor, eval_fun, &op, options, max_iterations);
        }
        else
        {
            let op = LinearBoundaryHierarchisationOperation;
            self.0.refine(functor, eval_fun, &op, options, max_iterations);
        }        
        self.base_mut().storage.generate_adjacency_data();
    }

    #[cfg(feature="rayon")]
    pub fn refine_parallel<F: RefinementFunctor<D, DIM_OUT>, EF: Fn(&[f64;D])->[f64; DIM_OUT] + Send + Sync>(&mut self, functor: &F, eval_fun: &EF, options: RefinementOptions, max_iterations: usize)
    {
        if !self.0.has_boundary()
        {
            let op = LinearHierarchisationOperation;  
            self.0.refine_parallel(functor, eval_fun, &op, options, max_iterations);
        }
        else
        {
            let op = LinearBoundaryHierarchisationOperation;
            self.0.refine_parallel(functor, eval_fun, &op, options, max_iterations);
        }        
        self.base_mut().storage.generate_adjacency_data();
    }

    pub fn update_refined_values(&mut self, values: &[[f64; DIM_OUT]], sort_data: bool)
    {
        let starting_index = self.0.values().len() - values.len();
        for (value, &new_value) in self.base_mut().values[starting_index..].iter_mut().zip(values)
        {
            *value = new_value;
        }
        self.hierarchize();
        if sort_data
        {            
            self.sort();
        }    
    }
    pub fn sort(&mut self) {
        self.0.sort();
        self.0.storage.generate_adjacency_data();
    }

    pub fn sparse_grid(&mut self, levels: [usize; D]) {
        self.0.sparse_grid(levels, &LinearGridGenerator);
    }

    pub fn full_grid(&mut self, level: usize) {
        self.0.full_grid(level, &LinearGridGenerator);
    }

    pub fn sparse_grid_with_boundaries(&mut self, levels: [usize; D]) {
        self.0.sparse_grid_with_boundaries(levels, &LinearGridGenerator);
    }

    pub fn full_grid_with_boundaries(&mut self, level: usize) {
        self.0.full_grid_with_boundaries(level, &LinearGridGenerator);
    }
    pub fn integrate(&self) -> [f64; DIM_OUT]
    {
        self.0.integrate_isotropic()
    }
    
    pub fn read<Reader: std::io::Read>(reader: Reader) -> Result<Self, SGError> where Self: Sized {
        Ok(Self(SparseGridBase::<D, DIM_OUT>::read(reader)?))
    }
    
    pub fn read_buffer(buffer: &[u8]) -> Result<Self, SGError> where Self: Sized {
        Ok(Self(SparseGridBase::<D, DIM_OUT>::read_buffer(buffer)?))
    }    

     /// Get the surplus coefficients for this grid.
    pub fn alpha(&self) -> &[[f64; DIM_OUT]]
    {
        &self.0.alpha
    }

    /// Get the surplus coefficients for this grid (mutable).
    pub fn alpha_mut(&mut self) -> &mut Vec<[f64; DIM_OUT]>
    {
        &mut self.base_mut().alpha
    }

    /// Get the bounding box for this grid.
    pub fn bounding_box(&self) -> &BoundingBox<D>
    {
        self.0.bounding_box()
    }

    /// Get the bounding box for this grid (mutable).
    pub fn bounding_box_mut(&mut self) -> &mut BoundingBox<D>
    {
        self.0.bounding_box_mut()
    }

    pub fn is_empty(&self) -> bool
    {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize
    {
        self.0.len()
    }

    /// Check if boundaries are enabled for this grid.
    pub fn has_boundary(&self) -> bool
    {
        self.0.has_boundary()
    }

    /// Retrieve the underlying storage 
    pub fn storage(&self) -> &SparseGridData<D>
    {
        &self.0.storage
    }

    /// Get copy of points that make up this grid.
    pub fn points(&self) -> PointIterator<D>
    {
        self.0.points()
    }

    /// Return reference to values on this grid.
    pub fn values(&self) -> &Vec<[f64;DIM_OUT]>
    {
        &self.0.values
    }

    /// Interpolate on grid (single point). Checks that point lies within bounding box.
    pub fn interpolate(&self, x: [f64; D]) -> Result<[f64; DIM_OUT], SGError>
    {
        self.0.interpolate(x)
    }
    
    #[cfg(feature="rayon")]
    /// Interpolate on grid (multiple point). Checks that each point lies within bounding box.
    pub fn interpolate_batch(&self, x: &[[f64; D]]) -> Vec<Result<[f64; DIM_OUT], SGError>>
    {
        self.0.interpolate_batch(x)
    }

    /// Interpolate on grid (single point). No bounding box check.
    pub fn interpolate_unchecked(&self, x: [f64; D]) -> Result<[f64; DIM_OUT], SGError>
    {
        self.0.interpolate_unchecked(x)
    }
    
    #[cfg(feature="rayon")]
    /// Interpolate on grid (multiple point). No bounding box check.
    pub fn interpolate_batch_unchecked(&self, x: &[[f64; D]]) -> Vec<Result<[f64; DIM_OUT], SGError>>
    {
        self.0.interpolate_batch_unchecked(x)
    }

    /// Set values on this grid.
    pub fn set_values(&mut self, values: Vec<[f64; DIM_OUT]>) -> Result<(), SGError>
    {
        self.base_mut().set_values(values)?;
        self.hierarchize();
        self.base_mut().storage.generate_adjacency_data();
        Ok(())
    }

    /// Set values using a given evaluation function...
    pub fn update_values<F: Fn(&[f64;D])->[f64; DIM_OUT]> (&mut self,  eval_fun: &F)
    {
        let mut values = Vec::with_capacity(self.len());
        for point in self.points()
        {
            values.push(eval_fun(&point));
        }
        self.set_values(values).expect("Failed to set values");
    }

    #[cfg(feature="rayon")]
    /// Set values by using evaluation function in parallel...
    pub fn update_values_parallel<EF: Fn(&[f64;D])->[f64; DIM_OUT] + Send + Sync>(&mut self, eval_fun: &EF)
    {
        use rayon::iter::{ParallelBridge, ParallelIterator};
        let mut values = vec![[0.0; DIM_OUT]; self.len()];
        self.points().zip(values.iter_mut()).par_bridge().for_each(
        |(point, value)|
        {
            *value = eval_fun(&point);
        });        
        self.set_values(values).expect("Failed to set values");
    }
    
    /// Coarsen grid based on functor `F`
    pub fn coarsen<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, threshold: f64) -> usize
    {
        self.0.coarsen(functor, true, threshold)
    }
    
    pub fn refine_iteration<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, options: RefinementOptions) -> Vec<[f64; D]>
    {
        self.0.refine_iteration(functor, options)
    }
    
    /// Save data to path
    pub fn save(&mut self, path: &str) -> Result<(), SGError>
    {
        self.0.save(path)
    }
    
}

#[test]
fn check_make_grid_1d()
{
    use crate::const_generic::storage::BoundingBox;
    let level = 8;
    let mut grid = LinearGrid::<1,1>::new();
    *grid.bounding_box_mut() = BoundingBox::new([0.0], [1.00]);
    grid.full_grid_with_boundaries(level);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        value[0] = point[0]*point[0];        
        //println!("{},{}", point[0], value[0]);
    }    
    grid.set_values(values.clone()).unwrap();
    println!("number of points={}", grid.len());
    println!("interpolated value={}", grid.interpolate([0.2]).unwrap()[0]);
    assert!((grid.interpolate([0.2]).unwrap()[0]-0.04).abs() < 1e-2);
    let start = std::time::Instant::now();
    for _i in 0..1e6 as usize
    {
        let _r = grid.interpolate([0.8]).unwrap();
    }
    
    println!("1e6 iterations in {} msec", std::time::Instant::now().duration_since(start).as_millis());
}



#[test]
fn check_make_grid_2d()
{
    let level = 8;
    let mut grid = LinearGrid::<2,1>::new();
    grid.full_grid(level);
    assert_eq!(grid.len(), (2_i32.pow(level as u32)-1).pow(2) as usize);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }    
    grid.set_values(values.clone()).unwrap();
    println!("interpolated value={}", grid.interpolate([0.2,0.2]).unwrap()[0]);
    assert!((grid.interpolate( [0.2,0.2]).unwrap()[0]-0.08).abs() < 1e-2);
    let start = std::time::Instant::now();
    for _i in 0..1e7 as usize
    {
        let _r = grid.interpolate([0.2,0.2]).unwrap();
    }
    println!("number of points={}", grid.len());
    println!("1e7 iterations in {} msec", std::time::Instant::now().duration_since(start).as_millis());
}

#[cfg(feature="rayon")]
#[test]
fn check_make_grid_with_boundaries_2d()
{
    let level = 6;
    let mut grid = LinearGrid::<2,1>::new();
    grid.full_grid_with_boundaries(level);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    let thresholds = RefinementOptions::new(1e-3);
    for (point, value) in points.zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();
    grid.hierarchize();
    println!("number of points={}", grid.len());
    println!("coarsening");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    grid.coarsen(&functor, 1e-5);
    println!("number of points after coarsening={}", grid.len());
    println!("interpolated_value={}", grid.interpolate([0.2,0.2]).unwrap()[0]);
    assert!((grid.interpolate([0.2,0.2]).unwrap()[0]-0.08).abs() < 1e-4);
    let start = std::time::Instant::now();
    let points = vec![[0.2,0.2]; 1e6 as usize];
    let _ = grid.interpolate_batch(&points);
    println!("number of points={}", grid.len());
    println!("1e6 iterations in {} msec", std::time::Instant::now().duration_since(start).as_millis());
}

#[test]
fn check_integration()
{
    let mut grid: LinearGrid<2, 1> = LinearGrid::<2,1>::new();
    grid.sparse_grid_with_boundaries([12,12]);

    let points: Vec<[f64; 2]> = grid.points().collect();
    let mut values = vec![[0.0; 1]; points.len()];
    for (point, value) in points.iter().zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();
    grid.hierarchize();
    // grid.refine(&SurplusRefinement(1e-6), &mut |point|
    // {
    //     [point[0]*point[0] + point[1]*point[1]]
    // }, 15 );
    println!("number of points={}", grid.len());
    println!("integral={:?}", grid.integrate());
    assert!((1.0-grid.integrate()[0]*3.0/2.0).abs() < 1e-6);
}

#[cfg(feature="rayon")]
#[test]
fn check_make_grid_with_boundaries_4d()
{
    let level = 5;
    let mut grid = LinearGrid::<4,1>::new();
    grid.full_grid_with_boundaries(level);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();
    println!("{}", grid.interpolate([0.2,0.2,0.2,0.2]).unwrap()[0]);
    assert!((grid.interpolate([0.2,0.2,0.2,0.2]).unwrap()[0]-0.08).abs() < 1e-3);
    let start = std::time::Instant::now();
    let points = vec![[0.2,0.2,0.2,0.2]; 1e6 as usize];
    let _ = grid.interpolate_batch(&points);
    println!("number of points={}", grid.len());
    println!("1e6 iterations in {} msec", std::time::Instant::now().duration_since(start).as_millis());
}

#[cfg(feature="rayon")]
#[test]
fn check_make_grid_with_boundaries_6d()
{
    let level = 8;
    let mut grid = LinearGrid::<6,1>::new();
    grid.sparse_grid_with_boundaries([level; 6]);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();
    println!("number of points={}", grid.len());
    println!("coarsening");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    grid.coarsen(&functor, 1e-4);
    println!("number of points after coarsening={}", grid.len());
    println!("{}", grid.interpolate([0.2,0.2,0.0,0.0,0.0,0.0]).unwrap()[0]);
    assert!((grid.interpolate([0.2,0.2,0.3,0.6,0.99,0.99]).unwrap()[0]-0.08).abs() < 1e-2);

    let start = std::time::Instant::now();
    let points= vec![[0.2,0.2,0.3,0.6,0.4,0.7]; 1e6 as usize];    
    let _ = grid.interpolate_batch(&points);
    println!("number of points={}", grid.len());
    println!("1e6 iterations in {} msec", std::time::Instant::now().duration_since(start).as_millis());
}


#[test]
fn check_grid_refinement()
{
    let level = 2;
    let mut grid = LinearGrid::<2,1>::new();
    let thresholds = RefinementOptions::new(1e-7);
    grid.full_grid_with_boundaries(level);
   // assert_eq!(grid.storage.len(), (2_i32.pow(level as u32)-1).pow(2) as usize);
    grid.update_values(&|point| [point[0]*point[0] + point[1]]);
    println!("---- After Refinement ----");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    grid.refine(&functor, &|x| [x[0]*x[0] + x[1]], thresholds, 10);
    // for point in grid.storage.iter()
    // {
    //     let c = point.unit_coordinate();
    //     println!("{},{}", c[0], c[1]);
    // }    
    println!("number of points={}", grid.len());   
    println!("{},{}",grid.interpolate([0.2,0.3]).unwrap()[0], (0.2*0.2+0.3));
    assert!((1.0 - grid.interpolate([0.2,0.3]).unwrap()[0]/(0.2*0.2+0.3)).abs() < 1e-6);
}

#[test]
fn check_grid_refinement_dimension_adaptive()
{
    let level = 2;
    let mut grid = LinearGrid::<2,1>::new();
    grid.full_grid_with_boundaries(level);
   // assert_eq!(grid.storage.len(), (2_i32.pow(level as u32)-1).pow(2) as usize);
    grid.update_values(&|point| [point[0]*point[0] + point[1]]);
    println!("---- After Refinement ----");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    let mut options = RefinementOptions::new(1e-7);
    options.refinement_mode = crate::const_generic::algorithms::refinement::RefinementMode::Anisotropic;
    //options.level_limits = Some(vec![12, 11]);
    grid.refine(&functor, &|x| [x[0]*x[0] + x[1]], options, 15);
    println!("number of points={}", grid.len());   
    println!("{},{}",grid.interpolate([0.2,0.3]).unwrap()[0], (0.2*0.2+0.3));
    assert!((1.0 - grid.interpolate([0.2,0.3]).unwrap()[0]/(0.2*0.2+0.3)).abs() < 1e-6);
}

#[cfg(feature="rayon")]
#[test]
fn check_grid_refinement_parallel()
{
    let level = 2;
    let mut grid = LinearGrid::<2,1>::new();
    let thresholds = RefinementOptions::new(1e-7);
    grid.full_grid_with_boundaries(level);
   // assert_eq!(grid.storage.len(), (2_i32.pow(level as u32)-1).pow(2) as usize);
    grid.update_values_parallel(&|point| [point[0]*point[0] + point[1]*point[1]]);
    println!("---- After Refinement ----");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    grid.refine_parallel(&functor, &|x| [x[0]*x[0] + x[1]*x[1]], thresholds, 20);
    println!("number of points={}", grid.len());   
    println!("{},{}",grid.interpolate([0.2,0.3]).unwrap()[0], (0.2*0.2+0.3*0.3));
    assert!((1.0 - grid.interpolate([0.2,0.3]).unwrap()[0]/(0.2*0.2+0.3*0.3)).abs() < 1e-6);
}

#[test]
fn check_grid_refinement_iteration()
{
    let level = 2;
    let mut grid = LinearGrid::<2,1>::new();
    grid.full_grid_with_boundaries(level);
    let thresholds = RefinementOptions::new(1e-7);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1];
    }
    grid.set_values(values.clone()).unwrap();

    for point in grid.storage().nodes().iter()
    {
        let c = point.unit_coordinate();
        println!("{},{}", c[0], c[1]);
    }
    println!("---- After Refinement ----");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    for _ in 0..20
    {
        let values: Vec<_> = grid.refine_iteration(&functor, thresholds.clone()).iter_mut().map(|point| [point[0]*point[0] + point[1]]).collect();
        grid.update_refined_values(&values, false);
    }
    grid.sort();
    println!("{},{}", grid.interpolate([0.2, 0.3]).unwrap()[0], (0.2*0.2+0.3));
    assert!((1.0 - grid.interpolate([0.2, 0.3]).unwrap()[0] / (0.2*0.2+0.3)).abs() < 1e-6);
    println!("number of points={}", grid.len());   
}

#[test]
fn check_grid_refinement_iteration_dimension_adaptive()
{
    let level = 2;
    let mut grid = LinearGrid::<2,1>::new();
    grid.full_grid_with_boundaries(level);
    let mut options = RefinementOptions::new(1e-7);
    options.refinement_mode = crate::const_generic::algorithms::refinement::RefinementMode::Anisotropic;
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1];
    }
    grid.set_values(values.clone()).unwrap();

    for point in grid.storage().nodes().iter()
    {
        let c = point.unit_coordinate();
        println!("{},{}", c[0], c[1]);
    }
    println!("---- After Refinement ----");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    for _ in 0..20
    {
        let values: Vec<_> = grid.refine_iteration(&functor, options.clone()).iter_mut().map(|point| [point[0]*point[0] + point[1]]).collect();
        grid.update_refined_values(&values, false);
    }
    grid.sort();
    println!("{},{}", grid.interpolate([0.2, 0.3]).unwrap()[0], (0.2*0.2+0.3));
    assert!((1.0 - grid.interpolate([0.2, 0.3]).unwrap()[0] / (0.2*0.2+0.3)).abs() < 1e-6);
    println!("number of points={}", grid.len());   
}


#[cfg(feature="rayon")]
#[test]
fn check_parallel_grid_refinement()
{
    let level = 3;
    let mut grid = LinearGrid::<5,1>::new();
    let thresholds = RefinementOptions::new(1e-6);
    grid.full_grid_with_boundaries(level);
   // assert_eq!(grid.storage.len(), (2_i32.pow(level as u32)-1).pow(2) as usize);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();

    let num_points_original = grid.len();
    println!("---- After Refinement ----");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    grid.refine_parallel(&functor, &|x| [x[0]*x[0] + x[1]*x[1]], thresholds, 20);
    println!("number of points={}", grid.len());   
    assert!(num_points_original < grid.len());
}

#[test]
fn fit_1d_gaussian_cdf()
{
    use crate::const_generic::refinement::surplus::SurplusRefinement;
    use std::f64::consts::SQRT_2;
    use libm::erf;
    let level = 2;
    let sigma = 0.1;
    let thresholds = RefinementOptions::new(1e-3);
    let mut grid = LinearGrid::<1, 1>::new();
    grid.sparse_grid_with_boundaries([level]);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        let x = point[0];
        value[0] = 0.5*(1.0+erf((x-0.5)/sigma/SQRT_2));
    }
    grid.set_values(values.clone()).unwrap();    
    grid.hierarchize();
    let functor = SurplusRefinement;
    grid.refine(&functor, &|x| [0.5*(1.0+erf((x[0]-0.5)/sigma/SQRT_2)); 1], thresholds, 20);
    let x = [0.323];
    let exact = 0.5*(1.0+erf((x[0]-0.5)/sigma/SQRT_2));
    assert!((grid.interpolate(x).unwrap()[0]-exact).abs() < 1e-3);
}