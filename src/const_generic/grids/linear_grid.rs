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
    fn regular(&self, storage: &mut SparseGridData<D>, levels: [usize; D], T :Option<f64>) -> Result<(), SGError> {
        regular(storage, levels, T)
    }

    #[allow(non_snake_case)]
    fn cliques(&self, storage: &mut SparseGridData<D>, levels: [usize; D], clique_size: usize, T :Option<f64>) -> Result<(), SGError> {
        cliques(storage, levels, clique_size, T)
    }

    fn full(&self, storage: &mut SparseGridData<D>, level: usize) -> Result<(), SGError>{
       full(storage, level)
    }
    fn full_with_boundaries(&self, storage: &mut SparseGridData<D>, level: usize) -> Result<(), SGError> {
        full_with_boundaries(storage, level)
    }

    #[allow(non_snake_case)]
    fn regular_with_boundaries(&self, storage: &mut SparseGridData<D>, levels: [usize; D], boundary_level: Option<usize>, T :Option<f64>) -> Result<(), SGError> {
        regular_with_boundaries(storage, levels, boundary_level, T)
    }
}

#[serde_as]
#[derive(Default, Serialize, Deserialize, Clone)]
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
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
    
    pub fn hierarchize(&mut self) -> Result<(), SGError>{
        if !self.0.has_boundary()
        {
            let op = LinearHierarchisationOperation;  
            self.0.hierarchize(&op)
        }
        else
        {
            let op = LinearBoundaryHierarchisationOperation;
            self.0.hierarchize(&op)  
        }
    }
    pub fn refine<F: RefinementFunctor<D, DIM_OUT>, EF: Fn(&[f64;D])->[f64; DIM_OUT]>(&mut self, functor: &F, eval_fun: &EF, options: RefinementOptions, max_iterations: usize) -> Result<(), SGError>
    {
        if !self.0.has_boundary()
        {
            let op = LinearHierarchisationOperation;  
            self.0.refine(functor, eval_fun, &op, options, max_iterations)?;
        }
        else
        {
            let op = LinearBoundaryHierarchisationOperation;
            self.0.refine(functor, eval_fun, &op, options, max_iterations)?;
        }        
        self.base_mut().storage.generate_adjacency_data();
        Ok(())
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

    pub fn update_refined_values(&mut self, values: &[[f64; DIM_OUT]], sort_data: bool) -> Result<(), SGError>
    {
        let starting_index = self.0.values().len() - values.len();
        for (value, &new_value) in self.base_mut().values[starting_index..].iter_mut().zip(values)
        {
            *value = new_value;
        }
        self.hierarchize()?;
        if sort_data
        {            
            self.sort();
        }    
        Ok(())
    }
    pub fn sort(&mut self) {
        self.0.sort();
        self.0.storage.generate_adjacency_data();
    }

    pub fn sparse_grid(&mut self, levels: [usize; D]) -> Result<(), SGError> {
        self.0.sparse_grid(levels, &LinearGridGenerator)
    }

    pub fn full_grid(&mut self, level: usize) -> Result<(), SGError> {
        self.0.full_grid(level, &LinearGridGenerator)
    }

    pub fn sparse_grid_with_boundaries(&mut self, levels: [usize; D]) -> Result<(), SGError>{
        self.0.sparse_grid_with_boundaries(levels, &LinearGridGenerator)
    }

    pub fn full_grid_with_boundaries(&mut self, level: usize) -> Result<(), SGError>{
        self.0.full_grid_with_boundaries(level, &LinearGridGenerator)
    }
    pub fn integrate(&self) -> [f64; DIM_OUT]
    {
        self.0.integrate_isotropic()
    }
    
    pub fn read<Reader: std::io::Read>(reader: Reader, format: crate::serialization::SerializationFormat) -> Result<Self, SGError> where Self: Sized {
        Ok(Self(SparseGridBase::<D, DIM_OUT>::read(reader, format)?))
    }
    
    pub fn read_buffer(buffer: &[u8], format: crate::serialization::SerializationFormat) -> Result<Self, SGError> where Self: Sized {
        Ok(Self(SparseGridBase::<D, DIM_OUT>::read_buffer(buffer, format)?))
    }    

    pub fn write_buffer(&self, format: crate::serialization::SerializationFormat) -> Result<Vec<u8>, SGError> where Self: Sized {
        self.base().write_buffer(format)
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
    pub fn points(&self) -> PointIterator<'_, D>
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
        self.hierarchize()?;
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
    
    /// Write data to path
    pub fn write(&mut self, path: &str, format: crate::serialization::SerializationFormat) -> Result<(), SGError>
    {
        self.0.write(path, format)
    }
    
}

#[test]
fn check_make_grid_1d()
{
    use crate::const_generic::storage::BoundingBox;
    let level = 8;
    let mut grid = LinearGrid::<1,1>::new();
    *grid.bounding_box_mut() = BoundingBox::new([0.0], [1.00]);
    grid.full_grid_with_boundaries(level).expect("Could not create grid.");
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
    grid.full_grid(level).expect("Could not create grid.");
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
    grid.full_grid_with_boundaries(level).expect("Could not create grid.");
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();
    grid.hierarchize().expect("Could not hierarchize grid");
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
    grid.sparse_grid_with_boundaries([12,12]).expect("Could not create grid.");

    let points: Vec<[f64; 2]> = grid.points().collect();
    let mut values = vec![[0.0; 1]; points.len()];
    for (point, value) in points.iter().zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();
    grid.hierarchize().expect("Could not generate grid.");
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
    grid.full_grid_with_boundaries(level).expect("Could not create grid.");
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
    let level = 4;
    let mut grid = LinearGrid::<6,1>::new();
    grid.sparse_grid_with_boundaries([level; 6]).expect("Could not create grid.");
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
    grid.full_grid_with_boundaries(level).expect("Could not create grid.");
   // assert_eq!(grid.storage.len(), (2_i32.pow(level as u32)-1).pow(2) as usize);
    grid.update_values(&|point| [point[0]*point[0] + point[1]]);
    println!("---- After Refinement ----");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    grid.refine(&functor, &|x| [x[0]*x[0] + x[1]], thresholds, 10).expect("Could not refine grid.");
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
    grid.full_grid_with_boundaries(level).expect("Could not create grid.");
   // assert_eq!(grid.storage.len(), (2_i32.pow(level as u32)-1).pow(2) as usize);
    grid.update_values(&|point| [point[0]*point[0] + point[1]]);
    println!("---- After Refinement ----");
    let functor = crate::const_generic::refinement::surplus::SurplusRefinement;
    let mut options = RefinementOptions::new(1e-7);
    options.refinement_mode = crate::const_generic::algorithms::refinement::RefinementMode::Anisotropic;
    //options.level_limits = Some(vec![12, 11]);
    grid.refine(&functor, &|x| [x[0]*x[0] + x[1]], options, 15).expect("Could not refine grid.");
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
    grid.full_grid_with_boundaries(level).expect("Could not create grid.");
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
    grid.full_grid_with_boundaries(level).expect("Could not create grid.");
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
        grid.update_refined_values(&values, false).expect("Could not update refined values.");
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
    grid.full_grid_with_boundaries(level).expect("Could not create grid.");
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
        grid.update_refined_values(&values, false).expect("Could not update refined values.");
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
    grid.full_grid_with_boundaries(level).expect("Could not create grid.");
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
    grid.sparse_grid_with_boundaries([level]).expect("Could not create grid.");
    let points = grid.points();
    let mut values = vec![[0.0; 1]; grid.len()];
    for (point, value) in points.zip(values.iter_mut())
    {
        let x = point[0];
        value[0] = 0.5*(1.0+erf((x-0.5)/sigma/SQRT_2));
    }
    grid.set_values(values.clone()).unwrap();    
    grid.hierarchize().expect("Could not hierarchize grid.");
    let functor = SurplusRefinement;
    grid.refine(&functor, &|x| [0.5*(1.0+erf((x[0]-0.5)/sigma/SQRT_2)); 1], thresholds, 20).expect("Could not refine grid.");
    let x = [0.323];
    let exact = 0.5*(1.0+erf((x[0]-0.5)/sigma/SQRT_2));
    assert!((grid.interpolate(x).unwrap()[0]-exact).abs() < 1e-3);
}


#[test]
fn compare_3d_sin_refinement_isotropic_vs_anisotropic()
{
    use crate::const_generic::refinement::surplus::SurplusRefinement;

    println!("\n=== Comparing 3D Sin Function Refinement ===\n");
    
    // Test function: f(x,y,z) = sin(2π*x) * sin(2π*y) * sin(π*z)
    // High frequency in x,y; lower frequency in z
    let eval_fn = |point: &[f64;3]| -> [f64;1] {
        [1.0/(0.5 - point[0].powi(4)- point[1].powi(4)).abs() + 0.1]
    };
    
    // Create two identical grids for isotropic and anisotropic refinement
    let mut grid_iso = LinearGrid::new();
    let mut grid_aniso = LinearGrid::new();
    
    // Initialize both grids with sparse grid at level [1,1,1] (smaller initial grid)
    grid_iso.full_grid_with_boundaries(2).expect("failed to create isotropic grid");
    grid_iso.update_values(&eval_fn);
    
    grid_aniso.full_grid_with_boundaries(2).expect("failed to create anisotropic grid");
    grid_aniso.update_values(&eval_fn);
    let initial_count = grid_iso.len();
    println!("Initial grid size: {} points\n", initial_count);
    
    // Create refinement options with different thresholds to show different behavior
    let mut options_iso = RefinementOptions::new(0.01);  // Refine everything (isotropic)
    options_iso.refinement_mode = crate::const_generic::algorithms::refinement::RefinementMode::Isotropic;
    
    let mut options_aniso = RefinementOptions::new(0.01);  // Refine everything (anisotropic)
    options_aniso.refinement_mode = crate::const_generic::algorithms::refinement::RefinementMode::Anisotropic;
    
    // Apply refinement with multiple iterations
    let functor = SurplusRefinement;
    grid_iso.refine(&functor, &eval_fn, options_iso, 6)
        .expect("isotropic refinement failed");
    grid_iso.coarsen(&functor, 0.01);
    grid_aniso.refine(&functor, &eval_fn, options_aniso, 6).unwrap();
    
    grid_aniso.coarsen(&functor, 0.01);
    // Get results
    let iso_count = grid_iso.len();
    let aniso_count = grid_aniso.len();
    
    println!("Refinement Results:");
    println!("  Isotropic:   {} points ({:+} new)", iso_count, iso_count - initial_count);
    println!("  Anisotropic: {} points ({:+} new)", aniso_count, aniso_count - initial_count);
    
    // Compute peak approximation error at test points
    let test_points: Vec<[f64; 3]> = (0..8)
        .flat_map(|i| {
            (0..8).map(move |j| {
                (0..8).map(move |k| {
                    [
                        (i as f64 + 0.5) / 8.0,
                        (j as f64 + 0.5) / 8.0,
                        (k as f64 + 0.5) / 8.0,
                    ]
                })
            })
        })
        .flatten()
        .collect();
    
    let mut iso_peak_error: f64 = 0.0;
    let mut aniso_peak_error: f64 = 0.0;
    
    for point in &test_points {
        let exact = eval_fn(point)[0];

        // Isotropic error
        if let Ok(iso_result) = grid_iso.interpolate(*point) {
            let error = (iso_result[0] - exact).abs();
            iso_peak_error = iso_peak_error.max(error);
        }

        // Anisotropic error
        if let Ok(aniso_result) = grid_aniso.interpolate(*point) {
            let error = (aniso_result[0] - exact).abs();
            aniso_peak_error = aniso_peak_error.max(error);
        }
    }
    
    println!("\nApproximation Errors:");
    println!("  Isotropic peak error:   {:.6e}", iso_peak_error);
    println!("  Anisotropic peak error: {:.6e}", aniso_peak_error);
    
    // Efficiency metrics
    let iso_eff = iso_peak_error / (iso_count as f64);
    let aniso_eff = aniso_peak_error / (aniso_count as f64);
    
    println!("\nEfficiency (error per point):");
    println!("  Isotropic:   {:.6e}", iso_eff);
    println!("  Anisotropic: {:.6e}", aniso_eff);
    
    // Analysis
    println!("\nAnalysis:");
    println!("  Test function has high frequency (sin(2π·x), sin(2π·y)) and low frequency (sin(π·z))");
    println!("  Isotropic refines all dimensions equally");
    println!("  Anisotropic refines selectively (more in under-refined dims)");
    println!("  Note: Dynamic implementation currently uses same refinement for both modes");
    println!("        See const_generic version for fully optimized anisotropic strategy");
    
    if iso_peak_error < aniso_peak_error {
        println!("\n  ✓ Isotropic achieves better absolute accuracy ({:.2}% lower error)",
                 ((aniso_peak_error - iso_peak_error) / aniso_peak_error) * 100.0);
    } else if aniso_peak_error < iso_peak_error {
        println!("\n  ✓ Anisotropic achieves better absolute accuracy ({:.2}% lower error)",
                 ((iso_peak_error - aniso_peak_error) / iso_peak_error) * 100.0);
    } else {
        println!("\n  ≈ Both modes achieve similar absolute accuracy");
    }
    
    if iso_eff < aniso_eff {
        println!("  ✓ Isotropic is more efficient ({:.2}% better error per point)",
                 ((aniso_eff - iso_eff) / aniso_eff) * 100.0);
    } else if aniso_eff < iso_eff {
        println!("  ✓ Anisotropic is more efficient ({:.2}% better error per point)",
                 ((iso_eff - aniso_eff) / iso_eff) * 100.0);
    } else {
        println!("  ≈ Both modes have similar efficiency");
    }
    
    // Verify both methods refine successfully
    assert!(iso_count > initial_count, "Isotropic should add points");
    assert!(aniso_count > initial_count, "Anisotropic should add points");
    
    println!("\n✓ 3D sin function refinement comparison complete");
    let start = std::time::Instant::now();
    for _ in 0..1e5 as usize
    {
        let _r = grid_aniso.interpolate([0.3,0.3,0.3]).unwrap();        
    }
    println!("1e5 iterations in {} msec", std::time::Instant::now().duration_since(start).as_millis());
}