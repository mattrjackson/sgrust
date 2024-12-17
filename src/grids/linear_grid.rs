use serde_with::serde_as;
use crate::algorithms::refinement::RefinementFunctor;
use crate::basis::linear::LinearBasis;
use crate::errors::SGError;
use crate::hierarchisation::{LinearBoundaryHierarchisationOperation, LinearHierarchisationOperation};
use crate::storage::linear_grid::SparseGridStorage;
use crate::generators::base::*;
use serde::{Serialize, Deserialize};

use super::sparse_grid::{SparseGrid, SparseGridBase};

#[derive(Default)]
pub struct LinearGridGenerator<const D: usize>;

impl<const D: usize> Generator<D> for LinearGridGenerator<D>
{
    #[allow(non_snake_case)]
    fn regular(&self, storage: &mut SparseGridStorage<D>, levels: [usize; D], T :Option<f64>) {
        regular(storage, levels, T);
    }

    #[allow(non_snake_case)]
    fn cliques(&self, storage: &mut SparseGridStorage<D>, levels: [usize; D], clique_size: usize, T :Option<f64>) {
        cliques(storage, levels, clique_size, T);
    }

    fn full(&self, storage: &mut SparseGridStorage<D>, level: usize) {
       full(storage, level);
    }
    fn full_with_boundaries(&self, storage: &mut SparseGridStorage<D>, level: usize) {
        full_with_boundaries(storage, level);
    }

    #[allow(non_snake_case)]
    fn regular_with_boundaries(&self, storage: &mut SparseGridStorage<D>, levels: [usize; D], boundary_level: Option<usize>, T :Option<f64>) {
        regular_with_boundaries(storage, levels, boundary_level, T);
    }
}

#[serde_as]
#[derive(Default, Serialize, Deserialize, Clone)]
pub struct LinearGrid<const D: usize, const DIM_OUT: usize>(pub(crate) SparseGridBase<D, DIM_OUT>);

impl<const D: usize, const DIM_OUT: usize> SparseGrid<D, DIM_OUT> for  LinearGrid<D, DIM_OUT>
{
    fn new() -> Self {
        Self(SparseGridBase::new())
    }

    fn base(&self) -> &SparseGridBase<D, DIM_OUT> {
        &self.0
    }

    fn base_mut(&mut self) -> &mut SparseGridBase<D, DIM_OUT> {
        &mut self.0
    }    
    
    fn hierarchize(&mut self) {
        if !self.has_boundary()
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
    fn refine<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, eval_fun: &mut dyn FnMut(&[f64;D])->[f64; DIM_OUT], max_iterations: usize)
    {
        if !self.has_boundary()
        {
            let op = LinearHierarchisationOperation;  
            self.0.refine(functor, eval_fun, &op, max_iterations);
        }
        else
        {
            let op = LinearBoundaryHierarchisationOperation;
            self.0.refine(functor, eval_fun, &op, max_iterations);
        }        
        self.update_runtime_data();
    }

    fn sparse_grid(&mut self, levels: [usize; D]) {
        self.0.sparse_grid(levels, &LinearGridGenerator);
    }

    fn full_grid(&mut self, level: usize) {
        self.0.full_grid(level, &LinearGridGenerator);
    }

    fn sparse_grid_with_boundaries(&mut self, levels: [usize; D]) {
        self.0.sparse_grid_with_boundaries(levels, &LinearGridGenerator);
    }

    fn full_grid_with_boundaries(&mut self, level: usize) {
        self.0.full_grid_with_boundaries(level, &LinearGridGenerator);
    }
    fn integrate(&self) -> [f64; DIM_OUT]
    {
        self.base().integrate_isotropic(&LinearBasis)
    }
    fn read_compact<Reader: std::io::Read>(reader: Reader)  -> Result<Self, SGError> where Self: Sized {
        Ok(Self(SparseGridBase::<D, DIM_OUT>::read_compact(reader)?))
    }

    fn read_full<Reader: std::io::Read>(reader: Reader)  -> Result<Self, SGError> where Self: Sized {
        Ok(Self(SparseGridBase::<D, DIM_OUT>::read_full(reader)?))
    }

    fn read_compact_buffer(buffer: &[u8]) -> Result<Self, SGError> where Self: Sized
    {
        Ok(Self(SparseGridBase::<D, DIM_OUT>::read_compact_buffer(buffer)?))
    }

    fn read_full_buffer(buffer: &[u8]) -> Result<Self, SGError> where Self: Sized
    {
        Ok(Self(SparseGridBase::<D, DIM_OUT>::read_full_buffer(buffer)?))
    }
}



#[test]
fn check_make_grid()
{
    let level = 6;
    let mut grid = LinearGrid::<2,1>::new();
    grid.full_grid(level);
    assert_eq!(grid.len(), (2_i32.pow(level as u32)-1).pow(2) as usize);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; points.len()];
    for (point, value) in points.iter().zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();
    assert!((grid.interpolate([0.2,0.2]).unwrap()[0]-0.08).abs() < 1e-4);
    let start = std::time::Instant::now();
    for _i in 0..1e6 as usize
    {
        let _r = grid.interpolate([0.2,0.2]).unwrap();
    }
    println!("number of points={}", grid.len());
    println!("1e6 iterations in {} msec", std::time::Instant::now().duration_since(start).as_millis());
}


#[test]
fn check_make_grid_with_boundaries()
{
    let level = 6;
    let mut grid = LinearGrid::<2,1>::new();
    grid.full_grid_with_boundaries(level);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; points.len()];
    for (point, value) in points.iter().zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();
    grid.hierarchize();
    println!("number of points={}", grid.len());
    println!("coarsening");
    let functor = crate::refinement::surplus::SurplusRefinement(1e-5);
    grid.coarsen(&functor);
    println!("number of points after coarsening={}", grid.len());
    println!("{}", grid.interpolate([0.2,0.2]).unwrap()[0]);
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

    let points: Vec<[f64; 2]> = grid.points();
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

#[test]
fn check_make_grid_with_boundaries_4d()
{
    let level = 5;
    let mut grid = LinearGrid::<4,1>::new();
    grid.full_grid_with_boundaries(level);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; points.len()];
    for (point, value) in points.iter().zip(values.iter_mut())
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


#[test]
fn check_make_grid_with_boundaries_6d()
{
    let level = 8;
    let mut grid = LinearGrid::<6,1>::new();
    grid.sparse_grid_with_boundaries([level; 6]);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; points.len()];
    for (point, value) in points.iter().zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();
    println!("number of points={}", grid.len());
    println!("coarsening");
    let functor = crate::refinement::surplus::SurplusRefinement(1e-4);
    grid.coarsen(&functor);
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
    grid.full_grid_with_boundaries(level);
   // assert_eq!(grid.storage.len(), (2_i32.pow(level as u32)-1).pow(2) as usize);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; points.len()];
    for (point, value) in points.iter().zip(values.iter_mut())
    {
        value[0] = point[0]*point[0] + point[1]*point[1];
    }
    grid.set_values(values.clone()).unwrap();

    for point in grid.storage().iter()
    {
        let c = point.unit_coordinate();
        println!("{},{}", c[0], c[1]);
    }
    println!("---- After Refinement ----");
    let functor = crate::refinement::surplus::SurplusRefinement(1e-6);
    grid.refine(&functor, &mut |x| [x[0]*x[0] + x[1]*x[1]], 20);
    // for point in grid.storage.iter()
    // {
    //     let c = point.unit_coordinate();
    //     println!("{},{}", c[0], c[1]);
    // }    
    println!("number of points={}", grid.len());   
}


#[test]
fn fit_1d_gaussian_cdf()
{
    use crate::refinement::surplus::SurplusRefinement;
    use std::f64::consts::SQRT_2;
    use libm::erf;
    let level = 2;
    let sigma = 0.1;
    let mut grid = LinearGrid::<1, 1>::new();
    grid.sparse_grid_with_boundaries([level]);
    let points = grid.points();
    let mut values = vec![[0.0; 1]; points.len()];
    for (point, value) in points.iter().zip(values.iter_mut())
    {
        let x = point[0];
        value[0] = 0.5*(1.0+erf((x-0.5)/sigma/SQRT_2));
    }
    grid.set_values(values.clone()).unwrap();    
    grid.hierarchize();
    let functor = SurplusRefinement(1e-3);
    grid.refine(&functor, &mut |x| [0.5*(1.0+erf((x[0]-0.5)/sigma/SQRT_2)); 1], 20);
    let x = [0.323];
    let exact = 0.5*(1.0+erf((x[0]-0.5)/sigma/SQRT_2));
    assert!((grid.interpolate(x).unwrap()[0]-exact).abs() < 1e-3);
}