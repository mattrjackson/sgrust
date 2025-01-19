use std::io::Write;

use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelBridge, ParallelIterator};
use serde_with::serde_as;
use crate::algorithms::basis_evaluation::BasisEvaluation;
use crate::algorithms::integration::{AnisotropicQuadrature, IsotropicQuadrature};
use crate::algorithms::refinement::{BaseRefinement, RefinementFunctor, SparseGridRefinement};
use crate::basis::linear::LinearBasis;
use crate::errors::SGError;
use crate::hierarchisation::HierarchisationOperation;
use crate::iterators::grid_iterator_cache::AdjacencyGridIterator;
use crate::storage::linear_grid::{BoundingBox, GridPoint, PointIterator, SparseGridData};
use crate::generators::base::*;
use crate::algorithms;
use serde::{Serialize,Deserialize};

pub trait SparseGrid<const D: usize, const DIM_OUT: usize>
{
    /// Create a new instant of the grid
    fn new() -> Self;

    /// Get the `SparseGridBase` for this grid.
    fn base(&self) -> &SparseGridBase<D, DIM_OUT>;

    /// Get the `SparseGridBase` for this grid (mutable).
    fn base_mut(&mut self) -> &mut SparseGridBase<D, DIM_OUT>;

    /// Get the surplus coefficients for this grid.
    fn alpha(&self) -> &[[f64; DIM_OUT]]
    {
        &self.base().alpha
    }

    /// Get the surplus coefficients for this grid (mutable).
    fn alpha_mut(&mut self) -> &mut Vec<[f64; DIM_OUT]>
    {
        &mut self.base_mut().alpha
    }

    /// Get the bounding box for this grid.
    fn bounding_box(&self) -> &BoundingBox<D>
    {
        self.base().bounding_box()
    }

    /// Get the bounding box for this grid (mutable).
    fn bounding_box_mut(&mut self) -> &mut BoundingBox<D>
    {
        self.base_mut().bounding_box_mut()
    }

    fn is_empty(&self) -> bool
    {
        self.base().is_empty()
    }

    fn len(&self) -> usize
    {
        self.base().len()
    }

    /// Check if boundaries are enabled for this grid.
    fn has_boundary(&self) -> bool
    {
        self.base().has_boundary()
    }

    ///
    /// Generate a regular sparse grid.
    /// 
    fn sparse_grid(&mut self, levels: [usize; D]);

    ///
    /// Generate a full grid.
    /// 
    fn full_grid(&mut self, level: usize);

    /// Generate regular sparse grid with boundaries.
    fn sparse_grid_with_boundaries(&mut self, levels: [usize; D]);
    
    /// Generate full sparse grid with boundaries.
    fn full_grid_with_boundaries(&mut self, level: usize);

    /// Generate grid hierarchy.
    fn hierarchize(&mut self);

    /// Interpolate on grid (single point). Checks that point lies within bounding box.
    fn interpolate(&self, x: [f64; D]) -> Result<[f64; DIM_OUT], SGError>
    {
        self.base().interpolate(x)
    }
    
    /// Interpolate on grid (multiple point). Checks that each point lies within bounding box.
    fn interpolate_batch(&self, x: &[[f64; D]]) -> Vec<Result<[f64; DIM_OUT], SGError>>
    {
        self.base().interpolate_batch(x)
    }

    /// Interpolate on grid (single point). No bounding box check.
    fn interpolate_unchecked(&self, x: [f64; D]) -> Result<[f64; DIM_OUT], SGError>
    {
        self.base().interpolate_unchecked(x)
    }
    
    /// Interpolate on grid (multiple point). No bounding box check.
    fn interpolate_batch_unchecked(&self, x: &[[f64; D]]) -> Vec<Result<[f64; DIM_OUT], SGError>>
    {
        self.base().interpolate_batch_unchecked(x)
    }

    /// Compute integral over grid.
    fn integrate(&self) -> [f64; DIM_OUT];

    /// Retrieve the underlying storage 
    fn storage(&self) -> &SparseGridData<D>
    {
        &self.base().storage
    }

    /// Get copy of points that make up this grid.
    fn points(&self) -> PointIterator<D>
    {
        self.base().points()
    }

    /// Return reference to values on this grid.
    fn values(&self) -> &Vec<[f64;DIM_OUT]>
    {
        &self.base().values
    }

    /// Set values on this grid.
    fn set_values(&mut self, values: Vec<[f64; DIM_OUT]>) -> Result<(), SGError>
    {
        self.base_mut().set_values(values)?;
        self.hierarchize();
        self.base_mut().storage.generate_adjacency_data();
        Ok(())
    }

    /// Set values using a given evaluation function...
    fn update_values(&mut self,  eval_fun: &mut dyn FnMut(&[f64;D])->[f64; DIM_OUT])
    {
        let mut values = Vec::with_capacity(self.len());
        for point in self.points()
        {
            values.push(eval_fun(&point));
        }
        self.set_values(values).expect("Failed to set values");
    }

    /// Set values by using evaluation function in parallel...
    fn update_values_parallel<EF: Fn(&[f64;D])->[f64; DIM_OUT] + Send + Sync>(&mut self, eval_fun: &EF)
    {
        let mut values = vec![[0.0; DIM_OUT]; self.len()];
        self.points().zip(values.iter_mut()).par_bridge().for_each(
        |(point, value)|
        {
            *value = eval_fun(&point);
        });        
        self.set_values(values).expect("Failed to set values");
    }

    /// Coarsen grid based on functor `F`
    fn coarsen<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F) -> usize
    {
        self.base_mut().coarsen(functor, true)
    }

    /// 
    /// Refine grid based on function `F`. Values are filled serially by calling `eval_fun`. 
    /// 
    fn refine<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, eval_fun: &mut dyn FnMut(&[f64;D])->[f64; DIM_OUT], max_iterations: usize);

    /// 
    /// Refine grid based on function `F`, but use `eval_fun` to fill values in parallel.  
    /// 
    fn refine_parallel<F: RefinementFunctor<D, DIM_OUT>, EF: Fn(&[f64;D])->[f64; DIM_OUT] + Send + Sync>(&mut self, functor: &F, eval_fun: &EF, max_iterations: usize);

    /// 
    /// Perform a single refinement iteration. Returns the slice of points that 
    /// 
    fn refine_iteration<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F) -> Vec<[f64; D]>
    {
        self.base_mut().refine_iteration(functor)
    }

    ///
    /// Update refined values for last iteration (used in conjunction with `refine_single_iteration`).
    /// For optimal performance, sort the data after your final iteration. Alternately, `sort` can be 
    /// called directly.
    /// 
    fn update_refined_values(&mut self, values: Vec<[f64; DIM_OUT]>, sort_data: bool);

    
    /// Save data to path
    fn save(&mut self, path: &str) -> Result<(), SGError>
    {
        self.base_mut().save(path)
    }

    /// Sort points... 
    fn sort(&mut self)
    {
        self.base_mut().sort();
    }
    /// Read data
    fn read<Reader: std::io::Read>(reader: Reader) -> Result<Self, SGError> where Self: Sized;    

    /// Read data from a buffer.
    fn read_buffer(buffer: &[u8]) -> Result<Self, SGError> where Self: Sized;

}

#[serde_as]
#[derive(Default, Serialize, Deserialize, Clone)]
pub struct SparseGridBase<const D: usize, const DIM_OUT: usize>
{
    pub(crate) storage: SparseGridData<D>,
    #[serde_as(as = "Vec<[_; DIM_OUT]>")]
    pub(crate) alpha: Vec<[f64; DIM_OUT]>,
    #[serde_as(as = "Vec<[_; DIM_OUT]>")]
    pub(crate) values: Vec<[f64; DIM_OUT]>,    
}


impl<const D: usize, const DIM_OUT: usize> SparseGridBase<D, DIM_OUT>
{
    pub fn new() -> Self
    {
        SparseGridBase { storage: SparseGridData::default(), alpha: Vec::new(), values: Vec::new() }
    }

    pub fn alpha(&self) -> &[[f64; DIM_OUT]]
    {
        &self.alpha
    }

    pub fn alpha_mut(&mut self) -> &mut Vec<[f64; DIM_OUT]>
    {
        &mut self.alpha
    }

    pub fn bounding_box(&self) -> &BoundingBox<D>
    {        
        self.storage.bounding_box()
    }

    pub fn bounding_box_mut(&mut self) -> &mut BoundingBox<D>
    {        
        self.storage.bounding_box_mut()
    }

    pub fn is_empty(&self) -> bool
    {
        self.storage.is_empty()
    }

    pub fn len(&self) -> usize
    {
        self.storage.len()
    }

    pub fn has_boundary(&self) -> bool
    {
        self.storage.has_boundary()
    }

    pub fn sparse_grid<GENERATOR: Generator<D>>(&mut self, levels: [usize; D], generator: &GENERATOR)
    {        
        generator.regular(&mut self.storage, levels, None);       
        self.sort();
    }    

    pub fn full_grid<GENERATOR: Generator<D>>(&mut self, level: usize, generator: &GENERATOR)
    {
        generator.full(&mut self.storage, level);
        self.sort();
    }

    pub fn sparse_grid_with_boundaries<GENERATOR: Generator<D>>(&mut self, levels: [usize; D], generator: &GENERATOR)
    {
        generator.regular_with_boundaries(&mut self.storage, levels, Some(1), None);
        self.storage.has_boundary = true;
        self.sort();
    }    

    pub fn full_grid_with_boundaries<GENERATOR: Generator<D>>(&mut self, level: usize, generator: &GENERATOR)
    {
        generator.full_with_boundaries(&mut self.storage, level);
        self.storage.has_boundary = true;
        self.sort();
    }
    pub fn hierarchize<OP: HierarchisationOperation<D, DIM_OUT>>(&mut self, op: &OP)
    {
        self.alpha.clone_from(&self.values);
        op.hierarchize(&mut self.alpha, &self.storage);        
    }

    pub fn argsort<T: Ord>(data: &[T]) -> Vec<usize> {
        let mut indices = (0..data.len()).collect::<Vec<_>>();
        indices.sort_by_key(|&i| &data[i]);
        indices
    }
    pub(crate) fn sort(&mut self)
    {   
        let mut indices = Self::argsort(self.storage.nodes());        
        let mut indices_rev = indices.clone();
        for i in 0..indices.len()
        {
            indices_rev[indices[i]] = i;
        }
        for i in 0..indices.len()
        {
            let mut current = i;
            while indices[current] != i {
                let next = indices[current];
                if self.alpha.len() == indices.len()
                {
                    self.alpha.swap(current, next);
                }
                if self.values.len() == indices.len()
                {
                    self.values.swap(current, next);
                }
                self.storage.nodes_mut().swap(current, next);
                
                indices[current] = current;                
                current = next;                
            }
            indices[current] = current;
        }
        for value in self.storage.map.values_mut()
        {            
            *value = indices_rev[*value];
        }
    }

    #[inline]
    pub fn interpolate(&self, x: [f64; D]) -> Result<[f64; DIM_OUT], SGError>
    {    
        if !self.bounding_box().contains(&x)
        {
            Err(SGError::OutOfDomain)
        }
        else 
        {
            self.interpolate_unchecked(x)    
        }
    }
    #[inline]
    pub fn interpolate_unchecked(&self, x: [f64; D]) -> Result<[f64; DIM_OUT], SGError>
    {
        use crate::algorithms::interpolation::InterpolationOperation;
        let iterator = &mut AdjacencyGridIterator::new(&self.storage);
        let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, [LinearBasis; D]));      
        op.interpolate(x, &self.alpha, iterator)       
    }

    #[inline]
    pub fn interpolate_batch(&self, x: &[[f64; D]]) -> Vec<Result<[f64; DIM_OUT], SGError>>
    {
       use crate::algorithms::interpolation::InterpolationOperation;
       let mut results = vec![Ok([0.0; DIM_OUT]); x.len()];
        x.par_iter().zip(results.par_iter_mut()).for_each(
            |(x, y)|
            {
                let iterator = &mut AdjacencyGridIterator::new(&self.storage);
                let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, [LinearBasis; D]));
                if !self.bounding_box().contains(x)
                {
                    *y = Err(SGError::OutOfDomain);
                }
                else
                {
                    *y = op.interpolate(*x, &self.alpha, iterator);
                }               
            }
        );       
       results
    }

    #[inline]
    pub fn interpolate_batch_unchecked(&self, x: &[[f64; D]]) -> Vec<Result<[f64; DIM_OUT], SGError>>
    {
        use crate::algorithms::interpolation::InterpolationOperation;
       
       let mut results = vec![Ok([0.0; DIM_OUT]); x.len()];
        x.par_iter().zip(results.par_iter_mut()).for_each(
            |(x, y)|
            {
                let iterator = &mut AdjacencyGridIterator::new(&self.storage);
                let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, [LinearBasis; D]));                
                *y = op.interpolate(*x, &self.alpha, iterator);
            }
        );       
       results
    }

    pub fn integrate_isotropic<QUADRATURE: IsotropicQuadrature<f64, D, DIM_OUT>>(&self, op: &QUADRATURE) -> [f64; DIM_OUT]
    {
        op.eval(self.storage(), &self.alpha)
    }

    pub fn eval_anisotropic_quadrature<QUADRATURE: AnisotropicQuadrature<D, DIM_OUT> + ?Sized>(&self, op: &QUADRATURE, index: usize, dim: usize) -> [f64; DIM_OUT]
    {
        op.eval(self.storage(), index, dim)
    }

    pub fn to_real_coordinate(&self, index: &GridPoint<D>) -> [f64; D]
    {
        let mut point = index.unit_coordinate();
        point = self.storage.bounding_box().to_real_coordinate(&point);
        point
    }
    pub fn storage(&self) -> &SparseGridData<D>
    {
        &self.storage
    }

    pub fn points(&self) -> PointIterator<'_, D>
    {
        self.storage.points()        
    }

    pub fn values(&self) -> &Vec<[f64;DIM_OUT]>
    {
        &self.values
    }

    pub fn set_values(&mut self, values: Vec<[f64; DIM_OUT]>) -> Result<(), SGError>
    {
        if values.len() != self.len()
        {
            Err(SGError::NumberOfPointsAndValuesMismatch)
        }
        else
        {
            self.values = values;
            Ok(())
        }
    }
    pub fn coarsen<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, update_iterator_data: bool) -> usize
    {
        let mut total_num_removed = 0;
        loop
        {
            let r = algorithms::coarsening::coarsen(&mut self.storage, functor, &self.alpha, &self.values);
            if r.len() == self.alpha.len()
            {
                break;
            }
            let mut new_alpha = Vec::with_capacity(r.len());
            let mut new_values = Vec::with_capacity(r.len());
            for i in r
            {
                new_alpha.push(self.alpha[i]);
                new_values.push(self.values[i]);                
            }
            total_num_removed += self.alpha.len() - new_alpha.len(); 
            self.alpha = new_alpha;
            self.values = new_values;
                
        }
        if update_iterator_data
        {
            self.storage.generate_map();
            self.storage.generate_adjacency_data();
        }        
        total_num_removed
        
    }

    pub fn refine_iteration(&mut self, functor: &impl RefinementFunctor<D, DIM_OUT>) -> Vec<[f64; D]>
    {
        let ref_op = BaseRefinement(self.storage.has_boundary());
        let indices = ref_op.refine(&mut self.storage, &self.alpha, &self.values, functor);
        self.values.resize(self.len(), [0.0; DIM_OUT]);            
        let mut points = Vec::with_capacity(indices.len());
        for &i in &indices
        {
            let mut point = self.storage[i].unit_coordinate();
            point = self.storage.bounding_box().to_real_coordinate(&point);            
            points.push(point);
        }           
        points
    }
    pub fn refine<F: RefinementFunctor<D, DIM_OUT>, OP: HierarchisationOperation<D, DIM_OUT>>(&mut self, functor: &F, eval_fun: &mut dyn FnMut(&[f64;D])->[f64; DIM_OUT], op: &OP, max_iterations: usize) 
    {
        let ref_op = BaseRefinement(self.storage.has_boundary());
        let mut iteration = 0;     
        loop
        {
            let indices = ref_op.refine(&mut self.storage, &self.alpha, &self.values, functor);
            if indices.is_empty()
            {
                break;
            }
            self.values.reserve(indices.len());            
            for &i in &indices
            {
                let mut point = self.storage[i].unit_coordinate();
                point = self.storage.bounding_box().to_real_coordinate(&point);
                self.values.push(eval_fun(&point));
            }            
            self.hierarchize(op);
            iteration += 1;            
            if max_iterations <= iteration
            {
                break;
            }
        }
        self.sort();        
    }

    pub fn refine_parallel<F: RefinementFunctor<D, DIM_OUT>, OP: HierarchisationOperation<D, DIM_OUT>, EF: Fn(&[f64;D]) -> [f64; DIM_OUT] + Send + Sync>(&mut self, functor: &F, 
        eval_fun: &EF, op: &OP, max_iterations: usize) 
    {
        let ref_op = BaseRefinement(self.storage.has_boundary());
        let mut iteration = 0;     
        loop
        {
            let indices = ref_op.refine(&mut self.storage, &self.alpha, &self.values, functor);
            if indices.is_empty()
            {
                break;
            }
            let mut temp_values = vec![[0.0; DIM_OUT]; indices.len()];
            
            indices.par_iter().zip(temp_values.par_iter_mut()).for_each(
            |(&index, value)|
            {
                let mut point = self.storage[index].unit_coordinate();
                point = self.storage.bounding_box().to_real_coordinate(&point);
                *value = eval_fun(&point);
            }
            );        
            self.values.extend(&temp_values);
            self.hierarchize(op);
            iteration += 1;            
            if max_iterations <= iteration
            {
                break;
            }
        }
        self.sort();        
    }

  
    ///
    /// Saves full grid information (including maps and iterator data).
    /// Larger than compact form but faster than regenerating (compressed using LZ4).
    /// 
    pub fn save(&mut self, path: &str) -> Result<(), SGError>
    {
        let mut file = std::io::BufWriter::new(std::fs::File::create(path).map_err(|_|SGError::FileIOError)?);        
        let buffer = lz4_flex::compress_prepend_size(&bincode::serialize(&self).map_err(|_|SGError::SerializationFailed)?);
        file.write_all(&buffer).map_err(|_|SGError::WriteBufferFailed)?;
        Ok(())
    }

    
    ///
    /// Reads full grid information (including maps and iterator data).
    /// Much larger than compact form but faster than regenerating.
    /// 
    pub fn read_buffer(buffer: &[u8]) -> Result<Self, SGError>
    {      
        let buffer = lz4_flex::decompress_size_prepended(buffer).map_err(|_|SGError::LZ4DecompressionFailed)?;
        bincode::deserialize(&buffer).map_err(|_|SGError::DeserializationFailed)
    }

    ///
    /// Reads full grid information (including maps and iterator data).
    /// Much larger than compact form but faster than regenerating.
    /// 
    pub fn read<Reader: std::io::Read>(mut reader: Reader)  -> Result<Self, SGError>
    {
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).map_err(|_|SGError::ReadBufferFailed)?;
        let buffer = lz4_flex::decompress_size_prepended(&bytes).map_err(|_|SGError::LZ4DecompressionFailed)?;
        bincode::deserialize(&buffer).map_err(|_|SGError::DeserializationFailed)
    }
}