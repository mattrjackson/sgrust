use std::io::Write;

use bincode::config::standard;

//use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
use serde_with::serde_as;
use crate::dynamic::algorithms::basis_evaluation::BasisEvaluation;
use crate::dynamic::algorithms::refinement::{BaseRefinement, RefinementFunctor, RefinementMode, RefinementOptions, SparseGridRefinement};
use crate::basis::linear::LinearBasis;
use crate::errors::SGError;
use crate::dynamic::algorithms::hierarchisation::HierarchisationOperation;
use crate::dynamic::iterators::grid_iterator_cache::AdjacencyGridIterator;
use crate::dynamic::storage::{BoundingBox, GridPoint, PointIterator, SparseGridData};
use crate::dynamic::generators::*;
use serde::{Serialize,Deserialize};
#[serde_as]
#[derive(Serialize, Deserialize, Clone)]
pub struct SparseGridBase
{
    pub(crate) storage: SparseGridData,
    pub(crate) alpha: Vec<f64>,
    pub(crate) values: Vec<f64>,    
}


impl SparseGridBase
{
    pub fn new(num_inputs: usize, num_outputs: usize) -> Self
    {
        SparseGridBase { storage: SparseGridData::new(num_inputs, num_outputs), alpha: Vec::new(), values: Vec::new() }
    }

    pub fn alpha(&self) -> &[f64]
    {
        &self.alpha
    }

    pub fn alpha_mut(&mut self) -> &mut [f64]
    {
        &mut self.alpha
    }

    pub fn bounding_box(&self) -> &BoundingBox
    {        
        self.storage.bounding_box()
    }

    pub fn bounding_box_mut(&mut self) -> &mut BoundingBox
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

    pub fn sparse_grid<GENERATOR: Generator>(&mut self, levels: &[usize], generator: &GENERATOR)
    {        
        generator.regular(&mut self.storage, levels, None);       
        self.sort();
    }    

    pub fn full_grid<GENERATOR: Generator>(&mut self, level: usize, generator: &GENERATOR)
    {
        generator.full(&mut self.storage, level);
        self.sort();
    }

    pub fn sparse_grid_with_boundaries<GENERATOR: Generator>(&mut self, levels: &[usize], generator: &GENERATOR)
    {
        generator.regular_with_boundaries(&mut self.storage, levels, Some(1), None);
        self.storage.has_boundary = true;
        self.sort();
    }    

    pub fn full_grid_with_boundaries<GENERATOR: Generator>(&mut self, level: usize, generator: &GENERATOR)
    {
        generator.full_with_boundaries(&mut self.storage, level);
        self.storage.has_boundary = true;
        self.sort();
    }
    pub fn hierarchize<OP: HierarchisationOperation>(&mut self, op: &OP)
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
        let num_inputs = self.storage.num_inputs;
        let num_outputs = self.storage.num_outputs;
        let node_vec: Vec<_> = self.storage.nodes().collect();
        let mut indices = Self::argsort(&node_vec);        
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
                for n in 0..num_outputs
                {
                    if self.alpha.len() == indices.len()
                    {
                        self.alpha.swap(current*num_outputs+n, next*num_outputs+n);
                    }
                    if self.values.len() == indices.len()
                    {
                        self.values.swap(current*num_outputs+n, next*num_outputs+n);
                    } 
                }
                for n in 0..num_inputs
                {                   
                    self.storage.index.swap(current*num_inputs+n, next*num_inputs+n);
                    self.storage.level.swap(current*num_inputs+n, next*num_inputs+n);                                        
                }
                self.storage.flags.swap(current, next);
                
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
    pub fn interpolate(&self, x: &[f64], result: &mut [f64] ) -> Result<(), SGError>
    {    
        if !self.bounding_box().contains(&x)
        {
            Err(SGError::OutOfDomain)
        }
        else 
        {
            self.interpolate_unchecked(x, result)    
        }
    }
    #[inline]
    pub fn interpolate_unchecked(&self, x: &[f64], result: &mut [f64]) -> Result<(), SGError>
    {
        use crate::dynamic::algorithms::interpolation::InterpolationOperation;
        let iterator = &mut AdjacencyGridIterator::new(&self.storage);
        let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, vec![LinearBasis; self.storage.num_inputs], self.storage.num_inputs, self.storage.num_outputs));      
        op.interpolate(x, &self.alpha, iterator, result)
    }

    #[cfg(feature="rayon")]
    #[inline]
    pub fn interpolate_batch(&self, x: &[f64]) -> Result<Vec<f64>, SGError>
    {
       use rayon::{iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator}, slice::{ParallelSlice, ParallelSliceMut}};
       use crate::dynamic::algorithms::interpolation::InterpolationOperation;
       let mut results = vec![0.0; x.len() / self.storage.num_inputs *self.storage.num_outputs];       
       let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, vec![LinearBasis; self.storage.num_inputs], self.storage.num_inputs, self.storage.num_outputs));
        x.par_chunks_exact(self.storage.num_inputs).zip(results.par_chunks_exact_mut(self.storage.num_outputs)).try_for_each(
            |(x, y)|
            {
                let iterator = &mut AdjacencyGridIterator::new(&self.storage);
                
                if !self.bounding_box().contains(x)
                {
                    Err(SGError::OutOfDomain)
                }
                else
                {
                    op.interpolate(x, &self.alpha, iterator, y)?;
                    Ok(())
                }               
            }
        )?;       
       Ok(results)
    }

    #[cfg(feature="rayon")]
    #[inline]
    pub fn interpolate_batch_unchecked(&self, x: &[f64]) -> Result<Vec<f64>, SGError>
    {
        use rayon::{iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator}, slice::{ParallelSlice, ParallelSliceMut}};
        use crate::dynamic::algorithms::interpolation::InterpolationOperation;
        let mut results = vec![0.0; x.len() / self.storage.num_inputs *self.storage.num_outputs];       
        let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, vec![LinearBasis; self.storage.num_inputs], self.storage.num_inputs, self.storage.num_outputs));        
        x.par_chunks_exact(self.storage.num_inputs).zip(results.par_chunks_exact_mut(self.storage.num_outputs)).for_each(
            |(x, y)|
            {
                let iterator = &mut AdjacencyGridIterator::new(&self.storage);        
                op.interpolate(x, &self.alpha, iterator, y).unwrap();
            }
        );       
       Ok(results)
    }

    pub fn integrate_isotropic(&self) -> Vec<f64>
    {
        let mut result = vec![0.0; self.storage.num_outputs];
        LinearBasis::eval_point_dynamic(self.storage(), &self.alpha, &mut result);
        result
    }

    pub fn to_real_coordinate(&self, index: &GridPoint) -> Vec<f64>
    {
        let mut point = index.unit_coordinate();
        point = self.storage.bounding_box().to_real_coordinate(&point);
        point
    }
    pub fn storage(&self) -> &SparseGridData
    {
        &self.storage
    }

    pub fn points(&self) -> PointIterator<'_>
    {
        self.storage.points()        
    }

    pub fn values(&self) -> &Vec<f64>
    {
        &self.values
    }

    pub fn set_values(&mut self, values: Vec<f64>) -> Result<(), SGError>
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
    pub fn coarsen<F: RefinementFunctor>(&mut self, functor: &F, update_iterator_data: bool, threshold: f64) -> usize
    {
        let mut total_num_removed = 0;
        let mut last_num_removed: usize;
        loop
        {            
            last_num_removed =  self.coarsen_iteration(functor, threshold);                
            total_num_removed += last_num_removed;
            if last_num_removed == 0
            {
                break;
            }
        }
        if update_iterator_data
        {
            self.storage.generate_map();
            self.storage.generate_adjacency_data();
        }        
        total_num_removed
        
    }

    fn coarsen_iteration<F: RefinementFunctor>(&mut self, functor: &F, threshold: f64) -> usize
    {
        let mut total_num_removed = 0;
  
        let r = crate::dynamic::algorithms::coarsening::coarsen(&mut self.storage, functor, &self.alpha, &self.values, threshold);
        if r.len() != self.alpha.len()
        {
            let mut new_alpha = Vec::with_capacity(r.len());
            let mut new_values = Vec::with_capacity(r.len());
            for i in r
            {
                new_alpha.push(self.alpha[i]);
                new_values.push(self.values[i]);                
            }
            total_num_removed = self.alpha.len() - new_alpha.len(); 
            self.alpha = new_alpha;
            self.values = new_values;
        }
        total_num_removed
        
    }

    pub fn refine_iteration(&mut self, functor: &impl RefinementFunctor, options: RefinementOptions) -> Vec<f64>
    {
        if options.refinement_mode == RefinementMode::Anisotropic
        {                 
            self.coarsen_iteration(functor, options.threshold);
        }
        let ref_op = BaseRefinement(self.storage.has_boundary());
        let indices = ref_op.refine(&mut self.storage, &self.alpha, &self.values, functor, options.clone());
        self.values.resize(self.len()*self.storage.num_outputs, 0.0);            
        let mut points = Vec::with_capacity(indices.len()*self.storage.num_inputs);
        for &i in &indices
        {
            let mut point = self.storage.unit_coordinate(i);
            self.storage.bounding_box().to_real_coordinate_in_place(&mut point);            
            points.extend(point);
        }           
        points
    }
    pub fn refine<F: RefinementFunctor, OP: HierarchisationOperation, EF: Fn(&[f64])->Vec<f64>>(&mut self, functor: &F, eval_fun: &EF, op: &OP, options: RefinementOptions, max_iterations: usize) 
    {
        let ref_op = BaseRefinement(self.storage.has_boundary());
        let mut iteration = 0;     
        loop
        {
            if iteration > 0 && options.refinement_mode == RefinementMode::Anisotropic
            {                 
                self.coarsen_iteration(functor, options.threshold);                
            }
            let indices = ref_op.refine(&mut self.storage, &self.alpha, &self.values, functor, options.clone());
            if indices.is_empty()
            {
                break;
            }
            self.values.reserve(indices.len());            
            for &i in &indices
            {
                let mut point = self.storage.unit_coordinate(i);
                point = self.storage.bounding_box().to_real_coordinate(&point);
                self.values.extend(eval_fun(&point));
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

    #[cfg(feature="rayon")]
    pub fn refine_parallel<F: RefinementFunctor, OP: HierarchisationOperation, EF: Fn(&[f64]) -> Vec<f64> + Send + Sync>(&mut self, functor: &F, 
        eval_fun: &EF, op: &OP, options: RefinementOptions, max_iterations: usize) 
    {
        use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
        let ref_op = BaseRefinement(self.storage.has_boundary());
        let mut iteration = 0;     
        loop
        {
            use rayon::slice::ParallelSliceMut;

            if iteration > 0 && options.refinement_mode == RefinementMode::Anisotropic
            {                 
                self.coarsen_iteration(functor, options.threshold);                
            }
            let indices = ref_op.refine(&mut self.storage, &self.alpha, &self.values, functor, options.clone());
            if indices.is_empty()
            {
                break;
            }
            let mut temp_values = vec![0.0; self.storage.num_outputs * indices.len()];
            
            indices.par_iter().zip(temp_values.par_chunks_exact_mut(self.storage.num_outputs)).for_each(
            |(&index, value)|
            {
                let mut point = self.storage.unit_coordinate(index);
                self.storage.bounding_box().to_real_coordinate_in_place(&mut point);
                value.copy_from_slice(&eval_fun(&point));
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
        let buffer = lz4_flex::compress_prepend_size(&bincode::serde::encode_to_vec(&self, standard()).map_err(|_|SGError::SerializationFailed)?);
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
        let (grid, _size) = bincode::serde::decode_from_slice(&buffer, standard()).map_err(|_|SGError::DeserializationFailed)?;
        Ok(grid)
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
        let (grid, _size) = bincode::serde::decode_from_slice(&buffer, standard()).map_err(|_|SGError::DeserializationFailed)?;
        Ok(grid)
    }
}