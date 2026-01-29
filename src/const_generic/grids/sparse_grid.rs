use std::io::Write;



use serde_with::serde_as;
use crate::const_generic::algorithms::basis_evaluation::BasisEvaluation;
use crate::const_generic::algorithms::refinement::{BaseRefinement, RefinementFunctor, RefinementMode, RefinementOptions};
use crate::basis::linear::LinearBasis;
use crate::errors::SGError;
use crate::const_generic::algorithms::hierarchisation::HierarchisationOperation;
use crate::const_generic::iterators::grid_iterator_cache::AdjacencyGridIterator;
use crate::const_generic::storage::{BoundingBox, GridPoint, PointIterator, SparseGridData};
use crate::const_generic::generators::*;
use crate::const_generic::algorithms;
use serde::{Serialize,Deserialize};
#[serde_as]
#[derive(Default, Serialize, Deserialize, Clone)]
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
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

    pub fn sparse_grid<GENERATOR: Generator<D>>(&mut self, levels: [usize; D], generator: &GENERATOR) -> Result<(), SGError>
    {        
        generator.regular(&mut self.storage, levels, None)?;     
        self.sort();
        Ok(())
    }    

    pub fn full_grid<GENERATOR: Generator<D>>(&mut self, level: usize, generator: &GENERATOR) -> Result<(), SGError>
    {
        generator.full(&mut self.storage, level)?;
        self.sort();
        Ok(())
    }

    pub fn sparse_grid_with_boundaries<GENERATOR: Generator<D>>(&mut self, levels: [usize; D], generator: &GENERATOR) -> Result<(), SGError>
    {
        generator.regular_with_boundaries(&mut self.storage, levels, Some(1), None)?;
        self.storage.has_boundary = true;
        self.sort();
        Ok(())
    }    

    pub fn full_grid_with_boundaries<GENERATOR: Generator<D>>(&mut self, level: usize, generator: &GENERATOR) -> Result<(), SGError>
    {
        generator.full_with_boundaries(&mut self.storage, level)?;
        self.storage.has_boundary = true;
        self.sort();
        Ok(())
    }
    pub fn hierarchize<OP: HierarchisationOperation<D, DIM_OUT>>(&mut self, op: &OP) -> Result<(), SGError>
    {
        self.alpha.clone_from(&self.values);
        op.hierarchize(&mut self.alpha, &self.storage)       
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
            *value = indices_rev[*value as usize] as u32;
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
        use crate::const_generic::algorithms::interpolation::InterpolationOperation;
        if self.values.len() == 1
        {
            return Ok(self.values[0]);
        }
        let iterator = &mut AdjacencyGridIterator::new(&self.storage);
        let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, [LinearBasis; D]));      
        op.interpolate(x, &self.alpha, iterator)       
    }

    #[cfg(feature="rayon")]
    #[inline]
    pub fn interpolate_batch(&self, x: &[[f64; D]]) -> Vec<Result<[f64; DIM_OUT], SGError>>
    {
        use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
       use crate::const_generic::algorithms::interpolation::InterpolationOperation;
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

    #[cfg(feature="rayon")]
    #[inline]
    pub fn interpolate_batch_unchecked(&self, x: &[[f64; D]]) -> Vec<Result<[f64; DIM_OUT], SGError>>
    {
        use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
        use crate::const_generic::algorithms::interpolation::InterpolationOperation;
       
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

    pub fn integrate_isotropic(&self) -> [f64; DIM_OUT]
    {
        LinearBasis::eval_point(self.storage(), &self.alpha)
    }

    // pub fn eval_anisotropic_quadrature<QUADRATURE: AnisotropicQuadrature<D, DIM_OUT> + ?Sized>(&self, op: &QUADRATURE, index: usize, dim: usize) -> [f64; DIM_OUT]
    // {
    //     op.eval(self.storage(), index, dim)
    // }

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
    pub fn coarsen<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, update_iterator_data: bool, threshold: f64) -> usize
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

    fn coarsen_iteration<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, threshold: f64) -> usize
    {
        let mut total_num_removed = 0;
  
        let r = algorithms::coarsening::coarsen(&mut self.storage, functor, &self.alpha, &self.values, threshold);
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

    pub fn refine_iteration(&mut self, functor: &impl RefinementFunctor<D, DIM_OUT>, options: RefinementOptions) -> Vec<[f64; D]>
    {
        if options.refinement_mode == RefinementMode::Anisotropic
        {                 
            self.coarsen_iteration(functor, options.threshold);
        }
        let ref_op = BaseRefinement(self.storage.has_boundary());
        let indices = ref_op.refine(&mut self.storage, &self.alpha, &self.values, functor, options.clone());
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
    pub fn refine<F: RefinementFunctor<D, DIM_OUT>, OP: HierarchisationOperation<D, DIM_OUT>, EF: Fn(&[f64;D])->[f64; DIM_OUT]>(&mut self, functor: &F, eval_fun: &EF, op: &OP, options: RefinementOptions, max_iterations: usize) -> Result<(), SGError>
    {
        let ref_op = BaseRefinement(self.storage.has_boundary());
        let mut iteration = 1;  
        loop
        {
            let indices = ref_op.refine(&mut self.storage, &self.alpha, &self.values, functor, options.clone());
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
            self.hierarchize(op)?;
            iteration += 1;            
            if max_iterations <= iteration
            {
                break;
            }
        }
        self.sort();        
        Ok(())
    }

    #[cfg(feature="rayon")]
    pub fn refine_parallel<F: RefinementFunctor<D, DIM_OUT>, OP: HierarchisationOperation<D, DIM_OUT>, EF: Fn(&[f64;D]) -> [f64; DIM_OUT] + Send + Sync>(&mut self, functor: &F, 
        eval_fun: &EF, op: &OP, options: RefinementOptions, max_iterations: usize) 
    {
        use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
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
            self.hierarchize(op).expect("Could not hierarchize grid.");
            iteration += 1;            
            if max_iterations <= iteration
            {
                break;
            }
        }
        self.sort();        
    }

  
    ///
    /// Writes full grid information (including maps and iterator data)    
    /// 
    pub fn write(&mut self, path: &str, format: crate::serialization::SerializationFormat) -> Result<(), SGError>
    {
        let mut file = std::io::BufWriter::new(std::fs::File::create(path).map_err(|_|SGError::FileIOError)?);        
        let buffer = crate::serialization::serialize(self, format)?;
        file.write_all(&buffer).map_err(|_|SGError::WriteBufferFailed)?;
        Ok(())
    }

    ///
    /// Write grid to buffer with the specified serialization format.
    /// 
    pub fn write_buffer(&self, format: crate::serialization::SerializationFormat) -> Result<Vec<u8>, SGError>
    {     
        crate::serialization::serialize(self, format)
    }
    
    ///
    /// Reads full grid information (including maps and iterator data).
    /// 
    pub fn read_buffer(buffer: &[u8], format: crate::serialization::SerializationFormat) -> Result<Self, SGError>
    {      
        crate::serialization::deserialize(buffer, format)
    }

    ///
    /// Reads full grid information from a reader.
    /// 
    pub fn read<Reader: std::io::Read>(mut reader: Reader, format: crate::serialization::SerializationFormat) -> Result<Self, SGError>
    {
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).map_err(|_|SGError::ReadBufferFailed)?;
        Self::read_buffer(&bytes, format)
    }
}