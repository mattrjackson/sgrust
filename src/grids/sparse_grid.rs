use std::io::Write;

use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
use serde_with::serde_as;
use crate::algorithms::basis_evaluation::BasisEvaluation;
use crate::algorithms::integration::{AnisotropicQuadrature, IsotropicQuadrature};
use crate::algorithms::refinement::{BaseRefinement, RefinementFunctor, SparseGridRefinement};
use crate::basis::linear::LinearBasis;
use crate::errors::SGError;
use crate::hierarchisation::HierarchisationOperation;
use crate::iterators::grid_iterator_cache::{GridIteratorData, GridIteratorWithCache};
use crate::storage::linear_grid::{BoundingBox, GridPoint, SparseGridData};
use crate::{base_generator::Generator, storage::linear_grid::SparseGridStorage};
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
    fn bounding_box(&self) -> Option<BoundingBox<D>>
    {
        self.base().bounding_box()
    }

    /// Get the bounding box for this grid (mutable).
    fn bounding_box_mut(&mut self) -> &mut Option<BoundingBox<D>>
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

    /// Generate the runtime data that is required for interpolation.
    fn update_runtime_data(&mut self)
    {
        self.base_mut().update_runtime_data();
    }

    /// Retrieve the underlying storage data.
    fn storage(&self) -> &SparseGridStorage<D>
    {
        &self.base().storage
    }

    /// Get copy of points that make up this grid.
    fn points(&self) -> Vec<[f64; D]>
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
        self.update_runtime_data();
        Ok(())
    }

    /// Coarsen grid based on functor `F`
    fn coarsen<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F) -> usize
    {
        self.base_mut().coarsen(functor, true)
    }

    /// Refine grid based on function `F`
    fn refine<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, eval_fun: &mut dyn FnMut(&[f64;D])->[f64; DIM_OUT], max_iterations: usize);

    /// Save compact form of data (no runtime iterator data is stored)
    fn save_compact(&self, path: &str) -> Result<(), SGError>
    {
        self.base().save_compact(path)
    }
    
    /// Save full form of data (runtime iterator data is stored)
    fn save_full(&mut self, path: &str) -> Result<(), SGError>
    {
        self.base_mut().save_full(path)
    }
    /// Read compact form of data.
    fn read_compact<Reader: std::io::Read>(reader: Reader) -> Result<Self, SGError> where Self: Sized;

    /// Read full form of data.
    fn read_full<Reader: std::io::Read>(reader: Reader) -> Result<Self, SGError> where Self: Sized;    

    /// Read compact form of data from a buffer.
    fn read_compact_buffer(buffer: &[u8]) -> Result<Self, SGError> where Self: Sized;

    /// Read full form of data from a buffer.
    fn read_full_buffer(buffer: &[u8]) -> Result<Self, SGError> where Self: Sized;

}

#[serde_as]
#[derive(Serialize, Deserialize)]
pub struct SparseGridDataAndValues<const D: usize, const DIM_OUT: usize>
{
    data: SparseGridData<D>,
    #[serde_as(as = "Vec<[_; DIM_OUT]>")]
    alpha: Vec<[f64; DIM_OUT]>,
    #[serde_as(as = "Vec<[_; DIM_OUT]>")]
    values: Vec<[f64; DIM_OUT]>,
}
#[serde_as]
#[derive(Default, Serialize, Deserialize, Clone)]
pub struct SparseGridBase<const D: usize, const DIM_OUT: usize>
{
    pub(crate) storage: SparseGridStorage<D>,
    pub(crate) iterator_data: Option<GridIteratorData<D>>,
    #[serde_as(as = "Vec<[_; DIM_OUT]>")]
    pub(crate) alpha: Vec<[f64; DIM_OUT]>,
    #[serde_as(as = "Vec<[_; DIM_OUT]>")]
    pub(crate) values: Vec<[f64; DIM_OUT]>,    
}
impl<const D: usize, const DIM_OUT: usize> From<SparseGridDataAndValues<D, DIM_OUT>> for SparseGridBase<D, DIM_OUT>
{
    fn from(value: SparseGridDataAndValues<D, DIM_OUT>) -> Self 
    {
        Self { storage: SparseGridStorage::new(value.data), iterator_data: None, alpha: value.alpha, values: value.values }
    }
}
impl<const D: usize, const DIM_OUT: usize> SparseGridBase<D, DIM_OUT>
{
    pub fn new() -> Self
    {
        SparseGridDataAndValues { data: SparseGridData::new(), alpha: Vec::new(), values: Vec::new() }.into()
    }

    pub fn alpha(&self) -> &[[f64; DIM_OUT]]
    {
        &self.alpha
    }

    pub fn alpha_mut(&mut self) -> &mut Vec<[f64; DIM_OUT]>
    {
        &mut self.alpha
    }

    pub fn bounding_box(&self) -> Option<BoundingBox<D>>
    {        
        self.storage.data.bounding_box
    }

    pub fn bounding_box_mut(&mut self) -> &mut Option<BoundingBox<D>>
    {        
        &mut self.storage.data.bounding_box
    }

    pub fn is_empty(&self) -> bool
    {
        self.storage.data.is_empty()
    }

    pub fn len(&self) -> usize
    {
        self.storage.data.len()
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
        self.storage.data.has_boundary = true;
        self.sort();
    }    

    pub fn full_grid_with_boundaries<GENERATOR: Generator<D>>(&mut self, level: usize, generator: &GENERATOR)
    {
        generator.full_with_boundaries(&mut self.storage, level);
        self.storage.data.has_boundary = true;
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
    fn sort(&mut self)
    {
        let mut indices = Self::argsort(self.storage.list());        
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
                self.storage.list_mut().swap(current, next);
                
                indices[current] = current;                
                current = next;                
            }
            indices[current] = current;
        }
        for (i, item) in self.storage.data.list.iter().enumerate()
        {            
            *self.storage.map.get_mut(item).unwrap() = i;
        }
    }

    pub fn update_runtime_data(&mut self)
    {
        self.iterator_data = Some(GridIteratorData::new(&self.storage));
    }
    #[inline]
    pub fn interpolate(&self, x: [f64; D]) -> Result<[f64; DIM_OUT], SGError>
    {
        if let Some(bbox) = self.bounding_box()
        {
            if !bbox.contains(&x)
            {
                return Err(SGError::OutOfDomain);
            }
        }
        else if !BoundingBox::default().contains(&x)
        {
            return Err(SGError::OutOfDomain);
        }
        self.interpolate_unchecked(x)
       
    }
    #[inline]
    pub fn interpolate_unchecked(&self, x: [f64; D]) -> Result<[f64; DIM_OUT], SGError>
    {
        use crate::algorithms::interpolation::InterpolationOperation;
        let data = self.iterator_data.as_ref().expect("Grid Iterator data must be computed before calling interpolate!");
        let iterator = &mut GridIteratorWithCache::new(data);
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
                let data = self.iterator_data.as_ref().expect("Grid Iterator data must be computed before calling interpolate!");
                let iterator = &mut GridIteratorWithCache::new(data);
                let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, [LinearBasis; D]));
                if let Some(bbox) = self.bounding_box()
                {
                    if !bbox.contains(x)
                    {
                        *y = Err(SGError::OutOfDomain);
                    }
                    else
                    {
                        *y = op.interpolate(*x, &self.alpha, iterator);
                    }
                }
                else if !BoundingBox::default().contains(x)
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
                let data = self.iterator_data.as_ref().expect("Grid Iterator data must be computed before calling interpolate!");
                let iterator = &mut GridIteratorWithCache::new(data);
                let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, [LinearBasis; D]));                
                *y = op.interpolate(*x, &self.alpha, iterator);
            }
        );       
       results
    }

    pub fn integrate_isotropic<QUADRATURE: IsotropicQuadrature<D, DIM_OUT>>(&self, op: &QUADRATURE) -> [f64; DIM_OUT]
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
        if let Some(bbox) = self.storage.bounding_box()
        {
            point = bbox.to_real_coordinate(&point);
        }
        point
    }
    pub fn storage(&self) -> &SparseGridStorage<D>
    {
        &self.storage
    }

    pub fn points(&self) -> Vec<[f64; D]>
    {
        self.storage.points()        
    }

    pub fn values(&self) -> &Vec<[f64;DIM_OUT]>
    {
        &self.values
    }

    pub fn set_values(&mut self, values: Vec<[f64; DIM_OUT]>) -> Result<(), SGError>
    {
        if values.len() != self.points().len()
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
            if r.is_empty()
            {
                break;
            }
            let mut new_alpha = Vec::with_capacity(self.alpha.len() - r.len());
            let mut new_values = Vec::with_capacity(self.alpha.len() - r.len());
            for i in 0..self.alpha.len()
            {
                if !r.contains(&i)
                {
                    new_alpha.push(self.alpha[i]);
                    new_values.push(self.values[i]);
                }
            }
            self.alpha = new_alpha;
            self.values = new_values;
            total_num_removed += r.len();     
        }
       // self.hierarchize();
        if update_iterator_data
        {
            self.update_runtime_data();
        }        
        total_num_removed
        
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
                let mut point = self.storage.list()[i].unit_coordinate();
                if let Some(bbox) = self.storage.bounding_box()
                {
                    point = bbox.to_real_coordinate(&point);
                }
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

    ///
    /// Save compact form of data. This requires
    /// regenerating iterator data, which can take seconds when grid
    /// is initially read.
    /// 
    pub fn save_compact(&self, path: &str) -> Result<(), SGError>
    {
        let file = std::io::BufWriter::new(std::fs::File::create(path).map_err(|_|SGError::FileIOError)?);        
        let storage = SparseGridDataAndValues{alpha: self.alpha.clone(), data: self.storage.data.clone(), values: self.values.clone()};
        bincode::serialize_into(file,   &storage).map_err(|_|SGError::SerializationFailed)?;
        Ok(())        
    }


    ///
    /// Saves full grid information (including maps and iterator data).
    /// Larger than compact form but faster than regenerating (compressed using LZ4).
    /// 
    pub fn save_full(&mut self, path: &str) -> Result<(), SGError>
    {
        // Save iterator data if it doesn't exist
        if self.iterator_data.is_none()
        {
            self.update_runtime_data();
        }
        let mut file = std::io::BufWriter::new(std::fs::File::create(path).map_err(|_|SGError::FileIOError)?);        
        let buffer = lz4_flex::compress_prepend_size(&bincode::serialize(&self).map_err(|_|SGError::SerializationFailed)?);
        file.write_all(&buffer).map_err(|_|SGError::WriteBufferFailed)?;
        Ok(())
    }

    ///
    /// Reads compact form of data. This requires
    /// regenerating iterator data, which can take seconds when grid
    /// is initially read.
    ///   
    pub fn read_compact_buffer(buffer: &[u8]) -> Result<Self, SGError>
    {      
        let data:  SparseGridDataAndValues<D, DIM_OUT> = bincode::deserialize_from(buffer).map_err(|_|SGError::DeserializationFailed)?;
        let mut grid: SparseGridBase<D, DIM_OUT> = data.into();
        grid.update_runtime_data();
        Ok(grid)
    }
    
    ///
    /// Reads full grid information (including maps and iterator data).
    /// Much larger than compact form but faster than regenerating.
    /// 
    pub fn read_full_buffer(buffer: &[u8]) -> Result<Self, SGError>
    {      
        let buffer = lz4_flex::decompress_size_prepended(buffer).map_err(|_|SGError::LZ4DecompressionFailed)?;
        bincode::deserialize(&buffer).map_err(|_|SGError::DeserializationFailed)
    }

    ///
    /// Reads full grid information (including maps and iterator data).
    /// Much larger than compact form but faster than regenerating.
    /// 
    pub fn read_full<Reader: std::io::Read>(mut reader: Reader)  -> Result<Self, SGError>
    {
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).map_err(|_|SGError::ReadBufferFailed)?;
        let buffer = lz4_flex::decompress_size_prepended(&bytes).map_err(|_|SGError::LZ4DecompressionFailed)?;
        bincode::deserialize(&buffer).map_err(|_|SGError::DeserializationFailed)
    }

    ///
    /// Reads compact form of data. This requires
    /// regenerating iterator data, which can take seconds when grid
    /// is initially read.
    /// 
    pub fn read_compact<Reader: std::io::Read>(reader: Reader)  -> Result<Self, SGError>
    {
        let data:  SparseGridDataAndValues<D, DIM_OUT> = bincode::deserialize_from(reader).map_err(|_|SGError::DeserializationFailed)?;
        let mut grid: SparseGridBase<D, DIM_OUT> = data.into();
        grid.update_runtime_data();
        Ok(grid)
    }
}