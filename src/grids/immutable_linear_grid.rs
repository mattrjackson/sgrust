use std::io::Write;

use num_traits::Float;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
use serde::Serialize;
use serde_with::serde_as;

use crate::{algorithms::{basis_evaluation::BasisEvaluation, integration::IsotropicQuadrature}, basis::linear::LinearBasis, errors::SGError, iterators::grid_iterator_cache::{GridIteratorData, GridIteratorWithCache}, storage::linear_grid::{BoundingBox, SparseGridData}};

use super::{linear_grid::LinearGrid, sparse_grid::SparseGrid};


///
/// Immutable variant of `LinearGrid` that only supports interpolation and integration. Smaller memory footprint, as we only store the data needed for
/// these two operations, as compared to `LinearGrid`, which supports refinement and coarsening operations.
/// 
#[serde_as]
#[derive(Default,Serialize, serde::Deserialize, Clone)]
pub struct ImmutableLinearGrid<T: Float + std::ops::AddAssign + serde::Serialize + for<'a> serde::Deserialize<'a> + Send + Sync, const D: usize, const DIM_OUT: usize>
{
    adjacency_data: GridIteratorData<D>,
    storage: SparseGridData<D>,
    #[serde_as(as = "Vec<[_; DIM_OUT]>")]
    alpha: Vec<[T; DIM_OUT]>,
    #[serde_as(as = "Vec<[_; DIM_OUT]>")]
    values: Vec<[T; DIM_OUT]>,     
}

impl<T: Float  + std::ops::AddAssign + serde::Serialize + for<'a> serde::Deserialize<'a> + Send + Sync, const D: usize, const DIM_OUT: usize> ImmutableLinearGrid<T, D, DIM_OUT>
{
    pub fn len(&self) ->usize
    {
        self.storage.len()
    }

    pub fn is_empty(&self) -> bool
    {
        self.storage.is_empty()
    }

    pub fn points(&self) -> Vec<[f64; D]>
    {
        let mut list = Vec::new();
        for index in self.storage.iter()
        {
            let point = index.unit_coordinate();
            list.push(self.bounding_box().to_real_coordinate(&point));
        }
        list
    }

    pub fn values(&self) -> &Vec<[T; DIM_OUT]>
    {
        &self.values
    }

    #[allow(unused)]
    pub(crate) fn new(alpha: Vec<[T; DIM_OUT]>, values: Vec<[T; DIM_OUT]>, adjacency_data: GridIteratorData<D>, storage: SparseGridData<D>) -> Self
    {
        Self { alpha, values, adjacency_data, storage }
    }
    #[allow(unused)]
    pub(crate) fn adjacency_data(&self) -> &GridIteratorData<D>
    {
        &self.adjacency_data
    }

    pub fn get_storage(&self) -> &SparseGridData<D>
    {
        &self.storage
    }

    pub fn bounding_box(&self) -> BoundingBox<D>
    {
        self.storage.bounding_box.unwrap_or_default()
    }

    pub fn get_num_points(&self) -> usize
    {
        self.storage.len()
    }

    #[inline]
    pub fn interpolate(&self, x: [f64; D]) -> Result<[T; DIM_OUT], SGError>
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
    pub fn interpolate_unchecked(&self, x: [f64; D]) -> Result<[T; DIM_OUT], SGError>
    {
        use crate::algorithms::interpolation::InterpolationOperation;        
        let iterator = &mut GridIteratorWithCache::new(&self.adjacency_data, &self.storage);
        let op = InterpolationOperation(self.storage.has_boundary, BasisEvaluation(&self.storage, [LinearBasis; D]));      
        op.interpolate(x, &self.alpha, iterator)       
    }

    #[inline]
    pub fn interpolate_batch(&self, x: &[[f64; D]]) -> Vec<Result<[T; DIM_OUT], SGError>>
    {
       use crate::algorithms::interpolation::InterpolationOperation;
       let mut results = vec![Ok([T::zero(); DIM_OUT]); x.len()];
        x.par_iter().zip(results.par_iter_mut()).for_each(
            |(x, y)|
            {
                let iterator = &mut GridIteratorWithCache::new(&self.adjacency_data, &self.storage);
                let op = InterpolationOperation(self.storage.has_boundary, BasisEvaluation(&self.storage, [LinearBasis; D]));
                
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
    pub fn interpolate_batch_unchecked(&self, x: &[[f64; D]]) -> Vec<Result<[T; DIM_OUT], SGError>>
    {
        use crate::algorithms::interpolation::InterpolationOperation;
       
       let mut results = vec![Ok([T::zero(); DIM_OUT]); x.len()];
        x.par_iter().zip(results.par_iter_mut()).for_each(
            |(x, y)|
            {
                let iterator: &mut GridIteratorWithCache<'_, '_, D> = &mut GridIteratorWithCache::new(&self.adjacency_data, &self.storage);
                let op = InterpolationOperation(self.storage.has_boundary, BasisEvaluation(&self.storage, [LinearBasis; D]));                
                *y = op.interpolate(*x, &self.alpha, iterator);
            }
        );       
       results
    }

    pub fn integrate(&self) -> [T; DIM_OUT]
    {
        LinearBasis.eval(&self.storage, &self.alpha)
    }

    ///
    /// Saves grid to file...
    /// 
    pub fn save(&self, path: &str) -> Result<(), SGError>
    {
        let mut file = std::io::BufWriter::new(std::fs::File::create(path).map_err(|_|SGError::FileIOError)?);        
        let buffer = lz4_flex::compress_prepend_size(&bincode::serialize(&self).map_err(|_|SGError::SerializationFailed)?);
        file.write_all(&buffer).map_err(|_|SGError::WriteBufferFailed)?;
        Ok(())
    }
    ///
    /// Reads full grid information from buffer
    /// 
    pub fn read_buffer(buffer: &[u8]) -> Result<Self, SGError>
    {      
        let buffer = lz4_flex::decompress_size_prepended(buffer).map_err(|_|SGError::LZ4DecompressionFailed)?;
        bincode::deserialize(&buffer).map_err(|_|SGError::DeserializationFailed)
    }

    ///
    /// Reads grid information
    /// 
    pub fn read<Reader: std::io::Read>(mut reader: Reader)  -> Result<Self, SGError>
    {
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).map_err(|_|SGError::ReadBufferFailed)?;
        let buffer = lz4_flex::decompress_size_prepended(&bytes).map_err(|_|SGError::LZ4DecompressionFailed)?;
        bincode::deserialize(&buffer).map_err(|_|SGError::DeserializationFailed)
    }
}

impl<T: Float  + std::ops::AddAssign + serde::Serialize + for<'a> serde::Deserialize<'a> + Send + Sync, 
    const D: usize, const DIM_OUT: usize> From<LinearGrid<D,DIM_OUT>> for ImmutableLinearGrid<T, D, DIM_OUT>
{
    fn from(value: LinearGrid<D,DIM_OUT>) -> Self {
        let adjacency_data = if let Some(data) = &value.0.iterator_data
        {
            data.clone()
        }
        else
        {
            GridIteratorData::new(value.storage())
        };
        let alpha: Vec<[T; DIM_OUT]> = value.alpha().iter().map(|alpha|
        {
            std::array::from_fn(|i|T::from(alpha[i]).unwrap())
        }).collect();

        let values: Vec<[T; DIM_OUT]> = value.values().iter().map(|value|
            {
                std::array::from_fn(|i|T::from(value[i]).unwrap())
            }).collect();
        Self { adjacency_data, storage: value.storage().data.clone(), alpha, values }
    }
}

#[test]
fn check_size_difference()
{
     // Build the 2D grid object, only one value per node.
     let mut grid = LinearGrid::<6,1>::default();
     // using a full grid here, but a sparse grid could be used instead...
     grid.full_grid_with_boundaries(3);
 
     let f = |x: &[f64; 6]|
     {
         let mut r = [0.0];
         (0..6).for_each(|i| {
             r[0] += x[i]*x[i]*x[i];
         });
         r
     };
 
     let values = grid.points().iter().map(f).collect();
     grid.set_values(values).unwrap();    

     grid.save_full("grid.bin").unwrap();
     let igrid: ImmutableLinearGrid<f32, 6,1> = grid.into();
     igrid.save("igrid.bin").unwrap();

     let grid_size = std::fs::metadata("grid.bin").unwrap().len();
     let igrid_size = std::fs::metadata("igrid.bin").unwrap().len();
     println!("size of LinearGrid={grid_size}");
     println!("size of ImmutableLinearGrid={igrid_size}");
     std::fs::remove_file("grid.bin").unwrap();
     std::fs::remove_file("igrid.bin").unwrap();
     assert!(igrid_size < grid_size);
}