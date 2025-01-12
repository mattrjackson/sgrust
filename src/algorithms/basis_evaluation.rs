use num_traits::Float;
use crate::{basis::base::Basis, errors::SGError, iterators::grid_iterator::GridIteratorT, storage::linear_grid::SparseGridData};


pub struct BasisEvaluation<'a, const D: usize, const DIM_OUT: usize, BASIS: Basis>(pub &'a SparseGridData<D>, pub [BASIS; D]);

impl <const D: usize, const DIM_OUT: usize, BASIS: Basis> BasisEvaluation<'_, D, DIM_OUT, BASIS>
{
    pub(crate) fn basis(&self, dim: usize) -> &BASIS
    {
        &self.1[dim]
    }
    fn recursive_eval<Iterator: GridIteratorT<D>, T: Float +std::ops::AddAssign> (&self, point: [f64; D], current_dim: usize, 
        value: T, iterator: &mut Iterator, source: [u32; D], alpha: &[[T; DIM_OUT]]) -> [T; DIM_OUT]
    {
        const MAX_LEVEL: u32 = 31;
        let idx = source[current_dim];
        let mut level = 1;
        let mut result = [T::zero(); DIM_OUT];
        
        loop
        {
            let index = iterator.node().index[current_dim];
            let val = T::from(self.basis(current_dim).eval(level, index, point[current_dim])).unwrap() * value;
            if current_dim == D - 1
            {
                (0..DIM_OUT).for_each(|d| {
                    result[d] += alpha[iterator.index()][d] * val;
                });
            }
            else 
            {
                let r_temp = self.recursive_eval(point, current_dim + 1, val, iterator, source, alpha);
                for i in 0..DIM_OUT
                {
                    result[i] += r_temp[i];
                }
            }
            if iterator.is_leaf()
            {
                break;
            }
            // this decides in which direction we should descend by evaluating
            // the corresponding bit
            // the bits are coded from left to right starting with level 1
            // being in position max_level
            let go_right = (idx & (1 << (MAX_LEVEL - level))) > 0;
            level += 1;

            if go_right 
            {
                if !iterator.right_child(current_dim)
                {
                    break;
                }
            } 
            else if !iterator.left_child(current_dim)
            {
                break;
            }            
        }        
        iterator.reset_to_level_one(current_dim);
        result
    }

    pub fn eval<Iterator: GridIteratorT<D>, T: Float +std::ops::AddAssign>(&self, point: [f64; D], alpha: &[[T; DIM_OUT]], 
        iterator: &mut Iterator) -> Result<[T; DIM_OUT], SGError> 
    {        
        let bits = std::mem::size_of::<u32>() * 8;
        let bbox: Option<crate::storage::linear_grid::BoundingBox<D>> = self.0.bounding_box;
        let unit_coord =
        if let Some(bbox) = bbox
        {
            if !bbox.contains(&point)
            {
                return Err(SGError::OutOfDomain);
            }
            bbox.to_unit_coordinate(&point)
        }
        else
        {
            point
        };
        let mut source = [0_u32; D];

        for d in 0..D
        {
            let val = (unit_coord[d]*(1 << (bits - 2)) as f64).floor() * 2.0;
            if unit_coord[d] == 1.0
            {
                source[d] = (val - 1.0) as u32;
            }
            else 
            {
                source[d] = (val + 1.0) as u32;
            }
        }
        Ok(self.recursive_eval(unit_coord, 0,  T::from(1.0).unwrap(), iterator, source, alpha))
    }
}