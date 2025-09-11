use crate::utilities::float::Float;
use crate::{basis::{base::Basis, linear::LinearBasis}, dynamic::{iterators::dynamic_grid_iterator::GridIteratorT, storage::SparseGridData}, errors::SGError};


pub struct BasisEvaluation<'a>(pub &'a SparseGridData, pub Vec<LinearBasis>, pub usize, pub usize);

impl BasisEvaluation<'_>
{
    #[inline]
    pub(crate) fn basis(&self, dim: usize) -> &LinearBasis
    {
        &self.1[dim]
    }
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn recursive_eval<Iterator: GridIteratorT, T: Float +std::ops::AddAssign> (&self, point: &[f64], current_dim: usize, 
        value: T, iterator: &mut Iterator, source: &Vec<u32>, alpha: &[T], result: &mut [T]) -> Result<(), SGError>
    {
        #[allow(non_snake_case)]
        let D = self.2;
        #[allow(non_snake_case)]
        let DIM_OUT = self.3;
        const MAX_LEVEL: u32 = 31;
        let idx = source[current_dim];
        let mut level = 1;      
        
        loop
        {
            let index = iterator.point_index(current_dim);
            let val = T::from(self.basis(current_dim).eval(level, index, point[current_dim])) * value;            
            if current_dim == D - 1
            {
                if iterator.index().is_none()
                {
                    return Err(SGError::InvalidIteratorSequence);
                }
                let node_index = iterator.index().ok_or_else(||SGError::InvalidIteratorSequence)?;
                #[allow(clippy::needless_range_loop)]
                for d in 0..DIM_OUT
                {
                    result[d] += alpha[node_index*DIM_OUT+d] * val;
                }
            }
            else 
            {
                self.recursive_eval(point, current_dim + 1, val, iterator, source, alpha, result)?;                
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
        Ok(())
    }
    #[inline]
    pub fn eval<Iterator: GridIteratorT, T: Float +std::ops::AddAssign>(&self, point: &[f64], alpha: &[T], 
        iterator: &mut Iterator, result: &mut [T]) -> Result<(), SGError> 
    {        
        #[allow(non_snake_case)]
        let D = self.2;
        #[allow(non_snake_case)]
        let _DIM_OUT = self.3;
        let bits = std::mem::size_of::<u32>() * 8;
        if !self.0.bounding_box.contains(&point)
        {
            return Err(SGError::OutOfDomain);
        }
        let unit_coord = self.0.bounding_box.to_unit_coordinate(&point);
        let mut source = vec![0_u32; D];

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
        self.recursive_eval(&unit_coord, 0,  T::from(1.0), iterator, &source, alpha, result)
    }
}