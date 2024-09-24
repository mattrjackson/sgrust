use crate::{basis::base::Basis, errors::SGError, iterators::grid_iterator::GridIterator, storage::linear_grid::SparseGridStorage};


pub struct BasisEvaluation<'a, const D: usize, const DIM_OUT: usize, BASIS: Basis>(pub &'a SparseGridStorage<D>, pub [BASIS; D]);

impl <'a, const D: usize, const DIM_OUT: usize, BASIS: Basis> BasisEvaluation<'a, D, DIM_OUT, BASIS>
{
    pub(crate) fn basis(&self, dim: usize) -> &BASIS
    {
        &self.1[dim]
    }
    fn recursive_eval(&self, point: [f64; D], current_dim: usize, value: f64, iterator: &mut GridIterator<D>, source: [u32; D], alpha: &[[f64; DIM_OUT]]) -> [f64; DIM_OUT]
    {
        const MAX_LEVEL: u32 = 31;
        let idx = source[current_dim];
        let mut level = 1;
        let mut result = [0.0; DIM_OUT];
        loop 
        {
            if let Some(seq) = iterator.seq()
            {
                let index = iterator.index().index[current_dim];
                let val = self.basis(current_dim).eval(level, index, point[current_dim]) * value;
                if current_dim == D - 1
                {
                    (0..DIM_OUT).for_each(|d| {
                        result[d] += alpha[seq][d] * val;
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
                if iterator.hint()
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
                    iterator.right_child(current_dim);
                } 
                else 
                {
                    iterator.left_child(current_dim);
                }
            }
            else
            {
                break;
            }
        }
        iterator.reset_to_level_one(current_dim);
        result
    }

    pub fn eval(&self, point: [f64; D], alpha: &[[f64; DIM_OUT]]) -> Result<[f64; DIM_OUT], SGError> 
    {
        let mut iterator = GridIterator::new(self.0);
        let bits = std::mem::size_of::<u32>() * 8;
        let bbox: Option<crate::storage::linear_grid::BoundingBox<D>> = self.0.bounding_box();
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
        Ok(self.recursive_eval(unit_coord, 0, 1.0, &mut iterator, source, alpha))
    }
}