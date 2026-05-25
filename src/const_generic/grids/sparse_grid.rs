use std::io::Write;



use serde_with::serde_as;
use crate::const_generic::algorithms::basis_evaluation::BasisEvaluation;
use crate::const_generic::algorithms::refinement::{BaseRefinement, RefinementFunctor, RefinementMode, RefinementOptions};
use crate::basis::linear::LinearBasis;
use crate::basis::linear::POW2_F64;
use crate::errors::SGError;
use crate::const_generic::algorithms::hierarchisation::HierarchisationOperation;
use crate::const_generic::iterators::grid_iterator_cache::AdjacencyGridIterator;
use crate::const_generic::storage::{BoundingBox, GridPoint, PointIterator, SparseGridData};
use crate::const_generic::generators::*;
use crate::const_generic::algorithms;
use serde::{Serialize,Deserialize};
#[cfg(feature = "rkyv")]
use rkyv::with::Skip;
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
    #[serde(skip_serializing, skip_deserializing)]
    #[cfg_attr(feature = "rkyv", rkyv(with = Skip))]
    interpolation_nodes_1d: Vec<f64>,
    #[serde(skip_serializing, skip_deserializing)]
    #[cfg_attr(feature = "rkyv", rkyv(with = Skip))]
    interpolation_order_1d: Vec<usize>,
    #[serde(skip_serializing, skip_deserializing)]
    #[cfg_attr(feature = "rkyv", rkyv(with = Skip))]
    interpolation_uniform_1d: bool,
    #[serde(skip_serializing, skip_deserializing)]
    #[cfg_attr(feature = "rkyv", rkyv(with = Skip))]
    interpolation_step_1d: f64,
}


impl<const D: usize, const DIM_OUT: usize> SparseGridBase<D, DIM_OUT>
{
    pub fn new() -> Self
    {
        SparseGridBase {
            storage: SparseGridData::default(),
            alpha: Vec::new(),
            values: Vec::new(),
            interpolation_nodes_1d: Vec::new(),
            interpolation_order_1d: Vec::new(),
            interpolation_uniform_1d: false,
            interpolation_step_1d: 0.0,
        }
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
        self.update_1d_interpolation_data();
    }

    pub(crate) fn update_1d_interpolation_data(&mut self)
    {
        if D > 1
        {
            return;
        }
        self.interpolation_nodes_1d.clear();
        self.interpolation_order_1d.clear();
        self.interpolation_uniform_1d = false;
        self.interpolation_step_1d = 0.0;

        let mut nodes = Vec::with_capacity(self.storage.len());
        for (seq, point) in self.storage.nodes().iter().enumerate()
        {
            let coord = point.index[0] as f64 / POW2_F64[point.level[0] as usize];
            nodes.push((coord, seq));
        }
        nodes.sort_by(|lhs, rhs| lhs.0.total_cmp(&rhs.0));

        self.interpolation_nodes_1d.reserve(nodes.len());
        self.interpolation_order_1d.reserve(nodes.len());
        for (coord, seq) in nodes
        {
            self.interpolation_nodes_1d.push(coord);
            self.interpolation_order_1d.push(seq);
        }
        self.update_uniform_1d_metadata();
    }

    fn update_uniform_1d_metadata(&mut self)
    {
        let len = self.interpolation_nodes_1d.len();
        if len == 0
        {
            return;
        }

        let denom = if self.storage.has_boundary()
        {
            if len < 2
            {
                return;
            }
            len - 1
        }
        else
        {
            len + 1
        };
        let step = 1.0 / denom as f64;
        let offset = usize::from(!self.storage.has_boundary());
        let is_uniform = self.interpolation_nodes_1d.iter().enumerate().all(|(i, &node)|
        {
            (node - (i + offset) as f64 * step).abs() < 1e-14
        });

        if is_uniform
        {
            self.interpolation_uniform_1d = true;
            self.interpolation_step_1d = step;
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
        if D == 1
        {
            return self.interpolate_1d_unchecked(x[0]);
        }
        let iterator = &mut AdjacencyGridIterator::new(&self.storage);
        let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, [LinearBasis; D]));      
        op.interpolate(x, &self.alpha, iterator)       
    }

    #[inline]
    fn linear_basis_value(level: u32, index: u32, x: f64) -> f64
    {
        if level == 0
        {
            (1.0 - x) * (1 - index) as f64 + x * index as f64
        }
        else
        {
            let dist = (x * POW2_F64[level as usize] - index as f64).abs();
            (1.0 - dist).max(0.0)
        }
    }

    #[inline]
    fn add_1d_contribution(&self, seq: usize, weight: f64, result: &mut [f64; DIM_OUT])
    {
        for output in 0..DIM_OUT
        {
            result[output] += self.alpha[seq][output] * weight;
        }
    }

    #[inline]
    fn descend_1d(&self, seq: usize, go_right: bool) -> Option<usize>
    {
        let adj = &self.storage.adjacency_data[seq];
        if go_right && adj.has_right_child()
        {
            Some((seq as i64 + adj.down_right()) as usize)
        }
        else if !go_right && adj.has_left_child()
        {
            Some((seq as i64 + adj.down_left()) as usize)
        }
        else
        {
            None
        }
    }

    fn interpolate_1d_unchecked(&self, x: f64) -> Result<[f64; DIM_OUT], SGError>
    {
        if self.storage.is_empty()
        {
            return Ok([0.0; DIM_OUT]);
        }
        let bbox = self.storage.bounding_box();
        let x = (x - bbox.lower[0]) / (bbox.upper[0] - bbox.lower[0]);
        if self.interpolation_order_1d.len() == self.storage.len()
            && self.interpolation_nodes_1d.len() == self.storage.len()
            && self.values.len() == self.storage.len()
        {
            return Ok(self.interpolate_1d_sorted_unit(x));
        }
        if self.storage.has_boundary()
        {
            Ok(self.interpolate_1d_boundary_unit(x))
        }
        else
        {
            Ok(self.interpolate_1d_inner_unit(x))
        }
    }

    fn interpolate_1d_sorted_unit(&self, x: f64) -> [f64; DIM_OUT]
    {
        if self.interpolation_uniform_1d
        {
            return self.interpolate_1d_uniform_unit(x);
        }

        match self.interpolation_nodes_1d.binary_search_by(|node| node.total_cmp(&x))
        {
            Ok(pos) => self.values[self.interpolation_order_1d[pos]],
            Err(pos) => {
                let mut result = [0.0; DIM_OUT];
                let left = if pos == 0 {
                    None
                } else {
                    Some((self.interpolation_nodes_1d[pos - 1], self.interpolation_order_1d[pos - 1]))
                };
                let right = if pos == self.interpolation_nodes_1d.len() {
                    None
                } else {
                    Some((self.interpolation_nodes_1d[pos], self.interpolation_order_1d[pos]))
                };

                let (left_x, right_x) = match (left, right)
                {
                    (Some((left_x, _)), Some((right_x, _))) => (left_x, right_x),
                    (None, Some((right_x, _))) if !self.storage.has_boundary() => (0.0, right_x),
                    (Some((left_x, _)), None) if !self.storage.has_boundary() => (left_x, 1.0),
                    (None, Some((_, right_seq))) => return self.values[right_seq],
                    (Some((_, left_seq)), None) => return self.values[left_seq],
                    (None, None) => return result,
                };

                let width = right_x - left_x;
                if width <= 0.0
                {
                    return result;
                }
                let right_weight = (x - left_x) / width;
                let left_weight = 1.0 - right_weight;

                if let Some((_, left_seq)) = left
                {
                    for output in 0..DIM_OUT
                    {
                        result[output] += self.values[left_seq][output] * left_weight;
                    }
                }
                if let Some((_, right_seq)) = right
                {
                    for output in 0..DIM_OUT
                    {
                        result[output] += self.values[right_seq][output] * right_weight;
                    }
                }
                result
            }
        }
    }

    fn interpolate_1d_uniform_unit(&self, x: f64) -> [f64; DIM_OUT]
    {
        let mut result = [0.0; DIM_OUT];
        let len = self.interpolation_order_1d.len();

        if self.storage.has_boundary()
        {
            if x <= 0.0
            {
                return self.values[self.interpolation_order_1d[0]];
            }
            if x >= 1.0
            {
                return self.values[self.interpolation_order_1d[len - 1]];
            }
        }
        else if x <= 0.0 || x >= 1.0
        {
            return result;
        }

        let scaled = x / self.interpolation_step_1d;
        let left_lattice = if self.storage.has_boundary()
        {
            (scaled.floor() as usize).min(len - 2)
        }
        else
        {
            (scaled.floor() as usize).min(len)
        };
        let right_weight = scaled - left_lattice as f64;
        let left_weight = 1.0 - right_weight;

        let left = if self.storage.has_boundary()
        {
            Some(self.interpolation_order_1d[left_lattice])
        }
        else if left_lattice == 0
        {
            None
        }
        else
        {
            Some(self.interpolation_order_1d[left_lattice - 1])
        };

        let right = if self.storage.has_boundary()
        {
            Some(self.interpolation_order_1d[left_lattice + 1])
        }
        else if left_lattice >= len
        {
            None
        }
        else
        {
            Some(self.interpolation_order_1d[left_lattice])
        };

        if let Some(left_seq) = left
        {
            for output in 0..DIM_OUT
            {
                result[output] += self.values[left_seq][output] * left_weight;
            }
        }
        if let Some(right_seq) = right
        {
            for output in 0..DIM_OUT
            {
                result[output] += self.values[right_seq][output] * right_weight;
            }
        }
        result
    }

    fn interpolate_1d_inner_unit(&self, x: f64) -> [f64; DIM_OUT]
    {
        const MAX_LEVEL: u32 = 31;
        let bits = std::mem::size_of::<u32>() * 8;
        let val = (x * (1 << (bits - 2)) as f64).floor() * 2.0;
        let source = if x == 1.0 { (val - 1.0) as u32 } else { (val + 1.0) as u32 };
        let mut result = [0.0; DIM_OUT];
        let mut seq = 0;
        let mut level = 1;
        loop
        {
            let index = self.storage[seq].index[0];
            self.add_1d_contribution(seq, Self::linear_basis_value(level, index, x), &mut result);
            if self.storage[seq].is_leaf()
            {
                break;
            }
            let go_right = (source & (1 << (MAX_LEVEL - level))) > 0;
            level += 1;
            if let Some(next) = self.descend_1d(seq, go_right)
            {
                seq = next;
            }
            else
            {
                break;
            }
        }
        result
    }

    fn interpolate_1d_boundary_unit(&self, x: f64) -> [f64; DIM_OUT]
    {
        let mut result = [0.0; DIM_OUT];
        let mut seq = self.storage.adjacency_data.zero_index;
        if seq == usize::MAX
        {
            return result;
        }
        let left = self.storage.adjacency_data.left_zero[seq];
        if left != u32::MAX
        {
            seq = left as usize;
            self.add_1d_contribution(seq, Self::linear_basis_value(0, 0, x), &mut result);
        }
        let right = self.storage.adjacency_data.right_zero[seq];
        if right != u32::MAX
        {
            seq = right as usize;
            self.add_1d_contribution(seq, Self::linear_basis_value(0, 1, x), &mut result);
        }
        let mut level = 0;
        loop
        {
            if self.storage[seq].is_leaf()
            {
                break;
            }
            if level > 0
            {
                let index = self.storage[seq].index[0];
                let xh = x * POW2_F64[level as usize];
                let node = index as f64;
                if (xh - node).abs() < 1e-15
                {
                    break;
                }
                if let Some(next) = self.descend_1d(seq, xh > node)
                {
                    seq = next;
                }
                else
                {
                    break;
                }
            }
            else
            {
                if x.abs() < 1e-15 || (x - 1.0).abs() < 1e-15
                {
                    break;
                }
                let level_one = self.storage.adjacency_data[seq].level_one();
                if level_one == u32::MAX
                {
                    break;
                }
                seq = level_one as usize;
            }
            level += 1;
            let index = self.storage[seq].index[0];
            self.add_1d_contribution(seq, Self::linear_basis_value(level, index, x), &mut result);
        }
        result
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
            self.update_1d_interpolation_data();
            Ok(())
        }
    }
    pub fn coarsen<OP: HierarchisationOperation<D, DIM_OUT>, F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, op: &OP, update_iterator_data: bool, threshold: f64, remove_boundary: bool) -> Result<usize, SGError>
    {
        let mut total_num_removed = 0;
        let mut last_num_removed: usize;
        loop
        {            
            last_num_removed = self.coarsen_iteration(functor, threshold, remove_boundary);
            total_num_removed += last_num_removed;
            if last_num_removed > 0
            {
                self.hierarchize(op)?;
            }
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
        self.update_1d_interpolation_data();
        Ok(total_num_removed)
        
    }

    fn coarsen_iteration<F: RefinementFunctor<D, DIM_OUT>>(&mut self, functor: &F, threshold: f64, remove_boundary: bool) -> usize
    {
        let mut total_num_removed = 0;

        let r = algorithms::coarsening::coarsen(&mut self.storage, functor, &self.alpha, &self.values, threshold, remove_boundary);
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
            // Don't remove boundary nodes during iterative refine/coarsen - they may be needed by future refinements
            let removed = self.coarsen_iteration(functor, options.threshold, false);
            if removed > 0
            {
                if self.storage.has_boundary()
                {
                    let op = crate::const_generic::algorithms::hierarchisation::LinearBoundaryHierarchisationOperation;
                    self.hierarchize(&op).expect("Could not hierarchize grid after anisotropic coarsening.");
                }
                else
                {
                    let op = crate::const_generic::algorithms::hierarchisation::LinearHierarchisationOperation;
                    self.hierarchize(&op).expect("Could not hierarchize grid after anisotropic coarsening.");
                }
            }
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
        self.update_1d_interpolation_data();
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
                let removed = self.coarsen_iteration(functor, options.threshold, true);
                if removed > 0
                {
                    self.hierarchize(op).expect("Could not hierarchize grid after anisotropic coarsening.");
                }
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
