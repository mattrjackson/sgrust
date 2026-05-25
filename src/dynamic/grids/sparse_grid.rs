use std::io::Write;



//use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
use serde_with::serde_as;
use crate::dynamic::algorithms::basis_evaluation::BasisEvaluation;
use crate::dynamic::algorithms::refinement::{BaseRefinement, RefinementFunctor, RefinementOptions};
use crate::basis::linear::LinearBasis;
use crate::basis::linear::POW2_F64;
use crate::errors::SGError;
use crate::dynamic::algorithms::hierarchisation::HierarchisationOperation;
use crate::dynamic::iterators::grid_iterator_cache::AdjacencyGridIterator;
use crate::dynamic::storage::{BoundingBox, GridPoint, PointIterator, SparseGridData};
use crate::dynamic::generators::*;
use serde::{Serialize,Deserialize};
#[cfg(feature = "rkyv")]
use rkyv::with::Skip;
#[serde_as]
#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
pub struct SparseGridBase
{
    pub(crate) storage: SparseGridData,
    pub(crate) alpha: Vec<f64>,
    pub(crate) values: Vec<f64>,
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


impl SparseGridBase
{
    pub fn new(num_inputs: usize, num_outputs: usize) -> Self
    {
        SparseGridBase {
            storage: SparseGridData::new(num_inputs, num_outputs),
            alpha: Vec::new(),
            values: Vec::new(),
            interpolation_nodes_1d: Vec::new(),
            interpolation_order_1d: Vec::new(),
            interpolation_uniform_1d: false,
            interpolation_step_1d: 0.0,
        }
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

    pub fn sparse_grid<GENERATOR: Generator>(&mut self, levels: &[usize], generator: &GENERATOR) -> Result<(), SGError>
    {        
        generator.regular(&mut self.storage, levels, None)?;       
        self.sort();
        Ok(())
    }    

    pub fn full_grid<GENERATOR: Generator>(&mut self, level: usize, generator: &GENERATOR) -> Result<(), SGError>
    {
        generator.full(&mut self.storage, level)?;
        self.sort();
        Ok(())
    }

    pub fn sparse_grid_with_boundaries<GENERATOR: Generator>(&mut self, levels: &[usize], generator: &GENERATOR) -> Result<(), SGError>
    {
        generator.regular_with_boundaries(&mut self.storage, levels, Some(1), None)?;
        self.storage.has_boundary = true;
        self.sort();
        Ok(())
    }    

    pub fn full_grid_with_boundaries<GENERATOR: Generator>(&mut self, level: usize, generator: &GENERATOR) -> Result<(), SGError>
    {
        generator.full_with_boundaries(&mut self.storage, level)?;
        self.storage.has_boundary = true;
        self.sort();
        Ok(())
    }
    pub fn hierarchize<OP: HierarchisationOperation>(&mut self, op: &OP) -> Result<(), SGError>
    {
        self.alpha.clone_from(&self.values);
        op.hierarchize(&mut self.alpha, &self.storage)?;     
        Ok(())   
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
            *value = indices_rev[*value as usize] as u32;
        }
        self.update_1d_interpolation_data();
    }

    pub(crate) fn update_1d_interpolation_data(&mut self)
    {

        if self.storage.num_inputs > 1
        {
            return;
        }
        self.interpolation_nodes_1d.clear();
        self.interpolation_order_1d.clear();
        self.interpolation_uniform_1d = false;
        self.interpolation_step_1d = 0.0;


        let mut nodes = Vec::with_capacity(self.storage.len());
        for seq in 0..self.storage.len()
        {
            let coord = self.storage.index[seq] as f64 / POW2_F64[self.storage.level[seq] as usize];
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
    pub fn interpolate(&self, x: &[f64], result: &mut [f64] ) -> Result<(), SGError>
    {    
        if !self.bounding_box().contains(x)
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
        if self.storage.num_inputs == 1
        {
            return self.interpolate_1d_unchecked(x[0], result);
        }
        use crate::dynamic::algorithms::interpolation::InterpolationOperation;
        let iterator = &mut AdjacencyGridIterator::new(&self.storage);
        let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, self.storage.num_inputs, self.storage.num_outputs));      
        op.interpolate(x, &self.alpha, iterator, result)
    }

    #[inline]
    fn interpolate_1d_unchecked(&self, x: f64, result: &mut [f64]) -> Result<(), SGError>
    {
        if self.storage.is_empty()
        {
            return Ok(());
        }
        let bbox = self.storage.bounding_box();
        let x = (x - bbox.lower[0]) / (bbox.upper[0] - bbox.lower[0]);
        if self.interpolation_order_1d.len() == self.storage.len()
            && self.interpolation_nodes_1d.len() == self.storage.len()
            && self.values.len() == self.storage.len() * self.storage.num_outputs
        {
            self.interpolate_1d_sorted_unit(x, result);
            return Ok(());
        }
        if self.storage.has_boundary()
        {
            self.interpolate_1d_boundary_unit(x, result)
        }
        else
        {
            self.interpolate_1d_inner_unit(x, result)
        }
    }

    fn overwrite_1d_value(&self, seq: usize, result: &mut [f64])
    {
        let offset = seq * self.storage.num_outputs;
        result[..self.storage.num_outputs].copy_from_slice(&self.values[offset..offset + self.storage.num_outputs]);
    }

    fn write_1d_blend(&self, left: Option<usize>, left_weight: f64, right: Option<usize>, right_weight: f64, result: &mut [f64])
    {
        let left_offset = left.map(|seq| seq * self.storage.num_outputs);
        let right_offset = right.map(|seq| seq * self.storage.num_outputs);
        for output in 0..self.storage.num_outputs
        {
            let mut value = 0.0;
            if let Some(offset) = left_offset
            {
                value += self.values[offset + output] * left_weight;
            }
            if let Some(offset) = right_offset
            {
                value += self.values[offset + output] * right_weight;
            }
            result[output] = value;
        }
    }

    fn interpolate_1d_sorted_unit(&self, x: f64, result: &mut [f64])
    {
        if self.interpolation_uniform_1d
        {
            self.interpolate_1d_uniform_unit(x, result);
            return;
        }

        match self.interpolation_nodes_1d.binary_search_by(|node| node.total_cmp(&x))
        {
            Ok(pos) => self.overwrite_1d_value(self.interpolation_order_1d[pos], result),
            Err(pos) => {
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
                    (None, Some((_, right_seq))) => {
                        self.overwrite_1d_value(right_seq, result);
                        return;
                    }
                    (Some((_, left_seq)), None) => {
                        self.overwrite_1d_value(left_seq, result);
                        return;
                    }
                    (None, None) => {
                        result[..self.storage.num_outputs].fill(0.0);
                        return;
                    }
                };

                let width = right_x - left_x;
                if width <= 0.0
                {
                    result[..self.storage.num_outputs].fill(0.0);
                    return;
                }
                let right_weight = (x - left_x) / width;
                let left_weight = 1.0 - right_weight;

                self.write_1d_blend(left.map(|(_, seq)| seq), left_weight, right.map(|(_, seq)| seq), right_weight, result);
            }
        }
    }

    fn interpolate_1d_uniform_unit(&self, x: f64, result: &mut [f64])
    {
        let len = self.interpolation_order_1d.len();

        if self.storage.has_boundary()
        {
            if x <= 0.0
            {
                self.overwrite_1d_value(self.interpolation_order_1d[0], result);
                return;
            }
            if x >= 1.0
            {
                self.overwrite_1d_value(self.interpolation_order_1d[len - 1], result);
                return;
            }
        }
        else if x <= 0.0 || x >= 1.0
        {
            result[..self.storage.num_outputs].fill(0.0);
            return;
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

        self.write_1d_blend(left, left_weight, right, right_weight, result);
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
    fn add_1d_contribution(&self, seq: usize, weight: f64, result: &mut [f64])
    {
        let offset = seq * self.storage.num_outputs;
        for output in 0..self.storage.num_outputs
        {
            result[output] += self.alpha[offset + output] * weight;
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

    fn interpolate_1d_inner_unit(&self, x: f64, result: &mut [f64]) -> Result<(), SGError>
    {
        const MAX_LEVEL: u32 = 31;
        let bits = std::mem::size_of::<u32>() * 8;
        let val = (x * (1 << (bits - 2)) as f64).floor() * 2.0;
        let source = if x == 1.0 { (val - 1.0) as u32 } else { (val + 1.0) as u32 };
        let mut seq = 0;
        let mut level = 1;
        loop
        {
            let index = self.storage.index[seq];
            self.add_1d_contribution(seq, Self::linear_basis_value(level, index, x), result);
            if self.storage.flags[seq].is_leaf()
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
        Ok(())
    }

    fn interpolate_1d_boundary_unit(&self, x: f64, result: &mut [f64]) -> Result<(), SGError>
    {
        let mut seq = self.storage.adjacency_data.zero_index;
        if seq == usize::MAX
        {
            return Ok(());
        }
        let left = self.storage.adjacency_data.left_zero[seq];
        if left != u32::MAX
        {
            seq = left as usize;
            self.add_1d_contribution(seq, Self::linear_basis_value(0, 0, x), result);
        }
        let right = self.storage.adjacency_data.right_zero[seq];
        if right != u32::MAX
        {
            seq = right as usize;
            self.add_1d_contribution(seq, Self::linear_basis_value(0, 1, x), result);
        }
        let mut level = 0;
        loop
        {
            if self.storage.flags[seq].is_leaf()
            {
                break;
            }
            if level > 0
            {
                let index = self.storage.index[seq];
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
            let index = self.storage.index[seq];
            self.add_1d_contribution(seq, Self::linear_basis_value(level, index, x), result);
        }
        Ok(())
    }

    /// Create a reusable interpolation state for zero-allocation repeated interpolation.
    /// 
    /// Use with `interpolate_with_state` for maximum performance when interpolating many points.
    #[inline]
    pub fn create_interpolation_state(&self) -> crate::dynamic::algorithms::interpolation_state::InterpolationState<'_> {
        crate::dynamic::algorithms::interpolation_state::InterpolationState::new(&self.storage)
    }

    /// Interpolate using a pre-created state to avoid per-call overhead.
    /// 
    /// This is the fastest interpolation method when processing many points.
    #[inline]
    pub fn interpolate_with_state(&self, x: &[f64], state: &mut crate::dynamic::algorithms::interpolation_state::InterpolationState<'_>, result: &mut [f64]) -> Result<(), SGError>
    {
        if self.storage.num_inputs == 1
        {
            if !self.bounding_box().contains(x)
            {
                return Err(SGError::OutOfDomain);
            }
            return self.interpolate_1d_unchecked(x[0], result);
        }
        use crate::dynamic::algorithms::interpolation::InterpolationOperation;        
        let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, self.storage.num_inputs, self.storage.num_outputs));      
        op.interpolate(x, &self.alpha, &mut state.iterator, result)
    }

    #[cfg(feature="rayon")]
    #[inline]
    pub fn interpolate_batch(&self, x: &[f64]) -> Result<Vec<f64>, SGError>
    {
       use rayon::{iter::{IndexedParallelIterator, ParallelIterator}, slice::{ParallelSlice, ParallelSliceMut}};
       use crate::dynamic::algorithms::interpolation::InterpolationOperation;
       let mut results = vec![0.0; x.len() / self.storage.num_inputs *self.storage.num_outputs];       
       let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, self.storage.num_inputs, self.storage.num_outputs));
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
        use rayon::{iter::{IndexedParallelIterator, ParallelIterator}, slice::{ParallelSlice, ParallelSliceMut}};
        use crate::dynamic::algorithms::interpolation::InterpolationOperation;
        let mut results = vec![0.0; x.len() / self.storage.num_inputs *self.storage.num_outputs];       
        let op = InterpolationOperation(self.storage.has_boundary(), BasisEvaluation(&self.storage, self.storage.num_inputs, self.storage.num_outputs));        
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
            self.update_1d_interpolation_data();
            Ok(())
        }
    }
    pub fn coarsen<OP: HierarchisationOperation, F: RefinementFunctor>(&mut self, functor: &F, op: &OP, update_iterator_data: bool, threshold: f64, remove_boundary: bool) -> Result<usize, SGError>
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

    fn coarsen_iteration<F: RefinementFunctor>(&mut self, functor: &F, threshold: f64, remove_boundary: bool) -> usize
    {
        let mut total_num_removed = 0;
  
        let r = crate::dynamic::algorithms::coarsening::coarsen(&mut self.storage, functor, &self.alpha, &self.values, threshold, remove_boundary);
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
        if options.refinement_mode == crate::dynamic::algorithms::refinement::RefinementMode::Anisotropic
        {
            let removed = self.coarsen_iteration(functor, options.threshold, false);
            if removed > 0
            {
                if self.storage.has_boundary()
                {
                    let op = crate::dynamic::algorithms::hierarchisation::LinearBoundaryHierarchisationOperation;
                    self.hierarchize(&op).expect("Could not hierarchize grid after anisotropic coarsening.");
                }
                else
                {
                    let op = crate::dynamic::algorithms::hierarchisation::LinearHierarchisationOperation;
                    self.hierarchize(&op).expect("Could not hierarchize grid after anisotropic coarsening.");
                }
            }
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
        self.update_1d_interpolation_data();
        points
    }
    pub fn refine<F: RefinementFunctor, OP: HierarchisationOperation, EF: Fn(&[f64])->Vec<f64>>(&mut self, functor: &F, eval_fun: &EF, op: &OP, options: RefinementOptions, max_iterations: usize) -> Result<(), SGError>
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
                let mut point = self.storage.unit_coordinate(i);
                point = self.storage.bounding_box().to_real_coordinate(&point);
                self.values.extend(eval_fun(&point));
            }            
            self.hierarchize(op)?;
            iteration += 1;            
   
            if iteration == max_iterations
            {
                break;
            }
        }
        self.sort(); 
        Ok(())       
    }

    #[cfg(feature="rayon")]
    pub fn refine_parallel<F: RefinementFunctor, OP: HierarchisationOperation, EF: Fn(&[f64]) -> Vec<f64> + Send + Sync>(&mut self, functor: &F, 
        eval_fun: &EF, op: &OP, options: RefinementOptions, max_iterations: usize) 
    {
        use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
        let ref_op = BaseRefinement(self.storage.has_boundary());
        let mut iteration = 0;     
        loop
        {
            use rayon::slice::ParallelSliceMut;

            if iteration > 0 && options.refinement_mode == crate::dynamic::algorithms::refinement::RefinementMode::Anisotropic
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
            self.hierarchize(op).expect("Could not hiearchize grid.");
            iteration += 1;            
            if max_iterations <= iteration
            {
                break;
            }
        }
        self.sort();        
    }

  
    ///
    /// Writes full grid information (including maps and iterator data).
    /// All formats use LZ4 compression.
    /// 
    pub fn write(&mut self, path: &str, format: crate::serialization::SerializationFormat) -> Result<(), SGError>
    {
        let mut file = std::io::BufWriter::new(std::fs::File::create(path).map_err(|_|SGError::FileIOError)?);        
        let buffer = crate::serialization::serialize(self, format)?;
        file.write_all(&buffer).map_err(|_|SGError::WriteBufferFailed)?;
        Ok(())
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
