use crate::{
    basis::base::Basis,
    basis::global::{GlobalBasis, GlobalBasisType},
    dynamic::{
        algorithms::lagrange::{lagrange_coeffs, lagrange_weights},
        iterators::advanced_step_iterator::{AdvancedStepIterator, GridType},
    },
    errors::SGError,
};
use ndarray::{ArrayView1, ArrayView2};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

const NODE_KEY_TOLERANCE: f64 = 1e-12;

fn cartesian_product(bounds: &[u32]) -> Vec<u32> {
    let ndim = bounds.len();
    let total_combinations = bounds.iter().product::<u32>() as usize;
    let mut multi_indices = vec![0; total_combinations * ndim];

    for (i, current) in multi_indices.chunks_exact_mut(ndim).enumerate() {
        let mut index = i as u32;
        for (j, &bound) in bounds.iter().rev().enumerate() {
            current[ndim - 1 - j] = index % bound;
            index /= bound;
        }
    }

    multi_indices
}

fn point_key(point: &[f64]) -> Vec<i64> {
    point
        .iter()
        .map(|&x| (x / NODE_KEY_TOLERANCE).round() as i64)
        .collect()
}

#[derive(Default, Serialize, Deserialize, Clone, Copy)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub enum TensorSelectionStrategy
{    
    /// Tensors selected based on sum of levels across the dimensions.
    #[default]
    Level,    
    /// Tensors selected based on exactness of quadrature rule across the dimensions.
    QuadratureExactness,
    /// Tensors selected based on exactness of interpolation rule across the dimensions.
    InterpolationExactness,
}

#[derive(Default, Serialize, Deserialize, Clone)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct GenerationOptions
{
    /// Tensor Selection Strategy
    pub tensor_strategy: TensorSelectionStrategy,
    /// Level limits in each dimension
    pub level_limits: Vec<u32>,
    /// Total polynomial exactness (only used if TensorSelectionStrategy is `Exactness`)
    pub exactness_limit: Option<u32>,
}

#[derive(Default, Serialize, Deserialize, Clone)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct DynamicBoundingBox {
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
}

impl DynamicBoundingBox {
    pub fn new(lower: &[f64], upper: &[f64]) -> Self {
        Self {
            lower: lower.to_owned(),
            upper: upper.to_owned(),
        }
    }

    #[inline]
    pub fn ndim(&self) -> usize {
        self.lower.len()
    }

    #[inline]
    pub fn width(&self, dim: usize) -> f64 {
        self.upper[dim] - self.lower[dim]
    }

    ///
    /// Volume of hypercube (width(dim1)*...*width(dim_n))
    /// 
    pub fn volume(&self) -> f64 {
        let mut volume = 1.0;
        for d in 0..self.ndim() {
            volume *= self.width(d);
        }
        volume
    }

    pub fn to_unit_coordinate(&self, point: &[f64]) -> Vec<f64> {
        let mut r = vec![0.0; self.ndim()];
        for i in 0..self.ndim() {
            r[i] = (point[i] - self.lower[i]) / (self.upper[i] - self.lower[i]);
        }
        r
    }

    pub fn to_real_coordinate(&self, point: &mut [f64]) {
        for i in 0..self.ndim() {
            point[i] = self.lower[i] + (self.upper[i] - self.lower[i]) * point[i];
        }
    }

    pub fn contains(&self, point: &[f64]) -> bool {
        for d in 0..self.ndim() {
            if self.lower[d] > point[d] || self.upper[d] < point[d] {
                return false;
            }
        }
        true
    }
}

#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct NodesAndCoefficientsForLevel {
    pub level: u32,
    pub nodes: Vec<f64>,
    pub coefficients: Vec<f64>,
}

#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct FullGrid {
    pub level: Vec<u32>,
    pub lweights: Vec<NodesAndCoefficientsForLevel>,
    pub weight: f64,
    index_map: Vec<u32>,
    combinations: Vec<u32>,
}

impl FullGrid {
    ///
    /// Number of dimensions
    /// 
    pub fn ndim(&self) -> usize {
        self.level.len()
    }

    ///
    /// Returns the point indices for the `basis` functions. 
    /// 
    pub fn points(&self, basis: &[GlobalBasis]) -> Vec<u32> {
        let num_points: Vec<u32> = self
            .level
            .iter()
            .zip(basis)
            .map(|(&level, basis)| basis.num_nodes(level) as u32)
            .collect();
        cartesian_product(&num_points)
    }

    ///
    /// Create a new `FullGrid` for a given `level`, using the specified `basis` functions, 
    /// and the global `weight`for the grid.
    /// 
    pub fn new(basis: &[GlobalBasis], level: &[u32], weight: f64) -> Self {
        // Determine the number of points in each dimension from the basis functions
        let npts: Vec<u32> = level
            .iter()
            .zip(basis)
            .map(|(&level, basis)| basis.num_nodes(level) as u32)
            .collect();
        // we store combinations in order to speed up interpolation later...
        let combinations = Self::generate_combinations(&npts);
        let mut grid = Self {
            level: level.to_owned(),
            lweights: Vec::new(),
            weight,
            index_map: Vec::new(),
            combinations,
        };

        // Compute lagrange coefficients for tensor
        let mut coeffs = Vec::new();
        for d in 0..level.len() {
            let nodes = basis[d].nodes(level[d]);
            coeffs.push(NodesAndCoefficientsForLevel {
                coefficients: lagrange_coeffs(&nodes),
                level: level[d],
                nodes,
            });
        }
        grid.lweights = coeffs;
        grid
    }

    #[inline]
    fn lagrange_weight_size(&self) -> usize {
        let mut size = 1;
        for item in &self.lweights {
            size *= item.nodes.len();
        }
        size
    }

    fn generate_combinations(npts: &[u32]) -> Vec<u32> {
        let mut combinations = Vec::new();
        let total_combinations = npts.iter().product::<u32>() as usize;
        let ndim = npts.len();
        let mut current = vec![0; ndim];
        // Iterate through all possible combinations
        for i in 0..total_combinations {
            let mut index = i as u32;
            for (j, &bound) in npts.iter().rev().enumerate() {
                current[ndim - 1 - j] = index % bound;
                index /= bound;
            }
            combinations.extend(&current);
            current.fill(0);
        }

        combinations
    }

    ///
    /// Computes the lagrange weights for a given `x` coordinate.
    /// Returns vector starting indices for each dimension and a vector 
    /// containing weights across all dimensions. 
    /// 
    fn lagrange_weights(&self, x: &[f64]) -> (Vec<usize>, Vec<f64>) {
        let mut weights = Vec::with_capacity(self.lagrange_weight_size());
        let mut indices = vec![0; self.ndim()];
        for dim in 0..self.ndim() {
            indices[dim] = weights.len();
            weights.extend(lagrange_weights(
                x[dim],
                &self.lweights[dim].coefficients,
                &self.lweights[dim].nodes,
            ));
        }
        (indices, weights)
    }
    ///
    /// Interpolate at `x` along a given grid using the specified `basis` functions. 
    /// 
    #[inline]
    pub fn interpolate(&self, x: &[f64], y: &mut [f64], values: &[f64]) {
        // Retrieve the lagrange weights for this coordinate
        let (indices, weights) = self.lagrange_weights(x);
        // Iterate through all possible combinations
        for (i, combination) in self.combinations.chunks_exact(self.ndim()).enumerate() {
            let mut combined_weight = 1.0;
            for dim in 0..self.ndim() {
                combined_weight *= weights[indices[dim] + combination[dim] as usize];
            }
            // Retrive the value(s) for this node by using the `index_map` to determine 
            // the global index from the local grid index
            let values = &values[y.len() * self.index_map[i] as usize..];
            // Add contribution in each dimension
            for output in 0..y.len() {
                y[output] += combined_weight * values[output] * self.weight;
            }
        }
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub enum TensorIndicatorNorm {
    #[default]
    MaxAbs,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct AdaptiveRefinementOptions {
    pub threshold: f64,
    pub level_limits: Vec<u32>,
    pub max_new_tensors: Option<usize>,
    pub indicator_norm: TensorIndicatorNorm,
}

impl AdaptiveRefinementOptions {
    pub fn new(threshold: f64, level_limits: Vec<u32>) -> Self {
        Self {
            threshold,
            level_limits,
            max_new_tensors: None,
            indicator_norm: TensorIndicatorNorm::MaxAbs,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct TensorCandidate {
    pub level: Vec<u32>,
    pub parent_level: Vec<u32>,
    pub dimension: usize,
    pub node_start: usize,
    pub node_len: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct RefinementProposal {
    pub candidates: Vec<TensorCandidate>,
    pub candidate_nodes: Vec<f64>,
    pub predicted_values: Vec<f64>,
    pub num_outputs: usize,
    threshold: f64,
    max_new_tensors: Option<usize>,
    indicator_norm: TensorIndicatorNorm,
    candidate_unit_nodes: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct RefinementUpdate {
    pub accepted_tensors: Vec<Vec<u32>>,
    pub accepted_scores: Vec<f64>,
    pub new_nodes: Vec<f64>,
    pub new_values: Vec<f64>,
}

///
/// Generates a sparse grid using the combination technique. Unlike the 
/// `LinearSparseGrid`, this grid is not presently adaptively refineable. 
/// Instead, it is primarily optimized for uncertainty quantification 
/// operations (roughly equivalent to the `GlobalGrid` concept found in
/// TASMANIAN).
/// 
#[derive(Clone, Default, Serialize, Deserialize)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct CombinationSparseGrid {
    bounding_box: Option<DynamicBoundingBox>,
    ndim: usize,
    basis: Vec<GlobalBasis>,
    grids: Vec<FullGrid>,
    nodes: Vec<f64>,
    /// a hash of grid levels we already have...
    grid_hash: FxHashSet<Vec<u32>>,
    index_set: Vec<Vec<u32>>,
    index_hash: FxHashSet<Vec<u32>>,
    /// quadrature weight for node...
    qweight: Vec<f64>,
}

impl CombinationSparseGrid {
    ///
    /// Create an empty Grid with dimension `ndim` and
    /// basis functions `basis`.
    /// 
    pub fn new(ndim: usize, basis: Vec<GlobalBasis>) -> Self {
        assert!(basis.len() == ndim);
        Self {
            ndim,
            basis,
            ..Default::default()
        }
    }

    pub fn len(&self) -> usize
    {
        if self.ndim == 0
        {
            0
        }
        else
        {
            self.nodes.len() / self.ndim
        }
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    ///
    /// Number of dimensions
    /// 
    pub fn ndim(&self) -> usize {
        self.ndim
    }

    ///
    /// Compute the quadrature weight for a given `level` and `index` using 
    /// the user-specified `basis` functions.
    /// 
    pub(crate) fn quadrature_weight(basis: &[GlobalBasis], index: &[u32], level: &[u32]) -> f64 {
        let mut weight = 1.0;
        for d in 0..basis.len() {
            weight *= basis[d].integral(level[d], index[d]);
        }
        weight
    }

    ///
    /// Compute the floating point coordinates for a given `level` and `index`
    /// 
    fn unit_coordinate(basis: &[GlobalBasis], level: &[u32], index: &[u32]) -> Vec<f64> {
        let mut coor = vec![0.0; basis.len()];
        for d in 0..basis.len() {
            coor[d] = basis[d].get_point(level[d], index[d]);
        }
        coor
    }

    fn supports_adaptive_refinement(&self) -> bool {
        self.basis.iter().all(|basis| {
            matches!(
                basis.basis_type,
                GlobalBasisType::ClenshawCurtis | GlobalBasisType::GaussPatterson
            )
        })
    }

    fn validate_values(&self, values: &[f64], num_outputs: usize) -> Result<(), SGError> {
        if values.len() != self.len() * num_outputs {
            return Err(SGError::NumberOfPointsAndValuesMismatch);
        }
        Ok(())
    }

    fn validate_level_limits(&self, level_limits: &[u32]) -> Result<(), SGError> {
        if level_limits.len() != self.ndim {
            return Err(SGError::InvalidIndex);
        }
        Ok(())
    }

    fn index_set_from_options(
        &self,
        grid_type: GridType,
        generation_parameters: &GenerationOptions,
    ) -> Vec<Vec<u32>> {
        let exactness_bound = match generation_parameters.tensor_strategy {
            TensorSelectionStrategy::Level => u32::MAX,
            TensorSelectionStrategy::QuadratureExactness => {
                let max_exactness = self
                    .basis
                    .iter()
                    .zip(&generation_parameters.level_limits)
                    .max_by(|a, b| {
                        a.0.quadrature_exactness(a.1 - 1)
                            .cmp(&b.0.quadrature_exactness(b.1 - 1))
                    })
                    .unwrap();
                max_exactness.0.quadrature_exactness(max_exactness.1 - 1) + 1
            }
            TensorSelectionStrategy::InterpolationExactness => {
                let max_exactness = self
                    .basis
                    .iter()
                    .zip(&generation_parameters.level_limits)
                    .max_by(|a, b| {
                        a.0.interpolation_exactness(a.1 - 1)
                            .cmp(&b.0.interpolation_exactness(b.1 - 1))
                    })
                    .unwrap();
                max_exactness.0.interpolation_exactness(max_exactness.1 - 1) + 1
            }
        };

        if self.ndim == 1 {
            return vec![generation_parameters.level_limits.clone()];
        }

        AdvancedStepIterator::new(
            &generation_parameters.level_limits,
            grid_type,
            self.basis.clone(),
            generation_parameters.tensor_strategy,
            generation_parameters
                .exactness_limit
                .unwrap_or(exactness_bound),
        )
        .collect()
    }

    fn add_tensor_to_index_set(&mut self, level: Vec<u32>) {
        if self.index_hash.insert(level.clone()) {
            self.index_set.push(level);
        }
    }

    fn node_map(&self) -> FxHashMap<Vec<i64>, usize> {
        let mut map = FxHashMap::default();
        for (idx, node) in self.nodes.chunks_exact(self.ndim).enumerate() {
            map.insert(point_key(node), idx);
        }
        map
    }

    fn rebuild_from_index_set(&mut self) -> Result<(), SGError> {
        if self.index_set.is_empty() {
            self.grids.clear();
            self.grid_hash.clear();
            self.qweight.clear();
            self.nodes.clear();
            return Ok(());
        }

        self.index_set.sort();
        let flattened: Vec<u32> = self
            .index_set
            .iter()
            .flat_map(|level| level.iter().copied())
            .collect();
        let weights =
            crate::utilities::multi_index_manipulation::weight_modifiers(&flattened, self.ndim)?;
        let mut node_lookup = self.node_map();
        let mut grids = Vec::new();
        let mut grid_hash = FxHashSet::default();
        let mut qweight = vec![0.0; self.len()];

        for (level_set, weight) in self.index_set.iter().zip(weights) {
            if weight == 0.0 {
                continue;
            }

            let mut grid = FullGrid::new(&self.basis, level_set, weight);
            for index in grid.points(&self.basis).chunks_exact(self.ndim) {
                let point = Self::unit_coordinate(&self.basis, &grid.level, index);
                let key = point_key(&point);
                let idx = if let Some(&idx) = node_lookup.get(&key) {
                    idx
                } else {
                    let idx = self.nodes.len() / self.ndim;
                    self.nodes.extend(&point);
                    qweight.push(0.0);
                    node_lookup.insert(key, idx);
                    idx
                };

                qweight[idx] += weight * Self::quadrature_weight(&self.basis, index, level_set);
                grid.index_map.push(idx as u32);
            }

            grid_hash.insert(level_set.clone());
            grids.push(grid);
        }

        self.grids = grids;
        self.grid_hash = grid_hash;
        self.qweight = qweight;
        Ok(())
    }

    fn real_nodes_from_unit(&self, unit_nodes: &[f64]) -> Vec<f64> {
        let mut nodes = unit_nodes.to_vec();
        if let Some(bbox) = &self.bounding_box {
            for point in nodes.chunks_exact_mut(self.ndim) {
                bbox.to_real_coordinate(point);
            }
        }
        nodes
    }

    fn interpolate_unit(&self, x: &[f64], values: &[f64], num_outputs: usize) -> Vec<f64> {
        let mut y = vec![0.0; num_outputs];
        if self.ndim == 1 {
            if let Some(grid) = self.grids.last() {
                let lweights = grid.lweights.last().unwrap();
                let weights = lagrange_weights(
                    x[0],
                    &lweights.coefficients,
                    &self.basis[0].nodes(grid.level[0]),
                );
                for (&weight, value) in weights.iter().zip(values.chunks_exact(num_outputs)) {
                    for i in 0..num_outputs {
                        y[i] += weight * value[i];
                    }
                }
            }
        } else {
            for grid in &self.grids {
                grid.interpolate(x, &mut y, values);
            }
        }
        y
    }

    fn candidate_score(
        predicted: &[f64],
        actual: &[f64],
        indicator_norm: TensorIndicatorNorm,
    ) -> f64 {
        match indicator_norm {
            TensorIndicatorNorm::MaxAbs => predicted
                .iter()
                .zip(actual)
                .map(|(lhs, rhs)| (lhs - rhs).abs())
                .fold(0.0, f64::max),
        }
    }

    fn is_admissible_candidate(&self, level: &[u32]) -> bool {
        if self.ndim == 1 {
            return true;
        }

        for dim in 0..self.ndim {
            if level[dim] == 0 {
                continue;
            }
            let mut predecessor = level.to_vec();
            predecessor[dim] -= 1;
            if !self.index_hash.contains(&predecessor) {
                return false;
            }
        }

        true
    }

    ///
    /// Generate a standard sparse grid using the combination technique.
    /// 
    pub fn sparse_grid(&mut self, generation_parameters: GenerationOptions) -> Result<(), SGError> {
        self.generate_grid(GridType::Sparse, generation_parameters)
    }
    ///
    /// Compute a full grid using the combination technique
    /// 
    pub fn full_grid(&mut self, generation_parameters: GenerationOptions) -> Result<(), SGError> {
        self.generate_grid(GridType::Full, generation_parameters)
    }

    fn generate_grid(
        &mut self,
        grid_type: GridType,
        generation_parameters: GenerationOptions,
    ) -> Result<(), SGError> {
        self.validate_level_limits(&generation_parameters.level_limits)?;
        for level_set in self.index_set_from_options(grid_type, &generation_parameters) {
            self.add_tensor_to_index_set(level_set);
        }
        self.rebuild_from_index_set()
    }
    ///
    /// Return reference to bounding box
    /// 
    pub fn bounding_box(&self) -> &Option<DynamicBoundingBox> {
        &self.bounding_box
    }

    ///
    /// Return mutable reference for bounding box
    /// 
    pub fn bounding_box_mut(&mut self) -> &mut Option<DynamicBoundingBox> {
        &mut self.bounding_box
    }
    ///
    /// Return vector containing real coordinates for grid (size = `ndim` * nodes().len()).
    /// 
    pub fn nodes(&self) -> Vec<f64> {
        self.real_nodes_from_unit(&self.nodes)
    }

    ///
    /// Compute integral over grid
    /// 
    pub fn integral(&self, values: &[f64], num_outputs: usize) -> Vec<f64> {
        let mut y = vec![0.0; num_outputs];
        for (&weight, value) in self.qweight.iter().zip(values.chunks_exact(num_outputs)) {
            for i in 0..num_outputs {
                y[i] += weight * value[i];
            }
        }
        y
    }
    ///
    /// Compute integral over grid. Values must be in column major format...
    /// 
    pub fn integral_fast(&self, values: &[f64], num_outputs: usize) -> Vec<f64>
    {
        
        if num_outputs == 1
        {
            let weights = ArrayView1::from(&self.qweight);
            let values = ArrayView1::from(&values);    
            return vec![weights.dot(&values)];
        }
        let weights = ArrayView1::from(&self.qweight);
        let values = ArrayView2::from_shape((self.len(), num_outputs), values).unwrap();
           
        let product = weights.dot(&values);        
        let mut y = Vec::from(product.as_slice().unwrap());       
        if let Some(bbox) = &self.bounding_box
        {
            let volume = bbox.volume();
            for value in &mut y {
                *value *= volume;
            }
        }

        y
    }

    ///
    /// Interpolate on grid using previously computed `values`
    /// 
    pub fn interpolate(&self, x: &[f64], values: &[f64], num_outputs: usize) -> Vec<f64> {
        self.validate_values(values, num_outputs)
            .expect("values length must match number of nodes");
        let x = if let Some(bbox) = &self.bounding_box {
            bbox.to_unit_coordinate(x)
        } else {
            x.to_vec()
        };
        self.interpolate_unit(&x, values, num_outputs)
    }

    pub fn propose_refinement(
        &self,
        values: &[f64],
        num_outputs: usize,
        options: &AdaptiveRefinementOptions,
    ) -> Result<RefinementProposal, SGError> {
        self.validate_values(values, num_outputs)?;
        self.validate_level_limits(&options.level_limits)?;
        if !self.supports_adaptive_refinement() {
            return Err(SGError::OperationNotSupported);
        }

        let current_node_map = self.node_map();
        let mut seen = FxHashSet::default();
        let mut candidates = Vec::new();
        let mut candidate_unit_nodes = Vec::new();
        let mut predicted_values = Vec::new();
        let base_levels = if self.index_set.is_empty() {
            vec![vec![0; self.ndim]]
        } else {
            self.index_set.clone()
        };

        for base_level in &base_levels {
            for dim in 0..self.ndim {
                let mut level = base_level.clone();
                if level[dim] >= options.level_limits[dim] {
                    continue;
                }
                level[dim] += 1;
                if self.index_hash.contains(&level)
                    || !seen.insert(level.clone())
                    || !self.is_admissible_candidate(&level)
                {
                    continue;
                }

                let grid = FullGrid::new(&self.basis, &level, 1.0);
                let node_start = candidate_unit_nodes.len() / self.ndim;
                for index in grid.points(&self.basis).chunks_exact(self.ndim) {
                    let point = Self::unit_coordinate(&self.basis, &level, index);
                    if current_node_map.contains_key(&point_key(&point)) {
                        continue;
                    }
                    predicted_values.extend(self.interpolate_unit(&point, values, num_outputs));
                    candidate_unit_nodes.extend(&point);
                }

                let node_len = candidate_unit_nodes.len() / self.ndim - node_start;
                if node_len == 0 {
                    continue;
                }

                let mut parent_level = level.clone();
                parent_level[dim] -= 1;
                candidates.push(TensorCandidate {
                    level,
                    parent_level,
                    dimension: dim,
                    node_start,
                    node_len,
                });
            }
        }

        Ok(RefinementProposal {
            candidates,
            candidate_nodes: self.real_nodes_from_unit(&candidate_unit_nodes),
            predicted_values,
            num_outputs,
            threshold: options.threshold,
            max_new_tensors: options.max_new_tensors,
            indicator_norm: options.indicator_norm,
            candidate_unit_nodes,
        })
    }

    pub fn apply_refinement(
        &mut self,
        proposal: RefinementProposal,
        candidate_values: &[f64],
    ) -> Result<RefinementUpdate, SGError> {
        if candidate_values.len() != proposal.predicted_values.len() {
            return Err(SGError::NumberOfPointsAndValuesMismatch);
        }
        if !self.supports_adaptive_refinement() {
            return Err(SGError::OperationNotSupported);
        }

        let mut scored_candidates = Vec::new();
        for candidate in &proposal.candidates {
            let start = candidate.node_start * proposal.num_outputs;
            let end = start + candidate.node_len * proposal.num_outputs;
            let score = Self::candidate_score(
                &proposal.predicted_values[start..end],
                &candidate_values[start..end],
                proposal.indicator_norm,
            );
            if score >= proposal.threshold {
                scored_candidates.push((score, candidate.clone()));
            }
        }

        scored_candidates.sort_by(|lhs, rhs| {
            rhs.0
                .partial_cmp(&lhs.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        if let Some(limit) = proposal.max_new_tensors {
            scored_candidates.truncate(limit);
        }

        if scored_candidates.is_empty() {
            return Ok(RefinementUpdate {
                accepted_tensors: Vec::new(),
                accepted_scores: Vec::new(),
                new_nodes: Vec::new(),
                new_values: Vec::new(),
            });
        }

        let old_len = self.len();
        let mut new_value_map = FxHashMap::default();
        let mut accepted_tensors = Vec::with_capacity(scored_candidates.len());
        let mut accepted_scores = Vec::with_capacity(scored_candidates.len());

        for (score, candidate) in scored_candidates {
            accepted_scores.push(score);
            accepted_tensors.push(candidate.level.clone());
            self.add_tensor_to_index_set(candidate.level);

            for point_idx in candidate.node_start..candidate.node_start + candidate.node_len {
                let node = &proposal.candidate_unit_nodes
                    [point_idx * self.ndim..(point_idx + 1) * self.ndim];
                let start = point_idx * proposal.num_outputs;
                let end = start + proposal.num_outputs;
                new_value_map.insert(point_key(node), candidate_values[start..end].to_vec());
            }
        }

        self.rebuild_from_index_set()?;

        let mut new_unit_nodes = Vec::new();
        let mut new_values = Vec::new();
        for node in self.nodes.chunks_exact(self.ndim).skip(old_len) {
            let value = new_value_map
                .get(&point_key(node))
                .ok_or(SGError::InvalidIndex)?;
            new_unit_nodes.extend(node);
            new_values.extend(value);
        }

        Ok(RefinementUpdate {
            accepted_tensors,
            accepted_scores,
            new_nodes: self.real_nodes_from_unit(&new_unit_nodes),
            new_values,
        })
    }

    pub fn refine_with<EF: Fn(&[f64]) -> Vec<f64>>(
        &mut self,
        values: &mut Vec<f64>,
        num_outputs: usize,
        options: &AdaptiveRefinementOptions,
        eval: EF,
    ) -> Result<usize, SGError> {
        let proposal = self.propose_refinement(values, num_outputs, options)?;
        if proposal.candidate_nodes.is_empty() {
            return Ok(0);
        }

        let mut candidate_values =
            Vec::with_capacity(proposal.candidate_nodes.len() / self.ndim * num_outputs);
        for point in proposal.candidate_nodes.chunks_exact(self.ndim) {
            let value = eval(point);
            if value.len() != num_outputs {
                return Err(SGError::NumberOfPointsAndValuesMismatch);
            }
            candidate_values.extend(value);
        }

        let update = self.apply_refinement(proposal, &candidate_values)?;
        values.extend(&update.new_values);
        Ok(update.accepted_tensors.len())
    }
    ///
    /// Compute the product of two sparse grids.
    /// 
    pub fn product(
        &self,
        other: &CombinationSparseGrid,
        mut options: GenerationOptions,
    ) -> Result<CombinationSparseGrid, SGError> {
        let mut basis_vectors = Vec::new();
        basis_vectors.extend(self.basis.clone());
        basis_vectors.extend(other.basis.clone());

        let mut lower = Vec::new();
        let mut upper = Vec::new();
        if let Some(bbox) = &self.bounding_box {
            lower.extend(bbox.lower.clone());
            upper.extend(bbox.upper.clone());
            if let Some(bbox2) = &other.bounding_box {
                lower.extend(bbox2.lower.clone());
                upper.extend(bbox2.upper.clone());
            } else {
                lower.extend(vec![0.0; other.basis.len()]);
                upper.extend(vec![1.0; other.basis.len()]);
            }
        } else if let Some(bbox2) = &other.bounding_box {
             // Only need to create bbox if the second grid has a bounding box defined.
            lower.extend(vec![0.0; self.basis.len()]);
            lower.extend(bbox2.lower.clone());
            upper.extend(vec![1.0; self.basis.len()]);
            upper.extend(bbox2.upper.clone());
        }

        let mut grid = CombinationSparseGrid::new(self.ndim + other.ndim, basis_vectors);
        grid.bounding_box = if lower.is_empty() {
            None
        } else {
            Some(DynamicBoundingBox::new(&lower, &upper))
        };
        options.level_limits.clear();
        options.level_limits.extend(self.max_limits());
        options.level_limits.extend(other.max_limits());
        grid.sparse_grid(options)?;
        Ok(grid)
    }

    ///
    /// Copy values existing on previous grid into new one. 
    /// 
    pub fn copy_values_to_grid<T: Clone + Default>(
        &self,
        offset: usize,
        source: &CombinationSparseGrid,
        source_values: &[T],
    ) -> Vec<T> {
        let mut v = vec![T::default(); self.len()];
        let source_node_map = source.node_map();
        let map = self.node_map();
        for (node, &index) in &map {
            let source_idx = &node[offset..offset + source.basis.len()];
            if let Some(&idx) = source_node_map.get(source_idx) {
                v[index] = source_values[idx].clone();
            }
        }
        v
    }
    /// Get basis functions associated with this grid
    pub fn basis(&self) -> &Vec<GlobalBasis> {
        &self.basis
    }
    /// Get level limits for each dimension 
    pub fn max_limits(&self) -> Vec<u32> {
        let mut limits = vec![0; self.ndim];
        for level in &self.index_set {
            for (i, &value) in level.iter().enumerate() {
                limits[i] = limits[i].max(value);
            }
        }
        limits
    }
}

#[test]
fn test_combination_grid() {
    let mut grid = CombinationSparseGrid::new(
        2,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::ClenshawCurtis,
                custom_rule: None
            };
            2
        ],
    );
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::Level,
        level_limits: vec![3, 3],
        ..Default::default()
    };
    grid.full_grid(options.clone()).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / 2);
    for x in grid.nodes.chunks_exact(2) {
        values.push(x[0] * x[0] + x[1] * x[1]);
    }
    assert_eq!(
        grid.len(),
        options
            .level_limits
            .iter()
            .map(|&l| (1 << l) + 1)
            .product::<usize>()
    );
    println!("integral={}", grid.integral(&values, 1)[0]);
}

#[test]
fn test_combination_grid_gauss_legendre() {
    let mut grid = CombinationSparseGrid::new(
        2,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::GaussLegendre,
                custom_rule: None
            };
            2
        ],
    );
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::Level,
        level_limits: vec![10, 10],
        ..Default::default()
    };
    grid.full_grid(options).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / 2);
    for x in grid.nodes.chunks_exact(2) {
        values.push(x[0] * x[0] + x[1] * x[1]);
    }
    assert_eq!(grid.len(), 100);
    println!("integral={}", grid.integral(&values, 1)[0]);
}

#[test]
fn test_combination_grid_integral_standard() {
    let mut grid = CombinationSparseGrid::new(
        2,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::ClenshawCurtis,
                custom_rule: None
            };
            2
        ],
    );
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::Level,
        level_limits: vec![3, 3],
        ..Default::default()
    };
    grid.full_grid(options).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / 2);
    for x in grid.nodes.chunks_exact(2) {
        values.push(x[0] * x[0] + x[1] * x[1]);
    }
    println!("integral={}", grid.integral(&values, 1)[0]);
}

#[test]
fn test_combination_grid_integral_ndarray() {
    let mut grid = CombinationSparseGrid::new(
        2,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::ClenshawCurtis,
                custom_rule: None
            };
            2
        ],
    );
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::Level,
        level_limits: vec![3, 3],
        ..Default::default()
    };
    grid.full_grid(options).unwrap();
    let mut values = vec![0.0; grid.len()];
    for (i, x) in grid.nodes.chunks_exact(2).enumerate() {
        values[i] = x[0] * x[0] + x[1] * x[1];
    }
    println!("integral={}", grid.integral_fast(&values, 1)[0]);
}

#[test]
fn test_1d() {
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::InterpolationExactness,
        level_limits: vec![5],
        ..Default::default()
    };
    let mut grid = CombinationSparseGrid::new(
        1,
        vec![GlobalBasis {
            basis_type: GlobalBasisType::ClenshawCurtis,
            custom_rule: None,
        }],
    );
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len());
    for x in grid.nodes.chunks_exact(1) {
        values.push(x[0] * x[0]);
    }
    assert!((1.0 - grid.interpolate(&[0.2], &values, 1)[0] / 0.04).abs() < 1e-12);
}

#[test]
fn test_1d_gauss_legendre() {
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::InterpolationExactness,
        level_limits: vec![10],
        ..Default::default()
    };
    let mut grid = CombinationSparseGrid::new(
        1,
        vec![GlobalBasis {
            basis_type: GlobalBasisType::GaussLegendre,
            custom_rule: None,
        }],
    );
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len());
    for x in grid.nodes.chunks_exact(1) {
        values.push(x[0] * x[0]);
    }
    assert_eq!(grid.len(), 10);
    assert!((1.0 - grid.interpolate(&[0.2], &values, 1)[0] / 0.04).abs() < 1e-12);
}

#[test]
fn test_2d() {
    let ndim = 2;
    let mut grid = CombinationSparseGrid::new(
        ndim,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::ClenshawCurtis,
                custom_rule: None
            };
            ndim
        ],
    );
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::Level,
        level_limits: vec![8; ndim],
        ..Default::default()
    };
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / ndim);
    for x in grid.nodes.chunks_exact(ndim) {
        values.push(x[0] * x[0] + x[1] * x[1]);
    }
    assert!((1.0 - grid.interpolate(&[0.2, 0.2], &values, 1)[0] / 0.08).abs() < 1e-12);
}

#[test]
fn integral_2d()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(2, vec![GlobalBasis{basis_type: GlobalBasisType::GaussPatterson, custom_rule: None}; 2]);
    let bbox = DynamicBoundingBox::new(&[-5.0,-5.0], &[5.0,5.0]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::QuadratureExactness, level_limits: vec![3;2], exactness_limit: Some(12)};
    grid.sparse_grid(options).unwrap();
    *grid.bounding_box_mut() = Some(bbox);
    let mut values = Vec::with_capacity(grid.nodes.len() / 2);
    for x in grid.nodes().chunks_exact(2) {
        values.push(x[0] * x[0] + x[1] * x[1]);
    }
    println!("nodes={}", grid.len());
    println!("integral={}", grid.integral(&values, 1)[0]);
}

#[test]
fn integral_2d_bbox() {
    let mut grid = CombinationSparseGrid::new(
        2,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::GaussPatterson,
                custom_rule: None
            };
            2
        ],
    );
    let bbox = DynamicBoundingBox::new(&[-5.0, -5.0], &[5.0, 5.0]);
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::InterpolationExactness,
        level_limits: vec![3; 2],
        ..Default::default()
    };
    grid.sparse_grid(options).unwrap();
    *grid.bounding_box_mut() = Some(bbox);
    let mut values = Vec::with_capacity(grid.nodes.len() / 2);
    for x in grid.nodes().chunks_exact(2) {
        values.push(x[0] * x[0] + x[1] * x[1]);
    }
    println!("integral={}", grid.integral(&values, 1)[0]);
}

#[test]
fn check_3d_grid() {
    let mut grid = CombinationSparseGrid::new(
        3,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::GaussPatterson,
                custom_rule: None
            };
            3
        ],
    );
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::QuadratureExactness,
        level_limits: vec![2; 3],
        ..Default::default()
    };
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::new();
    for node in grid.nodes.chunks_exact(grid.ndim()) {
        values.push(node[0] * node[0] + node[1] * node[1]);
    }
    println!("{}", grid.integral(&values, 1)[0]);
}

#[test]
fn check_4d_grid() {
    let mut grid = CombinationSparseGrid::new(
        4,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::GaussPatterson,
                custom_rule: None
            };
            4
        ],
    );
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::QuadratureExactness,
        level_limits: vec![2; 4],
        ..Default::default()
    };
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::new();
    for node in grid.nodes.chunks_exact(grid.ndim()) {
        values.push(node[0] * node[0] + node[1] * node[1]);
    }
    println!("{}", grid.integral(&values, 1)[0]);
}

#[test]
fn check_product_grid() {
    let mut grid1 = CombinationSparseGrid::new(
        2,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::GaussPatterson,
                custom_rule: None
            };
            2
        ],
    );
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::Level,
        level_limits: vec![2; 2],
        ..Default::default()
    };
    grid1.sparse_grid(options.clone()).unwrap();
    let mut values = Vec::new();
    for node in grid1.nodes.chunks_exact(grid1.ndim()) {
        values.push(node[0] * node[0] + node[1] * node[1]);
    }

    let mut grid2 = CombinationSparseGrid::new(
        2,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::GaussPatterson,
                custom_rule: None
            };
            2
        ],
    );
    grid2.sparse_grid(options.clone()).unwrap();
    let product_grid = grid1.product(&grid2, options).unwrap();
    let values2 = product_grid.copy_values_to_grid(0, &grid1, &values);

    let integral1 = grid1.integral(&values, 1)[0];
    let integral_product = product_grid.integral(&values2, 1)[0];
    assert!((1.0 - integral1 / integral_product).abs() < 1e-15);
}

#[test]
fn grid_3d_gp() {
    const NDIM: usize = 1;
    let mut grid = CombinationSparseGrid::new(
        NDIM,
        vec![
            GlobalBasis {
                basis_type: GlobalBasisType::GaussPatterson,
                custom_rule: None
            };
            NDIM
        ],
    );
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::QuadratureExactness,
        level_limits: vec![3; NDIM],
        ..Default::default()
    };
    grid.sparse_grid(options).unwrap();
    for node in grid.nodes.chunks_exact(grid.ndim()) {
        println!("{:?}", node);
    }
}

#[test]
fn combination_grid_incremental_generation_matches_rebuild() {
    let basis = vec![
        GlobalBasis {
            basis_type: GlobalBasisType::ClenshawCurtis,
            custom_rule: None
        };
        2
    ];
    let small = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::Level,
        level_limits: vec![2, 2],
        ..Default::default()
    };
    let big = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::Level,
        level_limits: vec![3, 3],
        ..Default::default()
    };

    let mut incremental = CombinationSparseGrid::new(2, basis.clone());
    incremental.sparse_grid(small).unwrap();
    incremental.sparse_grid(big.clone()).unwrap();

    let mut rebuilt = CombinationSparseGrid::new(2, basis);
    rebuilt.sparse_grid(big).unwrap();

    let mut incremental_nodes: Vec<Vec<i64>> =
        incremental.nodes().chunks_exact(2).map(point_key).collect();
    let mut rebuilt_nodes: Vec<Vec<i64>> = rebuilt.nodes().chunks_exact(2).map(point_key).collect();
    incremental_nodes.sort();
    rebuilt_nodes.sort();
    assert_eq!(incremental_nodes, rebuilt_nodes);

    let mut incremental_values = Vec::new();
    for point in incremental.nodes().chunks_exact(2) {
        incremental_values.push(point[0] * point[0] + point[1] * point[1]);
    }

    let mut rebuilt_values = Vec::new();
    for point in rebuilt.nodes().chunks_exact(2) {
        rebuilt_values.push(point[0] * point[0] + point[1] * point[1]);
    }

    assert!(
        (incremental.integral(&incremental_values, 1)[0] - rebuilt.integral(&rebuilt_values, 1)[0])
            .abs()
            < 1e-14
    );
    assert!(
        (incremental.interpolate(&[0.2, 0.2], &incremental_values, 1)[0]
            - rebuilt.interpolate(&[0.2, 0.2], &rebuilt_values, 1)[0])
            .abs()
            < 1e-14
    );
}

#[test]
fn combination_grid_bbox_interpolation_and_fast_integral_regression() {
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::InterpolationExactness,
        level_limits: vec![5],
        ..Default::default()
    };
    let mut grid = CombinationSparseGrid::new(
        1,
        vec![GlobalBasis {
            basis_type: GlobalBasisType::ClenshawCurtis,
            custom_rule: None,
        }],
    );
    grid.sparse_grid(options).unwrap();
    *grid.bounding_box_mut() = Some(DynamicBoundingBox::new(&[-5.0], &[5.0]));

    let mut values = Vec::new();
    for point in grid.nodes().chunks_exact(1) {
        values.push(point[0] * point[0]);
    }

    assert!((grid.interpolate(&[2.0], &values, 1)[0] - 4.0).abs() < 1e-8);
    assert!((grid.integral(&values, 1)[0] - grid.integral_fast(&values, 1)[0]).abs() < 1e-12);
}

#[test]
fn combination_grid_refinement_proposal_apply_matches_refine_with_1d() {
    let basis = vec![GlobalBasis {
        basis_type: GlobalBasisType::ClenshawCurtis,
        custom_rule: None,
    }];
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::InterpolationExactness,
        level_limits: vec![2],
        ..Default::default()
    };
    let refine_options = AdaptiveRefinementOptions {
        threshold: 1e-4,
        level_limits: vec![4],
        max_new_tensors: Some(1),
        indicator_norm: TensorIndicatorNorm::MaxAbs,
    };

    let mut proposal_grid = CombinationSparseGrid::new(1, basis.clone());
    proposal_grid.sparse_grid(options.clone()).unwrap();
    let mut proposal_values = Vec::new();
    for point in proposal_grid.nodes().chunks_exact(1) {
        proposal_values.push((5.0 * point[0]).sin());
    }

    let proposal = proposal_grid
        .propose_refinement(&proposal_values, 1, &refine_options)
        .unwrap();
    assert!(!proposal.candidate_nodes.is_empty());
    let mut candidate_values = Vec::new();
    for point in proposal.candidate_nodes.chunks_exact(1) {
        candidate_values.push((5.0 * point[0]).sin());
    }
    let update = proposal_grid
        .apply_refinement(proposal, &candidate_values)
        .unwrap();
    proposal_values.extend(&update.new_values);

    let mut refine_grid = CombinationSparseGrid::new(1, basis);
    refine_grid.sparse_grid(options).unwrap();
    let mut refine_values = Vec::new();
    for point in refine_grid.nodes().chunks_exact(1) {
        refine_values.push((5.0 * point[0]).sin());
    }
    let accepted = refine_grid
        .refine_with(&mut refine_values, 1, &refine_options, |point| {
            vec![(5.0 * point[0]).sin()]
        })
        .unwrap();

    assert_eq!(accepted, update.accepted_tensors.len());
    assert_eq!(proposal_grid.nodes(), refine_grid.nodes());
    assert_eq!(proposal_values, refine_values);
}

#[test]
fn combination_grid_refinement_proposal_apply_matches_refine_with_2d() {
    let basis = vec![
        GlobalBasis {
            basis_type: GlobalBasisType::GaussPatterson,
            custom_rule: None
        };
        2
    ];
    let options = GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::QuadratureExactness,
        level_limits: vec![1, 1],
        ..Default::default()
    };
    let refine_options = AdaptiveRefinementOptions {
        threshold: 1e-3,
        level_limits: vec![3, 3],
        max_new_tensors: Some(2),
        indicator_norm: TensorIndicatorNorm::MaxAbs,
    };
    let eval =
        |point: &[f64]| vec![((point[0] - 0.3) * 9.0).sin() * ((point[1] - 0.7) * 7.0).cos()];

    let mut proposal_grid = CombinationSparseGrid::new(2, basis.clone());
    proposal_grid.sparse_grid(options.clone()).unwrap();
    let mut proposal_values = Vec::new();
    for point in proposal_grid.nodes().chunks_exact(2) {
        proposal_values.extend(eval(point));
    }

    let proposal = proposal_grid
        .propose_refinement(&proposal_values, 1, &refine_options)
        .unwrap();
    assert!(!proposal.candidate_nodes.is_empty());
    let mut candidate_values = Vec::new();
    for point in proposal.candidate_nodes.chunks_exact(2) {
        candidate_values.extend(eval(point));
    }
    let update = proposal_grid
        .apply_refinement(proposal, &candidate_values)
        .unwrap();
    proposal_values.extend(&update.new_values);

    let mut refine_grid = CombinationSparseGrid::new(2, basis);
    refine_grid.sparse_grid(options).unwrap();
    let mut refine_values = Vec::new();
    for point in refine_grid.nodes().chunks_exact(2) {
        refine_values.extend(eval(point));
    }
    let accepted = refine_grid
        .refine_with(&mut refine_values, 1, &refine_options, eval)
        .unwrap();

    assert_eq!(accepted, update.accepted_tensors.len());
    assert_eq!(proposal_grid.nodes(), refine_grid.nodes());
    assert_eq!(proposal_values, refine_values);
}
