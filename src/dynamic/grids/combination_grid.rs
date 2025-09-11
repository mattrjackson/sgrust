
use crate::{basis::base::Basis, dynamic::{algorithms::lagrange::{lagrange_coeffs, lagrange_weights}, iterators::advanced_step_iterator::{AdvancedStepIterator, GridType}}, errors::SGError};
use crate::basis::global::GlobalBasis;
use kdtree::KdTree;
use kdtree::distance::squared_euclidean;
use ndarray::{ArrayView1, ArrayView2};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};


fn cartesian_product(bounds: &[u32]) -> Vec<u32> 
{
    let ndim = bounds.len();
    // Calculate the total number of combinations by multiplying all bounds
    let total_combinations = bounds.iter().product::<u32>() as usize;

    // Allocate a flat vector to hold all indices in a flattened form
    // (each multi-index will be a consecutive slice of length `bounds.len()` in this vector)
    let mut multi_indices = vec![0; total_combinations * bounds.len()];
    // Iterate through all possible combinations
    for (i, current) in multi_indices.chunks_exact_mut(ndim).enumerate() {
        let mut index = i as u32;
        for (j, &bound) in bounds.iter().rev().enumerate() {
            current[ndim - 1 - j] = index % (bound);
            index /= bound ;
        }
    }
    multi_indices
}

#[derive(Default, Serialize, Deserialize, Clone, Copy)]
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
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
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
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
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
pub struct DynamicBoundingBox
{
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
}
impl DynamicBoundingBox
{
    pub fn new(lower: &[f64], upper: &[f64]) -> Self
    {
        Self { lower: lower.to_owned(), upper: upper.to_owned() }
    }

    #[inline]
    pub fn ndim(&self) -> usize
    {
        self.lower.len()
    }

    #[inline]
    pub fn width(&self, dim: usize) -> f64
    {
        self.upper[dim] - self.lower[dim]
    }

    ///
    /// Volume of hypercube (width(dim1)*...*width(dim_n))
    /// 
    pub fn volume(&self) -> f64
    {
        let mut volume = 1.0;
        for d in 0..self.ndim()
        {
            volume *= self.width(d);
        }
        volume
    }
    pub fn to_unit_coordinate(&self, point: &[f64]) -> Vec<f64>
    {
        let mut r = vec![0.0; self.ndim()];
        for i in 0..self.ndim()
        {
            r[i] = (point[i] - self.lower[i])/(self.upper[i] - self.lower[i]);
        }
        r
    }
    pub fn to_real_coordinate(&self, point: &mut [f64])
    {        
        #[allow(clippy::needless_range_loop)]
        for i in 0..self.ndim()
        {
            point[i] = self.lower[i] + (self.upper[i] - self.lower[i]) * point[i];
        }
    }
    pub fn real_coordinate(&self, x: f64, dim: usize) -> f64
    {
        self.lower[dim] + (self.upper[dim] - self.lower[dim]) *x
    }   
    pub fn contains(&self, point: &[f64]) -> bool
    {
        #[allow(clippy::needless_range_loop)]
        for d in 0..self.ndim()
        {
            if self.lower[d] > point[d] || self.upper[d] < point[d]
            {
                return false;
            }
        }
        true
    }
}

#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
pub struct NodesAndCoefficientsForLevel
{
    pub level: u32,
    pub nodes: Vec<f64>,
    pub coefficients: Vec<f64>
}

#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
pub struct FullGrid
{
    pub level: Vec<u32>,    
    pub lweights: Vec<NodesAndCoefficientsForLevel>,
    pub weight: f64,
    index_map: Vec<u32>,
    combinations: Vec<u32>,
}

impl FullGrid
{
    ///
    /// Number of dimensions
    /// 
    pub fn ndim(&self) ->usize
    {
        self.level.len()
    }

    ///
    /// Returns the point indices for the `basis` functions. 
    /// 
    pub fn points(&self, basis: &Vec<GlobalBasis>) -> Vec<u32>
    {
        let num_points: Vec<u32> = self.level.iter().zip(basis).map(|(&level, basis)| basis.num_nodes(level) as u32).collect();
        cartesian_product(&num_points)
    }

    ///
    /// Create a new `FullGrid` for a given `level`, using the specified `basis` functions, 
    /// and the global `weight`for the grid.
    /// 
    pub fn new(basis: &[GlobalBasis], level: &[u32], weight: f64) -> Self
    {
        // Determine the number of points in each dimension from the basis functions
        let npts: Vec<u32> = level.iter().zip(basis).map(|(&level, basis)| 
          basis.num_nodes(level) as u32).collect();
        // we store combinations in order to speed up interpolation later...
        let combinations = Self::generate_combinations(&npts);
        let mut r = Self{level: level.to_owned(), lweights: Vec::new(), weight, index_map: Vec::new(), combinations };
        
        // Compute lagrange coefficients for tensor
        let mut coeffs = Vec::new();
        for d in 0..level.len()
        {
            let nodes = basis[d].nodes(level[d]);
            coeffs.push(NodesAndCoefficientsForLevel{ coefficients: lagrange_coeffs(&nodes), level: level[d], nodes });
        }
        r.lweights = coeffs;
        r
    }
    #[inline]
    fn langrange_weight_size(&self) -> usize
    {
        let mut r = 1;
        for item in &self.lweights
        {
            r *= item.nodes.len();
        }
        r
    }
    fn generate_combinations(npts: &[u32], ) -> Vec<u32>
    {
        let mut combinations = Vec::new();
        let total_combinations = npts.iter().product::<u32>() as usize;
        let ndim = npts.len();
        let mut current = vec![0; ndim];        
        // Iterate through all possible combinations
        (0..total_combinations).for_each(|i|
        {
            let mut index = i as u32;
            for (j, &bound) in npts.iter().rev().enumerate() {
                current[ndim - 1 - j] = index % (bound);
                index /= bound ;
            }
            combinations.extend(&current);

            current.fill(0);
        });
        combinations
    }

    ///
    /// Computes the lagrange weights for a given `x` coordinate.
    /// Returns vector starting indices for each dimension and a vector 
    /// containing weights across all dimensions. 
    /// 
    fn lagrange_weights(&self, x: &[f64]) -> (Vec<usize>, Vec<f64>)
    {
        let mut weights: Vec<f64> = Vec::with_capacity(self.langrange_weight_size());
        let mut indices: Vec<usize> = vec![0; self.ndim()];
        for dim in 0..self.ndim()
        {                        
            indices[dim] = weights.len();
            weights.extend(lagrange_weights(x[dim], &self.lweights[dim].coefficients,
                &self.lweights[dim].nodes));
        }
        (indices, weights)
    }
    ///
    /// Interpolate at `x` along a given grid using the specified `basis` functions. 
    /// 
    #[inline]
    pub fn interpolate(&self, x: &[f64], y: &mut [f64], values: &[f64])
    {
        // Retrieve the lagrange weights for this coordinate
        let (indices, weights) = self.lagrange_weights(x);
        // Iterate through all possible combinations
        for (i, combination) in self.combinations.chunks_exact(self.ndim()).enumerate()
        {
            let mut combined_weight = 1.0;            
            (0..self.ndim()).for_each(|dim|
            {
                let index = indices[dim];
                combined_weight *= weights[index + combination[dim] as usize];
            });
            // Retrive the value(s) for this node by using the `index_map` to determine 
            // the global index from the local grid index
            let values = &values[y.len()*(self.index_map[i] as usize)..];
            // Add contribution in each dimension.
            (0..y.len()).for_each(|i|
            {
                y[i] += combined_weight* values[i] * self.weight;
            });
        }
    }
}

///
/// Generates a sparse grid using the combination technique. Unlike the 
/// `LinearSparseGrid`, this grid is not presently adaptively refineable. 
/// Instead, it is primarily optimized for uncertainty quantification 
/// operations (roughly equivalent to the `GlobalGrid` concept found in
/// TASMANIAN).
/// 
#[derive(Clone, Default, Serialize, Deserialize)]
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
pub struct CombinationSparseGrid
{
    bounding_box: Option<DynamicBoundingBox>,
    ndim: usize,
    basis: Vec<GlobalBasis>,    
    grids: Vec<FullGrid>,    
    nodes: Vec<f64>,
    /// a hash of grid levels we already have...
    grid_hash: FxHashSet<Vec<u32>>,
    /// quadrature weight for node...
    qweight: Vec<f64>,
}
impl CombinationSparseGrid
{
    ///
    /// Create an empty Grid with dimension `ndim` and
    /// basis functions `basis`.
    /// 
    pub fn new(ndim: usize, basis: Vec<GlobalBasis>) -> Self
    {
        assert!(basis.len() == ndim);
        Self { ndim, basis, ..Default::default() }
    }

    pub fn len(&self) -> usize
    {
        self.nodes.len() / self.ndim
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    ///
    /// Number of dimensions
    /// 
    pub fn ndim(&self) -> usize
    {
        self.ndim
    }

    ///
    /// Compute the quadrature weight for a given `level` and `index` using 
    /// the user-specified `basis` functions.
    /// 
    pub(crate) fn quadrature_weight(basis: &[GlobalBasis], index: &[u32], level: &[u32]) -> f64
    {
        let mut weight = 1.0;
        for d in 0..basis.len()
        {
            weight *= basis[d].integral(level[d], index[d]);
        }
        weight
    }
    
    ///
    /// Compute the floating point coordinates for a given `level` and `index`
    /// 
    fn unit_coordinate(basis: &[GlobalBasis], level: &[u32], index: &[u32]) -> Vec<f64>
    {
        let ndim = basis.len();
        let mut coor = vec![0.0; ndim];
        for d in 0..ndim
        {
            coor[d] = basis[d].get_point(level[d], index[d]);
        }
        coor
    }    

    ///
    /// Builds a kd-tree from the existing nodes. Used internally
    /// to generate new grid points.
    /// 
    fn build_kdtree(nodes: &[f64], ndim: usize) -> Result<KdTree<f64, usize, Vec<f64>>, kdtree::ErrorKind>
    {
        let mut tree= KdTree::new(ndim);
        for (i, node) in nodes.chunks_exact(ndim).enumerate()
        {
            tree.add(node.to_vec(), i)?;
        }
        Ok(tree)
    }
    ///
    /// Generate a standard sparse grid using the combination technique.
    /// 
    pub fn sparse_grid(&mut self, generation_parameters: GenerationOptions) -> Result<(), SGError>
    {
        self.generate_grid(GridType::Sparse, generation_parameters)
    }

    ///
    /// Compute a full grid using the combination technique
    /// 
    pub fn full_grid(&mut self, generation_parameters: GenerationOptions) -> Result<(), SGError>
    {
        self.generate_grid( GridType::Full, generation_parameters)
    }
    fn generate_grid(&mut self, grid_type: GridType, generation_parameters: GenerationOptions) -> Result<(), SGError>
    {   
        let mut tree= Self::build_kdtree(&self.nodes, self.ndim).map_err(|_|SGError::KdTreeError)?;        
        let exactness_bound = match generation_parameters.tensor_strategy
        {
            TensorSelectionStrategy::Level => u32::MAX,
            TensorSelectionStrategy::QuadratureExactness => 
            {
                let max_exactness = self.basis.iter().zip(&generation_parameters.level_limits).                
                max_by(|a,b|a.0.quadrature_exactness(a.1 - 1).cmp(&b.0.quadrature_exactness(b.1 - 1))).unwrap();
                max_exactness.0.quadrature_exactness(max_exactness.1 - 1) + 1
            },
            TensorSelectionStrategy::InterpolationExactness =>
            {
                let max_exactness = self.basis.iter().zip(&generation_parameters.level_limits).
                max_by(|a,b|a.0.interpolation_exactness(a.1 - 1).cmp(&b.0.interpolation_exactness(b.1 - 1))).unwrap();
                max_exactness.0.interpolation_exactness(max_exactness.1 - 1) + 1
            }
        };
        let indices: Vec<_> = if self.ndim == 1 { generation_parameters.level_limits.clone() } else { AdvancedStepIterator::new(&generation_parameters.level_limits, grid_type, 
            self.basis.clone(), generation_parameters.tensor_strategy,  generation_parameters.exactness_limit.unwrap_or(exactness_bound)).flatten().collect() };
        let weights = crate::utilities::multi_index_manipulation::weight_modifiers(&indices, self.ndim)?;        
        for (level_set, weight) in indices.chunks_exact(self.ndim).zip(weights)
        {
            // Skip if weight is zero or we have already added this set
            if weight == 0.0 || self.grid_hash.contains(level_set)
            {
                continue;
            }
            let mut grid = FullGrid::new(&self.basis, level_set, weight);
                       
            // loop through indices of level_set...
            for index in grid.points(&self.basis).chunks_exact(self.ndim)
            {
                let point = Self::unit_coordinate(&self.basis, &grid.level, index);
                // add new nodes
                let nodes = tree.within(&point, 1e-15, &squared_euclidean).map_err(|_|SGError::KdTreeError)?;  
                if nodes.is_empty()
                {
                    let idx = self.nodes.len() / self.ndim;
                    tree.add(point.clone(), idx).unwrap();
                    self.nodes.extend(&point);
                    self.qweight.push(weight*Self::quadrature_weight(&self.basis, index, level_set));                    
                    grid.index_map.push(idx as u32);
                }                
                else
                {
                    self.qweight[*nodes[0].1] += weight*Self::quadrature_weight(&self.basis, index, level_set);
                    grid.index_map.push(*nodes[0].1 as u32);
                }
            }
            self.grids.push(grid);            
            self.grid_hash.insert(level_set.to_owned());
        }
        Ok(())
    }    

    ///
    /// For internal use. Generates a HashMap to map the floating point values stored
    /// to their index. 
    fn node_map(&self) -> FxHashMap<Vec<u64>, usize>
    {
        let mut map = FxHashMap::default();
        for (idx, node) in self.nodes.chunks_exact(self.ndim).enumerate()
        {
            let index: Vec<u64> = node.iter().map(|&x|x.to_bits()).collect();
            map.insert(index.clone(), idx);
        }
        map
    }

    ///
    /// Return reference to bounding box
    /// 
    pub fn bounding_box(&self) -> &Option<DynamicBoundingBox>
    {
        &self.bounding_box
    }

    ///
    /// Return mutable reference for bounding box
    /// 
    pub fn bounding_box_mut(&mut self) -> &mut Option<DynamicBoundingBox>
    {
        &mut self.bounding_box
    }
    ///
    /// Return vector containing real coordinates for grid (size = `ndim` * nodes().len()).
    /// 
    pub fn nodes(&self) -> Vec<f64>
    {
        
        let mut r = vec![0.0; self.nodes.len()];
        if let Some(bbox) = &self.bounding_box
        {
            for (node, point) in self.nodes.chunks_exact(self.ndim).zip(r.chunks_exact_mut(self.ndim))
            {
                point.copy_from_slice(node);
                bbox.to_real_coordinate(point);
            }
        }
        else
        {
            for (node, point) in self.nodes.chunks_exact(self.ndim).zip(r.chunks_exact_mut(self.ndim))
            {            
                point.copy_from_slice(node);
            }
        }
        r
    }

    ///
    /// Compute integral over grid
    /// 
    pub fn integral(&self, values: &[f64], num_outputs: usize) -> Vec<f64>
    {
        let mut y = vec![0.0; num_outputs];
        for (&weight,value) in self.qweight.iter().zip(values.chunks_exact(num_outputs))
        {
            for i in 0..num_outputs
            {
                y[i] += weight*value[i];
            }
        }
        if let Some(bbox) = &self.bounding_box
        {
            let volume = bbox.volume();
            #[allow(clippy::needless_range_loop)]
            for i in 0..num_outputs
            {
                y[i] *= volume;
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
        let values = ArrayView2::from_shape((self.len(), num_outputs), &values).unwrap();
           
        let product = weights.dot(&values);        
        let mut y = Vec::from(product.as_slice().unwrap());       
        if let Some(bbox) = &self.bounding_box
        {
            let volume = bbox.volume();
            #[allow(clippy::needless_range_loop)]
            for i in 0..num_outputs
            {
                y[i] *= volume;
            }
        }
        y
    } 

    ///
    /// Interpolate on grid using previously computed `values`
    /// 
    pub fn interpolate(&self, x: &[f64], values: &[f64], num_outputs: usize) -> Vec<f64>
    {
        let x = x.to_owned();
    
        if let Some(bbox) = &self.bounding_box
        {
            bbox.to_unit_coordinate(&x);
        }
        let mut y = vec![0.0; num_outputs];
        if self.ndim == 1
        {
            let grid = self.grids.last().unwrap();
            let lweights = grid.lweights.last().unwrap();
            let weights = lagrange_weights(x[0],&lweights.coefficients, &self.basis[0].nodes(grid.level[0]));
            
            for (&weight, value) in weights.iter().zip(values.chunks_exact(num_outputs))
            {                   
                for i in 0..num_outputs
                {
                    y[i] += weight * value[i];
                }
            }
        }
        else 
        {
            for grid in self.grids.iter()
            {
                grid.interpolate(&x, &mut y, values);                
            }
        }
        y
    }

    ///
    /// Compute the product of two sparse grids.
    /// 
    pub fn product(&self, other: &CombinationSparseGrid, mut options: GenerationOptions) -> Result<CombinationSparseGrid, SGError>
    {
        let mut basis_vectors = Vec::new();
        basis_vectors.extend(self.basis.clone());
        basis_vectors.extend(other.basis.clone());
        let mut lower = Vec::new();
        let mut upper = Vec::new();
        if let Some(bbox) = &self.bounding_box
        {
            lower.extend(bbox.lower.clone());
            upper.extend(bbox.upper.clone());
            if let Some(bbox2) = &other.bounding_box
            {
                lower.extend(bbox2.lower.clone());
                upper.extend(bbox2.upper.clone());
            }
            else 
            {
                lower.extend(vec![0.0; other.basis.len()]);
                upper.extend(vec![1.0; other.basis.len()]);
            }
        }
        else
        {
            // Only need to create bbox if the second grid has a bounding box defined.
            if let Some(bbox2) = &other.bounding_box
            {
                lower.extend(vec![0.0; self.basis.len()]);                
                lower.extend(bbox2.lower.clone());

                upper.extend(vec![1.0; self.basis.len()]);
                upper.extend(bbox2.upper.clone());
            }
        }
        let mut grid = CombinationSparseGrid::new(self.ndim + other.ndim, basis_vectors.clone());
        grid.bounding_box = if lower.is_empty() { None } else { Some(DynamicBoundingBox::new(&lower, &upper)) };
        let mut levels = vec![0; basis_vectors.len()];
        for level in &self.grid_hash
        {
            for i in 0..level.len()
            {
                levels[i] = levels[i].max(level[i]);
            }
        }
        let offset = self.basis.len();
        for level in &other.grid_hash
        {
            for i in 0..level.len()
            {
                levels[i + offset] = levels[i + offset].max(level[i]);
            }
        }
        options.level_limits.clear();
        options.level_limits.extend(self.max_limits());
        options.level_limits.extend(other.max_limits());
        grid.sparse_grid(options)?;
        Ok(grid)
    }

    ///
    /// Copy values existing on previous grid into new one. 
    /// 
    pub fn copy_values_to_grid<T: Clone+Default>(&self, offset: usize, source: &CombinationSparseGrid, source_values: &[T]) -> Vec<T>
    {
        let mut v = vec![T::default(); self.len()];
        let source_node_map = source.node_map();
        let map = self.node_map();       
        for (node, &index) in map.iter()
        {
            let source_idx = &node[offset..offset + source.basis.len()];

            if let Some(&idx) = source_node_map.get(source_idx)
            {
                v[index] = source_values[idx].clone();
            }
        }
        v
    }
    /// Get basis functions associated with this grid
    pub fn basis(&self) -> &Vec<GlobalBasis>
    {
        &self.basis
    }
    /// Get level limits for each dimension 
    pub fn max_limits(&self) -> Vec<u32>
    {
        let mut limits = vec![0; self.ndim];
        for grid in &self.grids
        {
            for (i, &level) in grid.level.iter().enumerate()
            {
                limits[i] = limits[i].max(level);
            }
        }
        limits
    }
}

#[test]
fn test_combination_grid()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(2, vec![GlobalBasis{basis_type: GlobalBasisType::ClenshawCurtis, custom_rule: None}; 2]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::Level, level_limits: vec![3,3], ..Default::default()};
    grid.full_grid( options.clone()).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / 2 );
    for x in grid.nodes.chunks_exact(2)
    {                
        values.push(x[0] * x[0] + x[1]*x[1]);
    } 
    assert!(grid.len() == options.level_limits.iter().map(|&l| (1 << l) + 1).product::<usize>());
    println!("integral={}", grid.integral(&values, 1)[0]);
}


#[test]
fn test_combination_grid_gauss_legendre()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(2, vec![GlobalBasis{basis_type: GlobalBasisType::GaussLegendre, custom_rule: None}; 2]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::Level, level_limits: vec![10,10], ..Default::default()};
    grid.full_grid( options.clone()).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / 2 );
    for x in grid.nodes.chunks_exact(2)
    {                
        values.push(x[0] * x[0] + x[1]*x[1]);
    } 
    assert!(grid.len() == 100);
    println!("integral={}", grid.integral(&values, 1)[0]);
}
#[test]
fn test_combination_grid_integral_standard()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(2, vec![GlobalBasis{basis_type: GlobalBasisType::ClenshawCurtis, custom_rule: None}; 2]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::Level, level_limits: vec![3,3], ..Default::default()};
    grid.full_grid( options).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / 2 );
    for x in grid.nodes.chunks_exact(2)
    {                
        values.push(x[0] * x[0] + x[1]*x[1]);
    } 
    let t1 = std::time::Instant::now();
    for _ in 0..1e6 as usize
    {
        let _ = grid.integral(&values, 1)[0];
    }
    println!("1e5 interpolation took {} msec", std::time::Instant::now().duration_since(t1).as_millis());
    println!("integral={}", grid.integral(&values, 1)[0]);
}


#[test]
fn test_combination_grid_integral_ndarray()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(2, vec![GlobalBasis{basis_type: GlobalBasisType::ClenshawCurtis, custom_rule: None}; 2]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::Level, level_limits: vec![3,3], ..Default::default()};
    grid.full_grid( options).unwrap();
    let mut values = vec![0.0; grid.len()];
    for (i,x) in grid.nodes.chunks_exact(2).enumerate()
    {   
         values[i] = x[0] * x[0] + x[1]*x[1];
        
    } 
    let t1 = std::time::Instant::now();
    for _ in 0..1e6 as usize
    {
        let _ = grid.integral_fast(&values, 1)[0];
    }
    println!("1e5 interpolation took {} msec", std::time::Instant::now().duration_since(t1).as_millis());
    println!("integral={}", grid.integral(&values, 1)[0]);
}

#[test]
fn test_1d()
{
    use crate::basis::global::GlobalBasisType;
    let ndim = 1;
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::InterpolationExactness, level_limits: vec![5], ..Default::default()};
    let mut grid = CombinationSparseGrid::new(1, vec![GlobalBasis{basis_type: GlobalBasisType::ClenshawCurtis, custom_rule: None}; 1]);
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / ndim );
    for x in grid.nodes.chunks_exact(ndim)
    {                
        values.push(x[0] * x[0]);
    } 
    println!("integral={}", grid.integral(&values, 1)[0]);
    let t1 = std::time::Instant::now();
    for _ in 0..1e6 as usize
    {
        let _ = grid.interpolate(&[0.2], &values, 1)[0];
    }
    println!("1e5 interpolation took {} msec", std::time::Instant::now().duration_since(t1).as_millis());
    println!("interpolated value @0.2={}", grid.interpolate(&[0.2], &values, 1)[0]);
    assert!((1.0-grid.interpolate(&[0.2], &values, 1)[0] / (0.04)).abs() < 1e-12);
}

#[test]
fn test_1d_gauss_legendre()
{
    use crate::basis::global::GlobalBasisType;
    let ndim = 1;
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::InterpolationExactness, level_limits: vec![10], ..Default::default()};
    let mut grid = CombinationSparseGrid::new(1, vec![GlobalBasis{basis_type: GlobalBasisType::GaussLegendre, custom_rule: None}; 1]);
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / ndim );
    for x in grid.nodes.chunks_exact(ndim)
    {                
        values.push(x[0] * x[0]);
    } 
    println!("integral={}", grid.integral(&values, 1)[0]);
    let t1 = std::time::Instant::now();
    for _ in 0..1e6 as usize
    {
        let _ = grid.interpolate(&[0.2], &values, 1)[0];
    }
    assert!(grid.len() == 10);
    println!("1e5 interpolation took {} msec", std::time::Instant::now().duration_since(t1).as_millis());
    println!("interpolated value @0.2={}", grid.interpolate(&[0.2], &values, 1)[0]);
    assert!((1.0-grid.interpolate(&[0.2], &values, 1)[0] / (0.04)).abs() < 1e-12);
}

#[test]
fn test_2d()
{
    use crate::basis::global::GlobalBasisType;
    let ndim = 2;
    let mut grid = CombinationSparseGrid::new(ndim, vec![
        GlobalBasis{basis_type: GlobalBasisType::ClenshawCurtis, custom_rule: None}; ndim]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::Level, level_limits: vec![8; ndim], ..Default::default()}; 
    grid.sparse_grid( options).unwrap();
    let mut values = Vec::with_capacity(grid.nodes.len() / ndim );
    for x in grid.nodes.chunks_exact(ndim)
    {                
        values.push(x[0] * x[0] + x[1]*x[1]);
    } 
    println!("integral={}", grid.integral(&values, 1)[0]);
    let t1 = std::time::Instant::now();
    for _ in 0..1e5 as usize
    {
        let _ = grid.interpolate(&[0.2, 0.2], &values, 1)[0];
    }
    println!("1e5 interpolation took {} msec", std::time::Instant::now().duration_since(t1).as_millis());
    println!("{}",grid.interpolate(&[0.2, 0.2], &values, 1)[0] );
    assert!((1.0-grid.interpolate(&[0.2, 0.2], &values, 1)[0] / (0.08)).abs() < 1e-12);
}
#[test]
fn integral_2d()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(2, vec![GlobalBasis{basis_type: GlobalBasisType::GaussPatterson, custom_rule: None}; 2]);
    let bbox = DynamicBoundingBox::new(&[-5.0,-5.0], &[5.0,5.0]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::QuadratureExactness, level_limits: vec![3;2], exactness_limit: Some(12), ..Default::default()};
    grid.sparse_grid(options).unwrap();
    *grid.bounding_box_mut() = Some(bbox);
    let mut values = Vec::with_capacity(grid.nodes.len() / 2 );
    for x in grid.nodes().chunks_exact(2)
    {                
        values.push(x[0] * x[0] + x[1]*x[1]);
    } 
    println!("nodes={}", grid.len());
    println!("integral={}", grid.integral(&values, 1)[0]);
}
#[test]
fn integral_2d_bbox()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(2, vec![GlobalBasis{basis_type: GlobalBasisType::GaussPatterson, custom_rule: None}; 2]);
    let bbox = DynamicBoundingBox::new(&[-5.0,-5.0], &[5.0,5.0]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::InterpolationExactness, level_limits: vec![3;2], ..Default::default()};
    grid.sparse_grid(options).unwrap();
    *grid.bounding_box_mut() = Some(bbox);
    let mut values = Vec::with_capacity(grid.nodes.len() / 2 );
    for x in grid.nodes().chunks_exact(2)
    {                
        values.push(x[0] * x[0] + x[1]*x[1]);
    } 
    println!("integral={}", grid.integral(&values, 1)[0]);
}

#[test]
fn check_3d_grid()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(3, vec![GlobalBasis{basis_type: GlobalBasisType::GaussPatterson, custom_rule: None}; 3]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::QuadratureExactness, level_limits: vec![2; 3], ..Default::default()};
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::new();
    for node in grid.nodes.chunks_exact(grid.ndim())
    {
        values.push(node[0]*node[0] + node[1]*node[1]);
    }
    let integral = grid.integral(&values, 1)[0];    
    println!("{integral}");

}

#[test]
fn check_4d_grid()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(4, vec![GlobalBasis{basis_type: GlobalBasisType::GaussPatterson, custom_rule: None}; 4]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::QuadratureExactness, level_limits: vec![2; 4], ..Default::default()};
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::new();
    for node in grid.nodes.chunks_exact(grid.ndim())
    {
        values.push(node[0]*node[0] + node[1]*node[1]);
    }
    let integral = grid.integral(&values, 1)[0];    
    println!("{integral}");

}

/// Checks the product of two grids.
#[test]
fn check_product_grid()
{
    use crate::basis::global::GlobalBasisType;
    let mut grid1 = CombinationSparseGrid::new(2, vec![GlobalBasis{basis_type: GlobalBasisType::GaussPatterson, custom_rule: None}; 2]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::Level, level_limits: vec![2; 2], ..Default::default()};
    grid1.sparse_grid(options.clone()).unwrap();
    let mut values = Vec::new();
    for node in grid1.nodes.chunks_exact(grid1.ndim())
    {
        values.push(node[0]*node[0] + node[1]*node[1]);
    }
    let mut grid2 = CombinationSparseGrid::new(2, vec![GlobalBasis{basis_type: GlobalBasisType::GaussPatterson, custom_rule: None}; 2]);
    grid2.sparse_grid(options.clone()).unwrap();
    let product_grid = grid1.product(&grid2, options).unwrap();
    let values2 = product_grid.copy_values_to_grid(0, &grid1, &values);
    
    let integral1 = grid1.integral(&values, 1)[0];
    let integral_product = product_grid.integral(&values2, 1)[0];
    assert!((1.0-integral1/integral_product).abs() < 1e-15);
}

#[test]
fn grid_3d_gp()
{
    const NDIM: usize = 1;
    use crate::basis::global::GlobalBasisType;
    let mut grid = CombinationSparseGrid::new(NDIM, vec![GlobalBasis{basis_type: GlobalBasisType::GaussPatterson, custom_rule: None}; NDIM]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::QuadratureExactness, level_limits: vec![3; NDIM], ..Default::default()};
    grid.sparse_grid(options).unwrap();

    for node in grid.nodes.chunks_exact(grid.ndim())
    {
        println!("{:?}", node);
    }

}