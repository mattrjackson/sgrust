use crate::dynamic::{algorithms::refinement::RefinementFunctor, storage::PointIterator};

pub struct SurplusRefinement(pub usize, pub usize);

impl RefinementFunctor for SurplusRefinement
{
    fn eval(&self, _points: PointIterator, alpha: &[f64], _values: &[f64]) -> Vec<f64>
    {
        alpha.chunks_exact(self.num_outputs()).map(|alpha_i|
        {
            let mut max = -1.0_f64;
            alpha_i.iter().for_each(|&val| max = max.max(val.abs()));
            max
        }).collect()
    }

    fn eval_per_dimension(&self, points: PointIterator, alpha: &[f64], _values: &[f64]) -> Vec<Vec<f64>>
    {
        let num_dims = self.num_inputs();
        let num_outputs = self.num_outputs();

        // Access the underlying storage through the iterator's internal state
        // We'll use the grid point structure from the storage
        let storage = points.storage;

        // For each point, estimate per-dimension error contribution
        storage.nodes().enumerate().map(|(idx, point)| {
            let alpha_i = &alpha[idx * num_outputs..(idx + 1) * num_outputs];
            let mut max_surplus = -1.0_f64;
            alpha_i.iter().for_each(|&val| max_surplus = max_surplus.max(val.abs()));

            // Estimate per-dimension contribution based on hierarchical level
            // Higher level in a dimension indicates more variation in that dimension
            // Use a heuristic: error_dim[d] = surplus * (1 + level[d] / max_level)
            let max_level = point.level_max() as f64;
            let min_level = point.level_min() as f64;
            let level_range = (max_level - min_level).max(1.0);

            let dim_errors: Vec<f64> = (0..num_dims).map(|d| {
                // Higher level in dimension d means more refinement needed there
                let level_factor = 1.0 + (point.level[d] as f64 - min_level) / level_range;
                max_surplus * level_factor
            }).collect();

            dim_errors
        }).collect()
    }

    #[inline]
    fn num_outputs(&self) -> usize {
        self.1
    }

    #[inline]
    fn num_inputs(&self) -> usize {
        self.0
    }
}