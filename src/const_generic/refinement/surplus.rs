use crate::const_generic::{algorithms::refinement::RefinementFunctor, storage::PointIterator};

pub struct SurplusRefinement<const D: usize, const DIM_OUT: usize>;

impl<const D: usize, const DIM_OUT: usize> RefinementFunctor<D, DIM_OUT> for SurplusRefinement<D, DIM_OUT>
{
    fn eval(&self, _points: PointIterator<D>, alpha: &[[f64; DIM_OUT]], _values: &[[f64; DIM_OUT]]) -> Vec<f64>
    {
        alpha.iter().map(|alpha_i|
        {
            let mut max = -1.0_f64;
            alpha_i.iter().for_each(|&val| max = max.max(val.abs()));
            max
        }).collect()
    }

    fn eval_per_dimension(&self, points: PointIterator<D>, alpha: &[[f64; DIM_OUT]], _values: &[[f64; DIM_OUT]]) -> Vec<Vec<f64>>
    {
        // Access the underlying grid points through the iterator's data field
        let grid_points = points.data;

        // For each point, estimate per-dimension error contribution
        grid_points.iter().enumerate().map(|(idx, point)| {
            let alpha_i = &alpha[idx];
            let mut max_surplus = -1.0_f64;
            alpha_i.iter().for_each(|&val| max_surplus = max_surplus.max(val.abs()));

            // Estimate per-dimension contribution based on hierarchical level
            // Higher level in a dimension indicates more variation in that dimension
            // Use a heuristic: error_dim[d] = surplus * (1 + level[d] / max_level)
            let max_level = point.level_max() as f64;
            let min_level = point.level_min() as f64;
            let level_range = (max_level - min_level).max(1.0);

            let dim_errors: Vec<f64> = (0..D).map(|d| {
                // Higher level in dimension d means more refinement needed there
                let level_factor = 1.0 + (point.level[d] as f64 - min_level) / level_range;
                max_surplus * level_factor
            }).collect();

            dim_errors
        }).collect()
    }
}