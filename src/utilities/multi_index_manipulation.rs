use std::cmp::Ordering;

use crate::errors::SGError;

///
/// Sot level sets lexographically.. Used by `weight_modifiers`.
/// 
fn sort_level_sets(level_sets: &[u32], ndim: usize) -> Result<(Vec<Vec<u32>>, Vec<Vec<u32>>), SGError> {
    let num_tensors = level_sets.len() / ndim;
    let num_levels = level_sets
        .chunks(ndim)
        .map(|idx| *idx.iter().max().unwrap_or(&0))
        .max()
        .unwrap_or(0)
        + 1;

    // Initialize maps and lines
    let mut sorted_maps: Vec<Vec<u32>> = vec![(0..num_tensors as u32).collect(); ndim];
    let mut lines_1d: Vec<Vec<u32>> = vec![Vec::with_capacity(num_levels as usize); ndim];

    // Helper closure to check equality outside a given dimension
    let match_outside_dim = |dim: usize, a: &[u32], b: &[u32]| -> bool {
        a.iter()
            .zip(b.iter())
            .enumerate()
            .all(|(j, (&x, &y))| j == dim || x == y)
    };

    for dim in 0..ndim {
        // Sort tensor indices lexicographically, ignoring the active `dim` dimension
        sorted_maps[dim].sort_by(|&a, &b| {
            let idx_a = level_sets.chunks(ndim).nth(a as usize).unwrap_or_default();
            let idx_b = level_sets.chunks(ndim).nth(b as usize).unwrap_or_default();

            idx_a.iter()
                .zip(idx_b.iter())
                .enumerate()
                .filter(|(j, _)| *j != dim) // Ignore the current dimension
                .find_map(|(_, (&v_a, &v_b))| match v_a.cmp(&v_b) {
                    Ordering::Equal => None,
                    other => Some(other),
                })
                .unwrap_or(Ordering::Equal)
        });

        // Identify boundary positions (segments) where levels differ outside `dim`
        let mut current_idx = level_sets.chunks(ndim).nth(sorted_maps[dim][0] as usize).ok_or_else(||SGError::InvalidIndex)?;
        lines_1d[dim].push(0);

        for (i, &tensor) in sorted_maps[dim].iter().enumerate().skip(1) {
            let next_idx = level_sets.chunks(ndim).nth(tensor as usize).ok_or_else(||SGError::InvalidIndex)?;
            if !match_outside_dim(dim, current_idx, next_idx) {
                lines_1d[dim].push(i as u32);
                current_idx = next_idx;
            }
        }

        lines_1d[dim].push(num_tensors as u32);
    }

    Ok((sorted_maps, lines_1d))
}

///
/// Compute tensor weight modifiers. Approach adapted from TASMANIAN. 
/// 
pub fn weight_modifiers(level_sets: &[u32], ndim: usize) -> Result<Vec<f64>, SGError> {
    let num_tensors = level_sets.len() / ndim;

    if ndim == 1 {
        return Ok(vec![1.0; num_tensors]);
    }

    let mut weights = vec![0.0; num_tensors];
    let (sorted_maps, lines_1d) = sort_level_sets(level_sets, ndim)?;

    // Step 1: Initialize weights for the last dimension
    let last_dim_lines = &lines_1d[ndim - 1];
    for i in last_dim_lines.windows(2) {
        weights[i[1] as usize - 1] = 1.0;
    }

    // Step 2: Process dimensions backwards, adjusting weights
    for dim in (0..ndim - 1).rev() {
        for segment in lines_1d[dim].windows(2) {
            let start = segment[0];
            let end = segment[1];

            for i in (start..end - 1).rev() {
                let mut val = weights[sorted_maps[dim][i as usize] as usize];
                for j in i + 1..end {
                    val -= weights[sorted_maps[dim][j as usize] as usize];
                }
                weights[sorted_maps[dim][i as usize] as usize] = val;
            }
        }
    }
    Ok(weights)
}