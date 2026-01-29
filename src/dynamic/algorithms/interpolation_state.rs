//! Interpolation state for zero-allocation repeated interpolation.
//!
//! This module provides `InterpolationState` which can be reused across multiple
//! interpolation calls to avoid per-call setup overhead.

use crate::dynamic::iterators::grid_iterator_cache::AdjacencyGridIterator;
use crate::dynamic::iterators::dynamic_grid_iterator::GridIteratorT;
use crate::dynamic::storage::SparseGridData;

/// Maximum supported dimensions for stack-based arrays
const MAX_DIM: usize = 32;

/// Reusable state for interpolation to eliminate per-call overhead.
/// 
/// Create once from a grid and reuse for many interpolation calls.
/// This avoids recreating the iterator and scratch arrays each call.
/// 
/// # Example
/// ```ignore
/// let state = grid.create_interpolation_state();
/// for point in points {
///     grid.interpolate_with_state(&mut state, &point, &mut result)?;
/// }
/// ```
pub struct InterpolationState<'a> {
    /// The iterator used for tree traversal
    pub(crate) iterator: AdjacencyGridIterator<'a>
}

impl<'a> InterpolationState<'a> {
    /// Create a new interpolation state from storage data
    #[inline]
    pub fn new(storage: &'a SparseGridData) -> Self {
        Self {
            iterator: AdjacencyGridIterator::new(storage),
        }
    }
    
    /// Reset the state for a new interpolation (called internally)
    #[inline]
    pub(crate) fn reset(&mut self) {
        self.iterator.reset_to_level_zero();
    }
}
