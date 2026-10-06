
use crate::const_generic::storage::{GridPoint, SparseGridData};
use super::grid_iterator::GridIteratorT;
///
/// This iterator uses array access to retrieve neighbors more quickly (~10x speedup relative to hash queries)
/// at the expense of pretty significant memory overhead. Eventually the data behind building the iterator
/// will be an optional format
/// 
pub(crate) struct AdjacencyGridIterator<'a, const D: usize>
{
    pub(crate) seq: usize,
    storage: &'a SparseGridData<D>,
}


impl<'a, const D: usize> AdjacencyGridIterator<'a, D>
{

    /// Compute the index of the left boundary node (e.g. the left zero-level node) for the given dimension.
    #[inline]
    pub fn compute_lzero(&self, dim: usize) -> Option<u32> {      
        let index = self.seq;
        let offset = self.offset(dim);  
        let index = self.storage.adjacency_data.left_zero[offset + index];
        if index == u32::MAX {
            None
        }
        else
        {
            Some(index)  
        }
    }
    
    /// Compute the index of the right boundary node (e.g. the right zero-level node) for the given dimension.
    #[inline]
    pub fn compute_rzero(&self, dim: usize) -> Option<u32> {
        let index = self.seq;
        let offset = self.offset(dim);  
        let index = self.storage.adjacency_data.right_zero[offset + index];
        if index == u32::MAX {
            None
        }
        else
        {
            Some(index)  
        }
    }
    pub(crate) fn new(storage: &'a SparseGridData<D>) -> Self
    {
        Self { seq: 0, storage }
    }
    #[inline(always)]
    fn offset(&self, dim: usize) -> usize
    {
        dim * self.storage.len()
    }
    

    
    #[inline]
    #[allow(unused)]
    pub(crate) fn has_left_leaf(&self, dim: usize) -> bool
    {
        self.storage.adjacency_data[self.offset(dim) + self.seq].has_left_child()
    }
    #[inline]
    #[allow(unused)]
    pub(crate) fn has_right_leaf(&self, dim: usize) -> bool
    {
        self.storage.adjacency_data[self.offset(dim) + self.seq].has_right_child()
    } 
  
}

impl<const D: usize> GridIteratorT<D> for AdjacencyGridIterator<'_, D>
{
    #[inline]
    fn reset_to_level_zero(&mut self) -> bool
    {
        self.seq = self.storage.adjacency_data.zero_index;
        true
    }
    #[inline]
    fn reset_to_left_level_zero(&mut self, dim: usize) -> bool
    {
        if let Some(index) = self.compute_lzero(dim)
        {
            self.seq = index as usize;
            true
        }
        else
        {
            false
        }
    }
    #[inline]
    fn reset_to_right_level_zero(&mut self, dim: usize) -> bool
    {
        if let Some(index) = self.compute_rzero(dim)
        {
            self.seq = index as usize;
            true
        }
        else
        {
            false
        }
    }

    #[inline]
    fn reset_to_level_one(&mut self, dim: usize) -> bool
    {
        let index = self.storage.adjacency_data[self.offset(dim) + self.seq].level_one();
        if index == u32::MAX
        {
            return false;
        }
        self.seq = index as usize;
        true
    }
    #[inline]
    fn left_child(&mut self, dim: usize) -> bool
    {
        let adj = &self.storage.adjacency_data[self.offset(dim) + self.seq];
        if adj.has_left_child()
        {
            let index = (self.seq as i64 + adj.down_left()) as usize;
            self.seq = index;
            true
        }
        else
        {
            false
        }
    }
    #[inline]
    fn right_child(&mut self, dim: usize) -> bool
    {
        let adj = &self.storage.adjacency_data[self.offset(dim) + self.seq];
        if adj.has_right_child()
        {
            let index = (self.seq as i64 + adj.down_right()) as usize;
            self.seq = index;
            true
        }
        else
        {
            false
        }
    }
   
    #[inline]
    fn is_leaf(&self) -> bool
    {
        // Leaf state belongs to the whole point, as in HashMapGridIterator.
        // Recursive evaluation changes dimensions on entry and return; the
        // last move's adjacency leaf flag cannot stop another dimension's traversal.
        // At level zero, traversal also uses level_one rather than child links.
        self.storage[self.seq].is_leaf()
    }
    
    fn index(&self) ->  Option<usize> {
        Some(self.seq)
    }
    
    fn up(&mut self, dim: usize) -> bool
    {
        let offset = self.offset(dim);
        if self.storage.adjacency_data[offset + self.seq].has_parent()
        {
            let index = (self.seq as i64 + self.storage.adjacency_data[offset + self.seq].up()) as usize;
            self.seq = index;
            true
        }
        else
        {
            false
        }
    }
    
    fn is_inner_point(&self) -> bool {
        self.storage[self.seq].is_inner_point()
    }
    
    fn point(&self) -> &GridPoint<D> {
        &self.storage[self.seq]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::const_generic::{
        algorithms::refinement::{RefinementFunctor, RefinementOptions},
        grids::linear_grid::LinearGrid,
        iterators::grid_iterator::HashMapGridIterator,
        storage::PointIterator,
    };

    struct RefineBoundary;
    impl RefinementFunctor<3, 1> for RefineBoundary {
        fn eval(&self, points: PointIterator<3>, _: &[[f64; 1]], _: &[[f64; 1]]) -> Vec<f64> {
            points
                .map(|p| if p == [0.25, 0.5, 1.0] { 1.0 } else { 0.0 })
                .collect()
        }
    }

    #[test]
    fn cached_moves_and_leaf_state_match_hash_iterator_after_asymmetric_refinement() {
        let mut grid = LinearGrid::<3, 1>::new();
        grid.full_grid_with_boundaries(2).unwrap();
        grid.update_values(&|p| [p.iter().map(|v| v * v).sum()]);
        grid.refine_iteration(&RefineBoundary, RefinementOptions::new(0.1));
        grid.update_values(&|p| [p.iter().map(|v| v * v).sum()]);
        let storage = grid.storage();
        // Check every cached adjacency against independent hash lookup. A failed
        // cached move leaves the iterator in place; a hash move can become invalid.
        for (seq, point) in storage.nodes().iter().enumerate() {
            for dim in 0..3 {
                for movement in 0..6 {
                    let mut cached = AdjacencyGridIterator::new(storage);
                    cached.seq = seq;
                    let mut hashed = HashMapGridIterator::new(storage);
                    hashed.set_index(*point);
                    assert_eq!(cached.is_leaf(), hashed.is_leaf());
                    let moved = match movement {
                        0 => (cached.left_child(dim), hashed.left_child(dim)),
                        1 => (cached.right_child(dim), hashed.right_child(dim)),
                        2 => (cached.up(dim), hashed.up(dim)),
                        3 => (
                            cached.reset_to_left_level_zero(dim),
                            hashed.reset_to_left_level_zero(dim),
                        ),
                        4 => (
                            cached.reset_to_right_level_zero(dim),
                            hashed.reset_to_right_level_zero(dim),
                        ),
                        _ => (
                            cached.reset_to_level_one(dim),
                            hashed.reset_to_level_one(dim),
                        ),
                    };
                    assert_eq!(
                        moved.0, moved.1,
                        "point {point:?}, dim {dim}, move {movement}"
                    );
                    if moved.0 {
                        assert_eq!(cached.index(), hashed.index());
                        assert_eq!(
                            cached.is_leaf(),
                            hashed.is_leaf(),
                            "point {point:?}, dim {dim}, move {movement}"
                        );
                    }
                }
            }
        }
        // Recreate the return from recursive dimension-2 evaluation: it resets
        // z to the left boundary, which has no z children but does have a y root.
        let seq = storage
            .points()
            .position(|p| p == [0.125, 1.0, 1.0])
            .unwrap();
        let mut cached = AdjacencyGridIterator::new(storage);
        cached.seq = seq;
        assert!(cached.reset_to_left_level_zero(2));
        let adj = &storage.adjacency_data[2 * storage.len() + cached.seq];
        assert!(adj.is_leaf());
        assert!(!cached.is_leaf());
        assert!(cached.reset_to_level_one(1));
        assert_eq!(cached.point().unit_coordinate(), [0.125, 0.5, 0.0]);
    }
}
