# Change log
<!-- Template:
## [version] — YYYY-MM-DD

### Security
### Added
### Changed
### Deprecated
### Removed
### Fixed
-->

## [0.6.0] — 2025-07-10

### Security
### Added
Created vector-based data structures for `LinearGrid`. This is generally (much) easier to integrate with other languages than the const generics. Note that the dynamic `LinearGrid` is quite a bit slower than the const generic version. Will be working on closing the performance gap, but as things currently stand, generally seeing a 30-50% performance hit from the dynamic `LinearGrid` compared to the const generic variant.
### Changed
The current const-generic based `LinearGrid` has been moved under the `const-generic` module. The combination grid has been moved under the `dynamic` module.
### Deprecated
### Removed
### Fixed

## [0.5.2] — 2025-06-17

### Security
### Added
### Changed
Added an additional parameter to the refinement function to pass in the grid points.
### Deprecated
### Removed
### Fixed


## [0.5.1] — 2025-06-15

### Security
### Added
### Changed
### Deprecated
### Removed
### Fixed
Removed debug section from `iterate_refinable_points` that prevented refinement above level 16.


## [0.5.0] — 2025-06-08

### Security
### Added
1. `RefinementMode` to allow switching between isotropic and anisotropic refinement. 
2. `RefinementOptions` to store refinement options. This includes the refinement threshold, the `RefinementMode`, and an optional parameter to specify the max level limits.
### Changed
The refinement function now takes a `RefinementOptions` input rather than a single floating point threshold.
### Deprecated
### Removed
### Fixed

## [0.4.6] — 2025-05-20

### Security
### Added
### Changed
Modified the refinement functor and removed the threshold method. We now pass in the threshold to the coarsen or refinement calls directly.
### Deprecated
### Removed
### Fixed
Corrected bug in `coarsen` function that on rare occasions could result in the zero index node being removed.

## [0.4.5] — 2025-05-18

### Security
### Added
### Changed
1. Upgraded to bincode 2.x. 
2. Minor performance improvements 
3. Modified refinement and coarsening functors to calculation values for all points in the refinement iteration, rather than individually. This can be significantly faster depending when calling the library outside of Rust.
### Deprecated
### Removed
### Fixed
Corrected bug in `coarsen` function that on rare occasions could result in the zero index node being removed.


## [0.4.4] — 2025-04-29

### Security
### Added
### Changed
### Deprecated
### Removed
### Fixed
Corrected bug in `coarsen` function that on rare occasions could result in the zero index node being removed.

## [0.4.3] — 2025-04-22

### Security
### Added
### Changed
### Deprecated
### Removed
Removed unneeded print statements that were used for debugging.
### Fixed


## [0.4.2] — 2025-02-08

### Security
### Added
### Changed
Modified the adjacency data to also store the left level zero, which can significantly speed up interpolation for smaller grids with boundaries.
### Deprecated
### Removed
### Fixed


## [0.4.1] — 2025-01-27

### Security
### Added
### Changed
Corrected a bug when deserializing grids.
### Deprecated
### Removed
### Fixed

## [0.4.0] — 2025-01-18

### Security
### Added
### Changed
Significant cleanup of data structures following 0.3.x enhancements. This introduces breaking changes with 0.3.x and ealier. Adjacency data is now serialized by default. 
### Deprecated
### Removed
### Fixed

## [0.3.3] — 2025-01-10

### Security
### Added
### Changed
Slightly improve performance for larger grids by removing needless allocations in boundary evalution. 
### Deprecated
### Removed
### Fixed

## [0.3.2] — 2024-12-25

### Security
### Added
Added a new test `check_grid_refinement_parallel` to ensure `update_values_parallel` works as intended.
### Changed
Fixed a bug in `update_values_parallel` that resulted in the newly computed values not being added to the grid. 
### Deprecated
### Removed
### Fixed


## [0.3.1] — 2024-12-25

### Security
### Added
### Changed
### Deprecated
### Removed
Removed asserts that are no longer needed (and were incorrect) in `AdjacencyGridIterator`. 
### Fixed

## [0.3.0] — 2024-12-25

### Security
### Added
### Changed
Reduced memory overhead of `AdjacencyGridIterator` further. This is a breaking change for serializing/deserializing a `LinearGrid`. Overall this is within 2% of the performance of the original methodology rather required 68 bytes of overhead * number of points * number of dimensions vs 20 bytes of overhead * number of points * number of dimensions.
### Deprecated
### Removed
### Fixed

## [0.2.2] — 2024-12-23

### Security
### Added
1. Created new methods, `update_refined_values` and `refine_iteration` to provide an alternate interface to iteratively call the refinement function and update the values.  
2. Created new methods `update_values` and `update_values_parallel` to provide an easier way to update the base grid values. 
### Changed
Reduced memory overhead of `AdjacencyGridIterator` by 2x. This also appears to improve performance by 10-20%. This is a breaking change for serializing/deserializing a `LinearGrid`.
### Deprecated
### Removed
### Fixed

## [0.2.1] — 2024-12-20

### Security
### Added
Created a new `refine_parallel` method for the `SparseGrid` trait. The results in values being filled for the refined nodes in parallel.
### Changed
`RefinementFunctor` must now implement Send + Sync.
### Deprecated
### Removed
Removed the useless sorting in the `refine` function. Presently we do not limit the number of nodes to refine in a given iteration.
### Fixed

## [0.2.0] — 2024-12-16

### Security
### Added
### Changed
Improved performance of `CombinationGrid` interpolation by ~40%. It is now on par with the performance of TASMANIAN.
### Deprecated
### Removed
### Fixed

## [0.1.1] — 2024-12-14

### Security
### Added
### Changed
### Deprecated
### Removed
Removed unused tests file.
### Fixed

## [0.1.0] — 2024-12-14

### Security
### Added
Initial commit of base capabilities
### Changed
### Deprecated
### Removed
### Fixed