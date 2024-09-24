# Sparse Grids for Rust (SGRust)

Traditionally full grids are used for interpolation or integration of multidimensional functions. Classical sparse grids are a multivariate discretization technique which selects a subset of the full tensors found in a full grid that contribute most to the accuracy of the solution. See https://andreasschaab.com/wp-content/uploads/2021/12/handbook_sparse_grids_in_econ.pdf for more details.

## Usage:

1. `LinearGrid`. This class utilizes linear basis functions, and supports interoplation, integration, as well as adaptive refinement and coarsening. 
2. `CombinationGrid`.  This class currently is primarily useful for UQ applications, but supports Lagrange interpolation and integration based on Clenshaw Curtis or Gauss Patterson. The combination grid also allows tensor selection based on total polynomial exactness, as well as the level sum. The former can allow additional flexibility for the quadrature rules, and can result in more tensors being selected than the simple level sum approach leverages.

## Future goals for this project:

1. Enable refinement/coarsening on the `CombinationGrid`. 
2. Improve interpolation performance on `CombinationGrid`.
3. Reduce memory usage for grid iterators used in the `LinearGrid`.
4. Revise data structures for `LinearGrid` to improve data locality.

## Credits:

Many of the algorithms were adapted from two excellent C++ codes: SG++ and TASMANIAN. The main motivation for creation of this library was to have a pure Rust implementation that didn't require the use of unsafe blocks to access external C++ libraries. For the low-dimensionality cases I work with (<10 dimensions), I have found the `LinearGrid` to be ~10x faster than corresponding capabilities in SG++ or TASMANIAN. This is largely due to a memory vs performance tradeoff (SG++ makes heavy use of hash lookups, while TASMANIAN utilizes a binary search, both of which are slower than storing the neighboring indices directly). Lagrange interpolation on the SGRust `CombinationGrid` is ~20-30% slower than TASMANIAN (likely because of a lack of parallelization of any of the underlying components in the Rust code). 

1. https://github.com/SGpp/SGpp
2. https://github.com/ORNL/TASMANIAN


