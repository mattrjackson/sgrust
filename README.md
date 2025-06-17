# Sparse Grids for Rust (SGRust)

Traditionally full grids are used for interpolation or integration of multidimensional functions. Classical sparse grids are a multivariate discretization technique which selects a subset of the full tensors found in a full grid that contribute most to the accuracy of the solution. See https://andreasschaab.com/wp-content/uploads/2021/12/handbook_sparse_grids_in_econ.pdf for more details.

## Usage:

1. `LinearGrid`. This class utilizes linear basis functions, and supports interpolation, integration, as well as adaptive refinement and coarsening. See more details below...
2. `CombinationGrid`.  This class currently is primarily useful for UQ applications (e.g. integration over global basis vectors), but supports Lagrange interpolation and integration based on Clenshaw Curtis or Gauss Patterson. The combination grid also allows tensor selection based on total polynomial exactness, as well as the level sum. The former can allow additional flexibility for the quadrature rules, and can result in more tensors being selected than the simple level sum approach leverages. 

For `LinearGrid`, there are effectively two ways to handle updating values and refining the grid based upon them:
1. Use Closures to set the base grid values, and update values during refinement. This is the most elegant way to approach the problem, but requires that the model creating the values to be linked to the Rust library (which may not be possible in all scenarios). 
2. Use Vectors to retrieve points (e.g. `points()`, `refine_iteration`), and then set values using `set_values` and ` update_refined_values`. This allows manually loading values from potentially long-running simulations that will not be directly linked to the sparse grid library.

## Future goals for this project:

1. Enable refinement/coarsening on the `CombinationGrid`. 
2. Revise data structures for `LinearGrid` to improve data locality.

## Credits:

Many of the algorithms were adapted from two excellent C++ codes: SG++ and TASMANIAN. The main motivation for creation of this library was to have a pure Rust implementation that didn't require the use of unsafe blocks to access external C++ libraries. 

1. https://github.com/SGpp/SGpp
2. https://github.com/ORNL/TASMANIAN




