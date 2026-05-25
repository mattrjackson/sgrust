# Sparse Grids for Rust (SGRust)

SGRust is a pure Rust sparse-grid library for interpolation, integration, adaptive refinement, and uncertainty-quantification style quadrature.

Traditional full grids become expensive quickly as dimension grows. Sparse grids select a smaller set of tensor-product nodes that usually captures most of the accuracy of a full grid with far fewer points. For more details, see the discussion in this paper: [Handbook on Sparse Grids in Economics](https://andreasschaab.com/wp-content/uploads/2021/12/handbook_sparse_grids_in_econ.pdf).

## Installation

Add SGRust to your `Cargo.toml`:

```toml
[dependencies]
sgrust = "0.8.6"
```

The default feature set includes `rkyv` serialization support. Parallel batch operations are available with the optional `rayon` feature:

```toml
sgrust = { version = "0.8.6", features = ["rayon"] }
```

## Which Grid Should I Use?

Use `const_generic::grids::linear_grid::LinearGrid<D, DIM_OUT>` when the input and output dimensions are known at compile time. This is usually the fastest and most ergonomic choice.

Use `dynamic::grids::linear_grid::LinearGrid` when the dimensions are only known at runtime.

Use `dynamic::grids::combination_grid::CombinationSparseGrid` for global-basis sparse grids, especially quadrature and UQ workflows using Clenshaw-Curtis, Gauss-Patterson, Gauss-Legendre, or custom global rules.

## Example: Const-Generic Linear Interpolation

This creates a 2D linear sparse grid, fits values for `f(x, y) = x^2 + y^2`, and interpolates one point.

```rust
use sgrust::const_generic::grids::linear_grid::LinearGrid;
use sgrust::errors::SGError;

fn main() -> Result<(), SGError> {
    let mut grid = LinearGrid::<2, 1>::new();
    grid.sparse_grid_with_boundaries([6, 6])?;

    grid.update_values(&|x| [x[0] * x[0] + x[1] * x[1]]);

    let y = grid.interpolate([0.25, 0.5])?;
    println!("interpolated value = {}", y[0]);
    println!("grid points = {}", grid.len());

    Ok(())
}
```

## Example: Dynamic Linear Interpolation

The dynamic API stores output values in flat buffers. It is useful when input/output dimensions are runtime values, but comes at the cost of ~30% performance hit.

```rust
use sgrust::dynamic::grids::linear_grid::LinearGrid;
use sgrust::errors::SGError;

fn main() -> Result<(), SGError> {
    let mut grid = LinearGrid::new(2, 1);
    grid.sparse_grid_with_boundaries(&[6, 6])?;

    grid.update_values(&|x| vec![x[0].sin() + x[1].cos()]);

    let mut y = [0.0];
    grid.interpolate(&[0.25, 0.5], &mut y)?;
    println!("interpolated value = {}", y[0]);

    Ok(())
}
```

## Example: Adaptive Refinement With a Closure

If the model function is available inside Rust, `refine` can evaluate new grid points automatically.

```rust
use sgrust::const_generic::algorithms::refinement::RefinementOptions;
use sgrust::const_generic::grids::linear_grid::LinearGrid;
use sgrust::const_generic::refinement::surplus::SurplusRefinement;
use sgrust::errors::SGError;

fn main() -> Result<(), SGError> {
    let mut grid = LinearGrid::<2, 1>::new();
    grid.full_grid_with_boundaries(2)?;

    let model = |x: &[f64; 2]| [1.0 / (0.1 + (x[0] - 0.4).powi(2) + (x[1] - 0.7).powi(2))];
    grid.update_values(&model);

    let functor = SurplusRefinement;
    let options = RefinementOptions::new(0.01);
    grid.refine(&functor, &model, options, 5)?;

    let y = grid.interpolate([0.4, 0.7])?;
    println!("refined grid points = {}", grid.len());
    println!("interpolated value = {}", y[0]);

    Ok(())
}
```

## Example: External or Long-Running Model Values

When values come from another process or a long-running simulation, use `points`, `refine_iteration`, and `update_refined_values`.

```rust
use sgrust::const_generic::algorithms::refinement::RefinementOptions;
use sgrust::const_generic::grids::linear_grid::LinearGrid;
use sgrust::const_generic::refinement::surplus::SurplusRefinement;
use sgrust::errors::SGError;

fn expensive_model(x: &[f64; 2]) -> [f64; 1] {
    [x[0] * x[0] + (4.0 * x[1]).sin()]
}

fn main() -> Result<(), SGError> {
    let mut grid = LinearGrid::<2, 1>::new();
    grid.sparse_grid_with_boundaries([3, 3])?;

    let initial_values = grid.points().map(|x| expensive_model(&x)).collect();
    grid.set_values(initial_values)?;

    let functor = SurplusRefinement;
    let candidates = grid.refine_iteration(&functor, RefinementOptions::new(0.001));

    let new_values: Vec<_> = candidates.iter().map(expensive_model).collect();
    grid.update_refined_values(&new_values, true)?;

    println!("grid points after one refinement pass = {}", grid.len());
    Ok(())
}
```

## Example: Integration With a LinearGrid

`LinearGrid::integrate` integrates the fitted interpolant over the grid domain.

```rust
use sgrust::const_generic::grids::linear_grid::LinearGrid;
use sgrust::errors::SGError;

fn main() -> Result<(), SGError> {
    let mut grid = LinearGrid::<2, 1>::new();
    grid.sparse_grid_with_boundaries([8, 8])?;
    grid.update_values(&|x| [x[0] * x[0] + x[1] * x[1]]);

    let integral = grid.integrate();
    println!("integral over [0, 1]^2 = {}", integral[0]);

    Ok(())
}
```

## Example: CombinationSparseGrid for Quadrature

Combination grids use global basis rules and are useful for quadrature and UQ-style workflows.

```rust
use sgrust::basis::global::{GlobalBasis, GlobalBasisType};
use sgrust::dynamic::grids::combination_grid::{
    CombinationSparseGrid, GenerationOptions, TensorSelectionStrategy,
};
use sgrust::errors::SGError;

fn main() -> Result<(), SGError> {
    let basis = vec![GlobalBasis::new(GlobalBasisType::ClenshawCurtis); 2];
    let mut grid = CombinationSparseGrid::new(2, basis);

    grid.sparse_grid(GenerationOptions {
        tensor_strategy: TensorSelectionStrategy::InterpolationExactness,
        level_limits: vec![5, 5],
        ..Default::default()
    })?;

    let values: Vec<f64> = grid
        .nodes()
        .chunks_exact(2)
        .map(|x| x[0] * x[0] + x[1] * x[1])
        .collect();

    let integral = grid.integral_fast(&values, 1);
    let interpolated = grid.interpolate(&[0.25, 0.5], &values, 1);

    println!("nodes = {}", grid.len());
    println!("integral = {}", integral[0]);
    println!("interpolated value = {}", interpolated[0]);

    Ok(())
}
```

## Example: Serialization

Linear grids can be serialized to buffers or files. `BitcodeLz4` is the default compact binary format.

```rust
use sgrust::const_generic::grids::linear_grid::LinearGrid;
use sgrust::errors::SGError;
use sgrust::serialization::SerializationFormat;

fn main() -> Result<(), SGError> {
    let mut grid = LinearGrid::<2, 1>::new();
    grid.sparse_grid_with_boundaries([4, 4])?;
    grid.update_values(&|x| [x[0] + x[1]]);

    let buffer = grid.write_buffer(SerializationFormat::BitcodeLz4)?;
    let restored = LinearGrid::<2, 1>::read_buffer(&buffer, SerializationFormat::BitcodeLz4)?;

    println!("restored grid points = {}", restored.len());
    Ok(())
}
```

## Notes

For one-dimensional `LinearGrid`s, interpolation uses a fast nodal interval path. Uniform 1D grids use direct interval lookup, while refined or otherwise non-uniform 1D grids use a coordinate-sorted binary-search fallback.

`CombinationSparseGrid` interpolation is global-basis interpolation, not piecewise-linear interpolation. Its 1D path therefore uses Lagrange weights rather than the `LinearGrid` interval shortcut.

## Future Goals

1. Continue improving data locality and cache behavior for `LinearGrid`.
2. Expand adaptive-refinement capabilities for `CombinationSparseGrid`.
3. Add more benchmark-driven examples for common interpolation and quadrature workloads.

## Credits

Many algorithms in this library were adapted from two excellent C++ codes: SG++ and TASMANIAN. The main motivation for this library was to provide a pure Rust implementation that does not require unsafe bindings to external C++ libraries. 

1. https://github.com/SGpp/SGpp
2. https://github.com/ORNL/TASMANIAN
