use criterion::{criterion_group, criterion_main, Criterion};
use sgrust::{basis::global::{GlobalBasis, GlobalBasisType}, dynamic::grids::combination_grid::{CombinationSparseGrid, GenerationOptions, TensorSelectionStrategy}};

fn build_grid() -> (CombinationSparseGrid, Vec<f64>)
{
    let ndim = 2;
    let mut grid = CombinationSparseGrid::new(ndim, vec![
        GlobalBasis{basis_type: GlobalBasisType::ClenshawCurtis, custom_rule: None}; ndim]);
    let options = GenerationOptions{tensor_strategy: TensorSelectionStrategy::Level, level_limits: vec![8; ndim], ..Default::default()}; 
    grid.sparse_grid(options).unwrap();
    let mut values = Vec::with_capacity(grid.len());
    for x in grid.nodes().chunks_exact(ndim)
    {                
        values.push(x[0] * x[0] + x[1]*x[1]);
    } 
    (grid, values)
}
fn interpolate(grid: &CombinationSparseGrid, values: &[f64])
{
    for _ in 0..1e3 as usize
    {
        let _ = grid.interpolate(&[0.2, 0.2], values, 1)[0];
    }
}

fn run_case(c: &mut Criterion)
{

    let (grid,values) = build_grid();
    c.bench_function("combigrid", |b|b.iter(||interpolate(&grid, &values)));
}

criterion_group!(benches, run_case);
criterion_main!(benches);