use criterion::{criterion_group, criterion_main, Criterion};
use sgrust::{errors::SGError, const_generic::grids::linear_grid::LinearGrid};

fn build_six_d_grid() -> Result<LinearGrid<6,1>, SGError>
{
    // Build the 2D grid object, only one value per node.
    let mut grid = LinearGrid::<6,1>::default();
    // using a full grid here, but a sparse grid could be used instead...
    grid.sparse_grid_with_boundaries([7;6]);

    let f = |x: [f64; 6]|
    {
        let mut r = [0.0];
        (0..6).for_each(|i| {
            r[0] += x[i]*x[i]*x[i];
        });
        r
    };

    let values = grid.points().map(f).collect();
    grid.set_values(values)?;    
    
    Ok(grid)
}

fn build_twelve_d_grid() -> Result<LinearGrid<12,1>, SGError>
{
    // Build the 2D grid object, only one value per node.
    let mut grid = LinearGrid::<12,1>::default();
    // using a full grid here, but a sparse grid could be used instead...
    grid.sparse_grid_with_boundaries([2;12]);

    let f = |x: [f64; 12]|
    {
        let mut r = [0.0];
        (0..12).for_each(|i| {
            r[0] += x[i]*x[i]*x[i];
        });
        r
    };

    let values = grid.points().map(f).collect();
    grid.set_values(values)?;    
    
    Ok(grid)
}

fn six_d(grid: &LinearGrid<6,1>)-> Result<(), SGError>
{
    // Compare the values...
    let x = [0.3, 0.1, 0.2, 0.1, 0.4, 0.7];
    for _ in 0..1000
    {
        let _value = grid.interpolate(x);
    }
    
    Ok(())
}

fn run_six_d(c: &mut Criterion)
{
    let grid = build_six_d_grid().unwrap();
    c.bench_function("6d", |b|b.iter(||six_d(&grid).unwrap()));
}
#[cfg(feature="rayon")]
fn six_d_immutable_graph(grid: &ImmutableLinearGrid<f32, 6,1>)-> Result<(), SGError>
{
    // Compare the values...
    let x = [[0.3, 0.1, 0.2, 0.1, 0.4, 0.7]; 1000];
    let _value = grid.interpolate_batch(&x);
    
    Ok(())
}
#[cfg(feature="rayon")]
fn twelve_d_immutable_graph(grid: &ImmutableLinearGrid<f32, 12,1>)-> Result<(), SGError>
{
    // Compare the values...
    let x = [[0.3, 0.1, 0.2, 0.1, 0.4, 0.7, 0.2, 0.4, 0.3, 0.6, 0.8, 0.1]; 1000];
    let _value = grid.interpolate_batch(&x);
    
    Ok(())
}
#[cfg(feature="rayon")]
fn run_six_d_immutable_graph(c: &mut Criterion)
{
    let grid: ImmutableLinearGrid<f32, 6, 1> = build_six_d_grid().unwrap().into();
    c.bench_function("6d_immmutable_grid", |b|b.iter(||six_d_immutable_graph(&grid).unwrap()));
}
#[cfg(feature="rayon")]
fn run_twelve_d_immutable_graph(c: &mut Criterion)
{
    let grid: ImmutableLinearGrid<f32, 12, 1> = build_twelve_d_grid().unwrap().into();
    c.bench_function("12d", |b|b.iter(||twelve_d_immutable_graph(&grid).unwrap()));
}
criterion_group!(benches_serial, run_six_d);
#[cfg(feature="rayon")]
criterion_group!(benches, run_six_d, run_six_d_immutable_graph, run_twelve_d_immutable_graph);

criterion_main!(benches_serial);
