//! Benchmark comparing const_generic vs dynamic LinearGrid interpolation performance
//!
//! Run with: cargo bench --bench dynamic_linear_grid --features rayon

use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

// Import const_generic LinearGrid
use sgrust::const_generic::grids::linear_grid::LinearGrid as ConstLinearGrid;
// Import dynamic LinearGrid  
use sgrust::dynamic::grids::linear_grid::LinearGrid as DynamicLinearGrid;
use sgrust::errors::SGError;

/// Build a 2D grid using const_generic version
fn build_const_2d_grid() -> Result<ConstLinearGrid<2, 1>, SGError> {
    let mut grid = ConstLinearGrid::<2, 1>::default();
    grid.sparse_grid_with_boundaries([6; 2])?;
    
    let f = |x: [f64; 2]| [x[0] * x[0] + x[1] * x[1]];
    let values = grid.points().map(f).collect();
    grid.set_values(values)?;
    
    Ok(grid)
}

/// Build a 2D grid using dynamic version
fn build_dynamic_2d_grid() -> Result<DynamicLinearGrid, SGError> {
    let mut grid = DynamicLinearGrid::new(2, 1);
    grid.sparse_grid_with_boundaries(&[6; 2])?;
    
    let f = |x: &[f64]| vec![x[0] * x[0] + x[1] * x[1]];
    grid.update_values(&f);
    
    Ok(grid)
}

/// Build a 4D grid using const_generic version
fn build_const_4d_grid() -> Result<ConstLinearGrid<4, 1>, SGError> {
    let mut grid = ConstLinearGrid::<4, 1>::default();
    grid.sparse_grid_with_boundaries([4; 4])?;
    
    let f = |x: [f64; 4]| [x.iter().map(|v| v * v).sum::<f64>()];
    let values = grid.points().map(f).collect();
    grid.set_values(values)?;
    
    Ok(grid)
}

/// Build a 4D grid using dynamic version
fn build_dynamic_4d_grid() -> Result<DynamicLinearGrid, SGError> {
    let mut grid = DynamicLinearGrid::new(4, 1);
    grid.sparse_grid_with_boundaries(&[4; 4])?;
    
    let f = |x: &[f64]| vec![x.iter().map(|v| v * v).sum::<f64>()];
    grid.update_values(&f);
    
    Ok(grid)
}

/// Build a 6D grid using const_generic version
fn build_const_6d_grid() -> Result<ConstLinearGrid<6, 1>, SGError> {
    let mut grid = ConstLinearGrid::<6, 1>::default();
    grid.sparse_grid_with_boundaries([4; 6])?;
    
    let f = |x: [f64; 6]| [x.iter().map(|v| v * v * v).sum::<f64>()];
    let values = grid.points().map(f).collect();
    grid.set_values(values)?;
    
    Ok(grid)
}

/// Build a 6D grid using dynamic version
fn build_dynamic_6d_grid() -> Result<DynamicLinearGrid, SGError> {
    let mut grid = DynamicLinearGrid::new(6, 1);
    grid.sparse_grid_with_boundaries(&[4; 6])?;
    
    let f = |x: &[f64]| vec![x.iter().map(|v| v * v * v).sum::<f64>()];
    grid.update_values(&f);
    
    Ok(grid)
}

/// Benchmark single-point interpolation for 2D grids
fn bench_2d_single_interpolation(c: &mut Criterion) {
    let const_grid = build_const_2d_grid().unwrap();
    let dynamic_grid = build_dynamic_2d_grid().unwrap();
    
    let x_const = [0.3, 0.5];
    let x_dynamic = [0.3, 0.5];
    let mut result = [0.0];
    
    let mut group = c.benchmark_group("2D Single Interpolation");
    group.bench_function("const_generic", |b| {
        b.iter(|| const_grid.interpolate(x_const).unwrap())
    });
    group.bench_function("dynamic", |b| {
        b.iter(|| dynamic_grid.interpolate(&x_dynamic, &mut result).unwrap())
    });
    group.finish();
    
    println!("2D Grid sizes: const={}, dynamic={}", const_grid.len(), dynamic_grid.len());
}

/// Benchmark single-point interpolation for 4D grids
fn bench_4d_single_interpolation(c: &mut Criterion) {
    let const_grid = build_const_4d_grid().unwrap();
    let dynamic_grid = build_dynamic_4d_grid().unwrap();
    
    let x_const = [0.3, 0.5, 0.2, 0.7];
    let x_dynamic = [0.3, 0.5, 0.2, 0.7];
    let mut result = [0.0];
    
    let mut group = c.benchmark_group("4D Single Interpolation");
    group.bench_function("const_generic", |b| {
        b.iter(|| const_grid.interpolate(x_const).unwrap())
    });
    group.bench_function("dynamic", |b| {
        b.iter(|| dynamic_grid.interpolate(&x_dynamic, &mut result).unwrap())
    });
    group.finish();
    
    println!("4D Grid sizes: const={}, dynamic={}", const_grid.len(), dynamic_grid.len());
}

/// Benchmark single-point interpolation for 6D grids
fn bench_6d_single_interpolation(c: &mut Criterion) {
    let const_grid = build_const_6d_grid().unwrap();
    let dynamic_grid = build_dynamic_6d_grid().unwrap();
    
    let x_const = [0.3, 0.5, 0.2, 0.7, 0.4, 0.8];
    let x_dynamic = [0.3, 0.5, 0.2, 0.7, 0.4, 0.8];
    let mut result = [0.0];
    
    let mut group = c.benchmark_group("6D Single Interpolation");
    group.bench_function("const_generic", |b| {
        b.iter(|| const_grid.interpolate(x_const).unwrap())
    });
    group.bench_function("dynamic", |b| {
        b.iter(|| dynamic_grid.interpolate(&x_dynamic, &mut result).unwrap())
    });
    group.finish();
    
    println!("6D Grid sizes: const={}, dynamic={}", const_grid.len(), dynamic_grid.len());
}

/// Benchmark multi-point interpolation (1000 iterations) for 2D grids
fn bench_2d_multi_interpolation(c: &mut Criterion) {
    let const_grid = build_const_2d_grid().unwrap();
    let dynamic_grid = build_dynamic_2d_grid().unwrap();
    
    let x_const = [0.3, 0.5];
    let x_dynamic = [0.3, 0.5];
    let mut result = [0.0];
    
    let mut group = c.benchmark_group("2D 1000x Interpolation");
    group.bench_function("const_generic", |b| {
        b.iter(|| {
            for _ in 0..1000 {
                let _ = const_grid.interpolate(x_const);
            }
        })
    });
    group.bench_function("dynamic", |b| {
        b.iter(|| {
            for _ in 0..1000 {
                let _ = dynamic_grid.interpolate(&x_dynamic, &mut result);
            }
        })
    });
    group.finish();
}

/// Benchmark multi-point interpolation (1000 iterations) for 6D grids
fn bench_6d_multi_interpolation(c: &mut Criterion) {
    let const_grid = build_const_6d_grid().unwrap();
    let dynamic_grid = build_dynamic_6d_grid().unwrap();
    
    let x_const = [0.3, 0.5, 0.2, 0.7, 0.4, 0.8];
    let x_dynamic = [0.3, 0.5, 0.2, 0.7, 0.4, 0.8];
    let mut result = [0.0];
    
    let mut group = c.benchmark_group("6D 1000x Interpolation");
    group.bench_function("const_generic", |b| {
        b.iter(|| {
            for _ in 0..1000 {
                let _ = const_grid.interpolate(x_const);
            }
        })
    });
    group.bench_function("dynamic", |b| {
        b.iter(|| {
            for _ in 0..1000 {
                let _ = dynamic_grid.interpolate(&x_dynamic, &mut result);
            }
        })
    });
    // Benchmark interpolate_with_state
    group.bench_function("dynamic_with_state", |b| {
        let mut state = dynamic_grid.create_interpolation_state();
        b.iter(|| {
            for _ in 0..1000 {
                let _ = dynamic_grid.interpolate_with_state(&x_dynamic, &mut state, &mut result);
            }
        })
    });
    group.finish();
}

/// Benchmark varying dimensions
fn bench_scaling_by_dimension(c: &mut Criterion) {
    let mut group = c.benchmark_group("Dimension Scaling");
    
    // 2D
    {
        let const_grid = build_const_2d_grid().unwrap();
        let dynamic_grid = build_dynamic_2d_grid().unwrap();
        let x_const = [0.3, 0.5];
        let x_dynamic = [0.3, 0.5];
        let mut result = [0.0];
        
        group.bench_with_input(BenchmarkId::new("const_generic", "2D"), &(), |b, _| {
            b.iter(|| const_grid.interpolate(x_const).unwrap())
        });
        group.bench_with_input(BenchmarkId::new("dynamic", "2D"), &(), |b, _| {
            b.iter(|| dynamic_grid.interpolate(&x_dynamic, &mut result).unwrap())
        });
    }
    
    // 4D
    {
        let const_grid = build_const_4d_grid().unwrap();
        let dynamic_grid = build_dynamic_4d_grid().unwrap();
        let x_const = [0.3, 0.5, 0.2, 0.7];
        let x_dynamic = [0.3, 0.5, 0.2, 0.7];
        let mut result = [0.0];
        
        group.bench_with_input(BenchmarkId::new("const_generic", "4D"), &(), |b, _| {
            b.iter(|| const_grid.interpolate(x_const).unwrap())
        });
        group.bench_with_input(BenchmarkId::new("dynamic", "4D"), &(), |b, _| {
            b.iter(|| dynamic_grid.interpolate(&x_dynamic, &mut result).unwrap())
        });
    }
    
    // 6D
    {
        let const_grid = build_const_6d_grid().unwrap();
        let dynamic_grid = build_dynamic_6d_grid().unwrap();
        let x_const = [0.3, 0.5, 0.2, 0.7, 0.4, 0.8];
        let x_dynamic = [0.3, 0.5, 0.2, 0.7, 0.4, 0.8];
        let mut result = [0.0];
        
        group.bench_with_input(BenchmarkId::new("const_generic", "6D"), &(), |b, _| {
            b.iter(|| const_grid.interpolate(x_const).unwrap())
        });
        group.bench_with_input(BenchmarkId::new("dynamic", "6D"), &(), |b, _| {
            b.iter(|| dynamic_grid.interpolate(&x_dynamic, &mut result).unwrap())
        });
    }
    
    group.finish();
}

criterion_group!(
    benches,
    bench_2d_single_interpolation,
    bench_4d_single_interpolation,
    bench_6d_single_interpolation,
    bench_2d_multi_interpolation,
    bench_6d_multi_interpolation,
    bench_scaling_by_dimension
);
criterion_main!(benches);
