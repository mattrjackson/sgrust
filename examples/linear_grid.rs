
use sgrust::{errors::SGError, grids::{linear_grid::LinearGrid, sparse_grid::SparseGrid}, refinement::surplus::SurplusRefinement, storage::linear_grid::BoundingBox};

fn one_d()-> Result<(), SGError>
{
    println!("\nRunning \"one_d\" example\n");
    // Build the 1D grid object, only one value per node.
    let mut grid = LinearGrid::<1,1>::default();
    // build a sparse grid in 1D (this is identical to a full grid in 1D)
    grid.sparse_grid([5]);

    // create our function for evaluating values on the grid
    let f = |x: &[f64; 1]|[x[0].powi(2)];
    // Calculate the values on our grid
    grid.update_values(&f);
    // Compare the values...
    let x = [0.3];
    let mut error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);
    // Let's set our objective to be 1/10 of the current error.
    let refinement = SurplusRefinement;
    println!("Number of points: {}", grid.len());
    grid.refine(&refinement, &f, error/10.0, 10);
    println!("Number of points after refinement: {}", grid.len());
    error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);

    Ok(())
}

///
///  An example of a 2D grid. Here we use refinement to get a more accurate interpolation
/// 
fn two_d()-> Result<(), SGError>
{
    println!("\nRunning \"two_d\" example\n");
    // Build the 2D grid object, only one value per node.
    let mut grid = LinearGrid::<2,1>::default();
    // We are using a 
    grid.full_grid_with_boundaries(5);

    let f = |x: &[f64; 2]|[x[0].powi(2) + x[1].powi(2)];

    // Calculate the values on our grid
    grid.update_values(&f);
    
    // Compare the values...
    let x = [0.3, 0.1];
    let mut error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);
    // Let's set our objective to be 1/10 of the current error.
    let refinement = SurplusRefinement;
    println!("Number of points: {}", grid.len());
    grid.refine(&refinement, &f, error/10.0, 10);
    println!("Number of points after refinement: {}", grid.len());
    error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);
    // now let's coarsen and check the final value...    
    grid.coarsen(&SurplusRefinement, 1e-8);
    println!("Number of points after coarsening: {}", grid.len());
    error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);
    Ok(())
}


///
///  An example of a 6D grid. Here we use refinement to get a more accurate interpolation
/// 
fn six_d()-> Result<(), SGError>
{
    println!("\nRunning \"six_d\" example\n");
    // Build the 2D grid object, only one value per node.
    let mut grid = LinearGrid::<6,1>::default();
    // using a full grid here, but a sparse grid could be used instead...
    grid.full_grid_with_boundaries(2);

    let f = |x: &[f64; 6]|
    {
        let mut r = [0.0];
        (0..6).for_each(|i| {
            r[0] += x[i]*x[i]*x[i];
        });
        r
    };

    grid.update_values(&f);

    // Compare the values...
    let x = [0.3, 0.1, 0.2, 0.1, 0.4, 0.7];
    let mut error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);
    // Let's set our objective to be 1/100 of the current error.
    let refinement = SurplusRefinement;
    println!("Number of points: {}", grid.len());
    grid.refine(&refinement, &f, error/100.0, 10);
    println!("Number of points after refinement: {}", grid.len());
    error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);
    // now let's coarsen and check the final value...    
    grid.coarsen(&SurplusRefinement, 1e-8);
    println!("Number of points after coarsening: {}", grid.len());
    error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);
    Ok(())
}

fn one_d_grid_with_bounding_box()-> Result<(), SGError>
{
    println!("\nRunning \"one_d_grid_with_bounding_box\"\n");
    // Build the 1D grid object, only one value per node.
    let mut grid = LinearGrid::<1,1>::default();
    // build a sparse grid in 1D (this is identical to a full grid in 1D)
    grid.sparse_grid([5]);
    *grid.bounding_box_mut() = BoundingBox::new([0.0], [5.0]);

    // create our function for evaluating values on the grid
    let f = |x: &[f64; 1]|[x[0].powi(2)];

    // Calculate the values on our grid
    grid.update_values(&f);
    // Compare the values...
    let x = [3.0];
    let mut error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);
    // Let's set our objective to be 1/1000 of the current error.
    let refinement = SurplusRefinement;
    println!("Number of points: {}", grid.len());
    grid.refine(&refinement, &f, error/1000.0, 10);
    println!("Number of points after refinement: {}", grid.len());
    error = (grid.interpolate(x)?[0] - f(&x)[0]).abs();
    println!("x={:?}, calculated {}, expected {}. Error={error}", x, grid.interpolate(x)?[0], f(&x)[0]);

    Ok(())
}



fn main()
{
    one_d().unwrap();
    two_d().unwrap();
    six_d().unwrap();
    one_d_grid_with_bounding_box().unwrap();
}