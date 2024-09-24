use sgrust::{basis::global::GlobalBasis, grids::combination_grid::{CombinationSparseGrid, GenerationOptions}};

fn integration_example()
{
    let mut grid = CombinationSparseGrid::new(3, vec![GlobalBasis::new(sgrust::basis::global::GlobalBasisType::GaussPatterson); 3]);
    grid.sparse_grid(GenerationOptions{tensor_strategy: sgrust::grids::combination_grid::TensorSelectionStrategy::Level, 
        level_limits: vec![4; 3], ..Default::default()}).unwrap();
    
    let values: Vec<_> = grid.nodes().chunks_exact(3).map(|x|
    {
        x[0]*x[0] + x[1]*x[1] + x[2]*x[2]
    }).collect();

    let integral = grid.integral(&values, 1)[0];
    println!("expected = 1.0, value = {integral}");
}


fn main()
{
    integration_example();
}