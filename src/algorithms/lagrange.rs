// Function to compute the Lagrange basis polynomials at a given point x
// points: interpolation points, values: values at those points
#[inline]
pub fn lagrange_weights(x: f64, coeffs: &[f64], points: &[f64]) -> Vec<f64> 
{        
    let mut weights = vec![0.0; points.len()];    
    let mut normalization_factor = 0.0;

    // handle case where point coincides with one of our nodes...
    for (&point, weight) in points.iter().zip(weights.iter_mut())
    {
        if (point-x).abs() < f64::EPSILON
        {
            *weight = 1.0;
            return weights;
        }
    }
    coeffs.iter().zip(points).zip(weights.iter_mut()).for_each(|((&coeff, &xi), weight)| 
    {                        
        *weight = coeff / (x-xi);
        normalization_factor += *weight;
    });    
    weights.iter_mut().for_each(|w|*w  /= normalization_factor);
    weights
}


/// A more efficient means to compute the lagrange weights that minimizes allocations...
#[allow(unused)]
#[inline]
pub fn lagrange_weights_cache(x: f64, coeffs: &[f64], points: &[f64], weights: &mut [f64]) -> usize 
{        
    
    weights.fill(0.0);
    
    let mut normalization_factor = 0.0;

    // handle case where point coincides with one of our nodes...
    for (&point, weight) in points.iter().zip(weights.iter_mut())
    {
        if (point-x).abs() < f64::EPSILON
        {
            *weight = 1.0;
            return points.len();
        }
    }
    coeffs.iter().zip(points).zip(weights.iter_mut()).for_each(|((&coeff, &xi), weight)| 
    {                        
        *weight = coeff / (x-xi);
        normalization_factor += *weight;
    });    
    weights.iter_mut().for_each(|w|*w  /= normalization_factor);
    points.len()
}

pub fn lagrange_coeffs(points: &[f64]) -> Vec<f64>
{
    let mut coeffs = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let mut li = 1.0;
        for j in 0..points.len() {
            if i != j {
                if points[j] == points[i]
                {
                    continue;
                }
                li *= points[i] - points[j];
            }
        }
        coeffs.push(1.0 / li);
    }
    coeffs
}

// Function to compute the Jacobian weights (derivatives of the Lagrange basis functions)
#[allow(unused)]
pub fn jacobian_lagrange_weights<const D: usize>( x: [f64; D], points: &[[f64; D]]) -> Vec<[f64; D]> {
    let n = points.len();
    let mut jacobian_weights = vec![[0.0; D]; n];

    // Loop over each Lagrange polynomial (one per interpolation point)
    for i in 0..n {
        #[allow(clippy::needless_range_loop)]
        for d in 0..D {
            let mut derivative = 0.0;

            // Apply the product rule for the derivative with respect to x_d
            for k in 0..n {
                if i != k {
                    let mut product = 1.0;

                    // Compute the product, skipping j = i and j = k
                    for j in 0..n {
                        if j != i && j != k {
                            product *= (x[d] - points[j][d]) / (points[i][d] - points[j][d]);
                        }
                    }

                    // Apply the derivative term
                    derivative += product / (points[i][d] - points[k][d]);
                }
            }

            jacobian_weights[i][d] = derivative;
        }
    }
    jacobian_weights
}


#[test]
fn test_lagrange_weights()
{
    use crate::tables::clenshaw_curtis_table::cc_nodes;
    let points = cc_nodes(8);
    let weights= lagrange_weights(0.2, &lagrange_coeffs(&points), &points);
    println!("sum={}",weights.iter().zip(&points).map(|(&w,x)|w*x*x).sum::<f64>());
    assert!((1.0-weights.iter().zip(points).map(|(&w,x)|w*x*x).sum::<f64>()/(0.2*0.2)).abs() < 1e-14);    
}

#[test]
fn test_lagrange_derivative()
{
    use crate::tables::clenshaw_curtis_table::cc_nodes;
    let points = cc_nodes(2).iter().map(|x|[*x;1]).collect::<Vec<_>>();
    let derivatives= jacobian_lagrange_weights([0.2], &points);
    assert!((1.0-derivatives.iter().zip(points).map(|(&w,x)|w[0]*x[0]*x[0]).sum::<f64>()/(0.4)).abs() < 1e-14);    
}

// #[test]
// fn test_lagrange_weights_2d()
// {
//     // f(x,y)= x^2+y^2
//     let points = crate::one_dimensional_nodes::clenshaw_curtis_nodes(2).iter().map(|x|[*x;1]).collect::<Vec<_>>();
    
//     let weights= lagrange_weights([0.2, 0.2], &points);
//     assert!((1.0-weights.iter().zip(points).map(|(&w,x)|w*x[0]*x[0]).sum::<f64>()/(0.2*0.2)).abs() < 1e-14);    
// }

// #[test]
// fn test_lagrange_derivative_2d()
// {
//     // f(x,y)= x^2+y^2
//     // f('x) = 2x
//     // f'(y) = 2y
//     let points = crate::one_dimensional_nodes::clenshaw_curtis_nodes(2).iter().map(|x|[*x;1]).collect::<Vec<_>>();
//     let derivatives= jacobian_lagrange_weights([0.2], &points);
//     assert!((1.0-derivatives.iter().zip(points).map(|(&w,x)|w[0]*x[0]*x[0]).sum::<f64>()/(0.4)).abs() < 1e-14);    
// }