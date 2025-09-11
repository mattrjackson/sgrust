use std::f64::consts::PI;

use rustfft::{FftPlanner, num_complex::Complex};
///
/// Compute Clenshaw-Curtis weights using approach by J. Waldvogel (2013)
/// "Fast construction of the Fejer and Clenshaw-Curtis quadrature rules". 
/// 
pub(crate) fn cc_weights(level: u32) -> Vec<f64> {
    let n = (1_u32 << level) as usize;
    if level == 0
    {
        return vec![1.0];
    }
    let mut n_vals = Vec::new();
    for i in 0..n
    {
        let val = 2 * i + 1;
        if val >= n
        {
            break;
        }
        n_vals.push(val as f64);
    }
    let l = n_vals.len();
    let  m = n - l;    
    let mut v0: Vec<f64> = n_vals.iter()
        .map(|&i| 2.0 / (i * (i - 2.0)))
        .collect();
    v0.push(1.0 / n_vals[n_vals.len()-1]);
    v0.extend(vec![0.0; m]);

    let mut v2 = vec![0.0; v0.len()-1];
    let end = v0.len();
    for i in 0..v2.len()
    {
        v2[i] = -v0[i] - v0[end - i - 1];
    }

    let mut planner = FftPlanner::new();
    let mut weights = vec![Complex::new(-1.0, 0.0); n];
    weights[l] += n as f64;
    weights[m] += n as f64;

    let g_scale = (n * n - 1 + n % 2) as f64;
    // Perform inverse FFT using rustfft
    for (g,w) in weights.iter_mut().zip(v2)
    {
        g.re = g.re / g_scale + w;
    }    
    let fft =   planner.plan_fft_inverse(n);
    fft.process(&mut weights);
    for w in weights.iter_mut()
    {
        w.re *= 0.5 / n as f64;
    }    
    weights.push(weights[0]);
    weights.iter().map(|x|x.re).collect()
}

/// Return Clenshaw-Curtis nodes over a (0,1) interval
pub(crate) fn cc_nodes(level: u32) -> Vec<f64>
{
    let n =  if level == 0 {1} else { 2_u32.pow(level) + 1 } ;
    let mut nodes = vec![0.0; n as usize];    
    if level > 0
    {
        #[allow(clippy::needless_range_loop)]
        for i in 0..n as usize
        {
            nodes[i] = 0.5*(1.0 + f64::cos(PI*(n as usize - 1 - i) as f64 / (n - 1) as f64));
        }
    }
    else 
    {
        nodes[0] = 0.5;
    }
    nodes
}
pub(crate) fn cc_num_nodes(level: u32) -> u32
{
    if level == 0 {1} else { 2_u32.pow(level) + 1 }
}

#[test]
fn check_weights() 
{
    for i in 0..2
    {
        let _ = cc_weights(i);
    }
    let level = 2;
    let cc_weights = cc_weights(level);
    
    println!("Number of Curtis Weights: {}", cc_weights.len());
    println!("Clenshaw-Curtis Weights: {:?}", cc_weights);
    assert!((1.0-cc_weights.iter().sum::<f64>()/1.0).abs() < 1e-15);

   

    // Weights computed from CLENSHAW_CURTIS_RULE by J. Burkardt. These are computed over a (-1,+1) interval, so we multiply
    // these values by 0.5 to correct weighting for our (0,1) interval.
    let weights5 = [0.06666666666666668, 0.5333333333333333, 0.7999999999999999, 0.5333333333333334, 0.06666666666666668];
    let nodes  = cc_nodes(level);
    for i in 0..5
    {
        println!("{},{}",  nodes[i], cc_weights[i]);
        assert!((1.0-cc_weights[i]/(0.5*weights5[i])).abs() < 1e-15);
    }
    
    
}


