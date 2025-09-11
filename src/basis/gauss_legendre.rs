use std::f64::consts::PI;

use crate::basis::base::Basis;

// Include the pre-computed Gauss-Legendre tables (generated data)
include!("gl_tables.rs");

/// Maximum level of pre-computed Gauss-Legendre data
pub const GL_MAX_LEVEL: u32 = 64;

/// Compute the Legendre polynomial P_n(x) and its derivative using recurrence
fn legendre_and_derivative(n: usize, x: f64) -> (f64, f64) {
    let mut p0 = 1.0;
    let mut p1 = x;
    let mut dp0 = 0.0;
    let mut dp1 = 1.0;

    for k in 2..=n {
        let kf = k as f64;
        let pk = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
        let dpk = ((2.0 * kf - 1.0) * (p1 + x * dp1) - (kf - 1.0) * dp0) / kf;

        p0 = p1;
        p1 = pk;
        dp0 = dp1;
        dp1 = dpk;
    }

    (p1, dp1)
}

/// Compute Gauss-Legendre nodes and weights on the interval (0, 1)
/// Only used for levels beyond GL_MAX_LEVEL
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = Vec::with_capacity(n);
    let mut weights = Vec::with_capacity(n);
    let eps = 1e-14;

    for i in 0..n {
        // Initial guess: Chebyshev nodes (good approximation)
        let theta = PI * (i as f64 + 0.75) / (n as f64 + 0.5);
        let mut x = theta.cos();

        // Newton-Raphson refinement
        for _ in 0..100 {
            let (p, dp) = legendre_and_derivative(n, x);
            let dx = -p / dp;
            x += dx;
            if dx.abs() < eps {
                break;
            }
        }

        // Compute weight
        let (_, dp) = legendre_and_derivative(n, x);
        let w = 2.0 / ((1.0 - x * x) * dp * dp);

        // Map from (-1, 1) to (0, 1)
        let x_mapped = 0.5 * (x + 1.0);
        let w_mapped = 0.5 * w;

        nodes.push(x_mapped);
        weights.push(w_mapped);
    }

    // Sort nodes and weights (not strictly necessary)
    let mut pairs: Vec<(f64, f64)> = nodes.into_iter().zip(weights).collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let (nodes_sorted, weights_sorted): (Vec<_>, Vec<_>) = pairs.into_iter().unzip();
    (nodes_sorted, weights_sorted)
}

/// Access pre-computed Gauss-Legendre data
pub struct GaussLegendreCache;

impl GaussLegendreCache {
    /// Get the node at the given level and index
    #[inline]
    pub fn node(level: u32, index: u32) -> f64 {
        debug_assert!(level <= GL_MAX_LEVEL, "Level {} exceeds GL_MAX_LEVEL {}", level, GL_MAX_LEVEL);
        GL_NODES[GL_OFFSETS[level as usize] + index as usize]
    }

    /// Get the weight at the given level and index
    #[inline]
    pub fn weight(level: u32, index: u32) -> f64 {
        debug_assert!(level <= GL_MAX_LEVEL, "Level {} exceeds GL_MAX_LEVEL {}", level, GL_MAX_LEVEL);
        GL_WEIGHTS[GL_OFFSETS[level as usize] + index as usize]
    }

    /// Get all weights for a level
    #[inline]
    pub fn weights_for_level(level: u32) -> &'static [f64] {
        debug_assert!(level <= GL_MAX_LEVEL, "Level {} exceeds GL_MAX_LEVEL {}", level, GL_MAX_LEVEL);
        let start = GL_OFFSETS[level as usize];
        let end = GL_OFFSETS[level as usize + 1];
        &GL_WEIGHTS[start..end]
    }

    /// Get all nodes for a level
    #[inline]
    pub fn nodes_for_level(level: u32) -> &'static [f64] {
        debug_assert!(level <= GL_MAX_LEVEL, "Level {} exceeds GL_MAX_LEVEL {}", level, GL_MAX_LEVEL);
        let start = GL_OFFSETS[level as usize];
        let end = GL_OFFSETS[level as usize + 1];
        &GL_NODES[start..end]
    }

    /// Get the number of nodes at this level
    #[inline]
    pub fn num_nodes(level: u32) -> usize {
        level as usize
    }
}

#[derive(Copy, Clone)]
pub struct GaussLegendre;

impl Basis for GaussLegendre
{
    fn eval(&self, level: u32, index: u32, x: f64) -> f64 {
        let xl = GaussLegendreCache::node(level, index  - 1);        
        let xr = GaussLegendreCache::node(level, index + 1);
        let xm = GaussLegendreCache::node(level, index);
        if x < xm
        {
            f64::max(0.0, 1.0-(xm-x)/(xm-xl))
        }
        else
        {
            f64::max(0.0, (xr-x)/(xr-xm))
        }
    }

    fn eval_deriv(&self, _level: u32, _index: u32, _x: f64) -> f64 {
        panic!("Derivative not implemented for GaussLegendre Basis.")
    }

    fn degree(&self) -> usize {
        1
    }

    fn integral(&self, level: u32, index: u32) -> f64 {
        GaussLegendreCache::weight(level, index)
    }

    fn basis_type(&self) -> super::base::BasisFunction {
        super::base::BasisFunction::GaussLegendre
    }

    fn node(&self, level: u32, index: u32) -> f64 {
        GaussLegendreCache::node(level, index)
    }

    fn num_nodes(&self, level: u32) -> usize {
        level as usize
    }

    fn get_point(&self, level: u32, index: u32) -> f64 {
        GaussLegendreCache::node(level, index)
    }

    fn weights(&self, level: u32) -> Vec<f64> {
        GaussLegendreCache::weights_for_level(level).to_vec()
    }

    fn nodes(&self, level: u32) -> Vec<f64> {
        GaussLegendreCache::nodes_for_level(level).to_vec()
    }

    fn interpolation_exactness(&self, level: u32) -> u32 {
        if level == 0
        {
            1
        }
        else
        {
            2_u32.pow(level) 
        }
    }

    fn quadrature_exactness(&self, level: u32) -> u32 {
        if level == 0
        {
            1
        }
        else
        {
            2_u32.pow(level) - 1
        }
    }
}

#[test]
fn test_gauss_legendre() {
    let (nodes, weights) = (GaussLegendre.nodes(10), GaussLegendre.weights(10));
    let expected_nodes = vec![0.0130467357414145,0.067468316655508,0.160295215850488,0.283302302935377,0.425562830509185,0.574437169490815,0.716697697064624,0.839704784149512,0.932531683344492,0.986953264258586];
    let expected_weights = vec![0.033335672154344,0.07472567457529,0.109543181257991,0.134633359654998,0.147762112357376,0.147762112357376,0.134633359654998,0.109543181257991,0.07472567457529,0.033335672154344];

    for (n1, n2) in nodes.iter().zip(expected_nodes.iter()) {
        assert!((n1 - n2).abs() < 1e-12);
    }
    for (w1, w2) in weights.iter().zip(expected_weights.iter()) {
        assert!((w1 - w2).abs() < 1e-12);
    }
}