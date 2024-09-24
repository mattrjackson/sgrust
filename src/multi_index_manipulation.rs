use std::cmp::Ordering;


///
/// This function sorts the tensors in each dimension. Adapted from TASMANIAN code. 
/// 
fn sort_level_sets(level_sets: &[u32], ndim: usize) -> (Vec<Vec<u32>>, Vec<Vec<u32>>)
{
    let num_levels = level_sets.chunks_exact(ndim).map(|idx|*idx.iter().max().unwrap()).max().unwrap() + 1;
    let num_tensors = level_sets.len() / ndim;
    let mut map: Vec<Vec<u32>> = vec![(0..num_tensors as u32).collect(); ndim];
    let mut lines1d: Vec<Vec<u32>> = vec![Vec::with_capacity(num_levels as usize); ndim];
    let match_outside_dim = |d, a: &[u32], b: &[u32]|
    {
        for j in 0..ndim
        {
            if j != d && a[j] != b[j] { return false; }
        }
        true
    };
    for d in 0..ndim
    {
        map[d].sort_by(|&a, &b|
        {
            let idxa = &level_sets.chunks_exact(ndim).nth(a as usize).unwrap();
            let idxb = &level_sets.chunks_exact(ndim).nth(b as usize).unwrap();
            for j in 0..ndim
            {                    
                if j != d
                {
                    match idxa[j].cmp(&idxb[j])
                    {
                        Ordering::Less => return Ordering::Less,
                        Ordering::Equal => {},
                        Ordering::Greater => return Ordering::Greater,
                    }
                }
            } 
            std::cmp::Ordering::Equal
        });
        if ndim == 1
        {
            lines1d[d].push(0);
            lines1d[d].push(num_tensors as u32);
        }
        else
        {
            let mut c_index = level_sets.chunks_exact(ndim).nth(map[d][0] as usize).unwrap();
            lines1d[d].push(0);
            for i in 1..num_tensors
            {
                let index = level_sets.chunks_exact(ndim).nth(map[d][i] as usize).unwrap();
                if !match_outside_dim(d, c_index, index)
                {
                    lines1d[d].push(i as u32);
                    c_index = index;
                }
            }
            lines1d[d].push(num_tensors as u32);
        }
    }
    (map, lines1d)
}

/// compute the weight modifiers the level sets...
pub(crate) fn weight_modifiers(level_sets: &[u32], ndim: usize) -> Vec<f64>
{
    if ndim == 1
    {
        return vec![1.0; level_sets.len() / ndim];
    }
    let mut weights = vec![0.0; level_sets.len() / ndim];
    let (map, lines1d) = sort_level_sets(level_sets, ndim);
    let last_jobs = &lines1d[ndim-1];
    for i in 0..last_jobs.len()-1
    {
        weights[last_jobs[i+1] as usize - 1] = 1.0;
    }
    for d in (0..=ndim-2).rev()
    {        
        for job in 0..lines1d[d].len() - 1
        {
            for i in (lines1d[d][job]..=lines1d[d][job+1]-2).rev() {
                let mut val = weights[map[d][i as usize] as usize];
                for j in i+1..lines1d[d][job+1]
                {
                    val -= weights[map[d][j as usize] as usize];
                }
                weights[map[d][i as usize] as usize] = val;
            }
        }
    }
    weights
}