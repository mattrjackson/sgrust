#[derive(Copy, Clone, Debug)]
pub enum SparseGridRule
{    
    /// \brief Classic nested rule using Chebyshev nodes with very low Lebesgue constant.
    ClenshawCurtis,   
    /// \brief Nested rule that is optimized for integration, probably the best integration rule in more than 2 dimensions.
    GaussPatterson,   
}
