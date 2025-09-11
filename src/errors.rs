use std::fmt::Display;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum SGError
{
    OutOfDomain,
    NumberOfPointsAndValuesMismatch,
    KdTreeError,
    LZ4DecompressionFailed,
    ReadBufferFailed,
    WriteBufferFailed,
    SerializationFailed,
    DeserializationFailed,    
    FileIOError,
    InvalidIndex,
    InvalidIteratorSequence
}
impl std::error::Error for SGError {}

impl Display for SGError
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", *self)
    }
}