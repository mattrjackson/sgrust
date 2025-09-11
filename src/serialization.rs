/// Serialization format options for sparse grid data.
/// 
/// Each format has both compressed (Lz4) and uncompressed variants.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SerializationFormat {
    /// JSON format - human readable, larger size, widest compatibility
    Json,
    /// JSON format with LZ4 compression
    JsonLz4,
    /// Bitcode format - compact binary, good performance
    Bitcode,
    /// Bitcode format with LZ4 compression (default, best balance of size and speed)
    #[default]
    BitcodeLz4,
    /// rkyv format - zero-copy deserialization, best performance
    /// Requires the `rkyv` feature to be enabled
    #[cfg(feature = "rkyv")]
    Rkyv,
    /// rkyv format with LZ4 compression
    /// Requires the `rkyv` feature to be enabled
    #[cfg(feature = "rkyv")]
    RkyvLz4,
}

impl SerializationFormat {
    /// Returns true if this format uses LZ4 compression
    pub fn is_compressed(&self) -> bool {
        match self {
            SerializationFormat::JsonLz4 | SerializationFormat::BitcodeLz4 => true,
            #[cfg(feature = "rkyv")]
            SerializationFormat::RkyvLz4 => true,
            _ => false,
        }
    }

    /// Returns true if this format uses rkyv
    #[cfg(feature = "rkyv")]
    pub fn is_rkyv(&self) -> bool {
        matches!(self, SerializationFormat::Rkyv | SerializationFormat::RkyvLz4)
    }
}

use crate::errors::SGError;
use serde::{de::DeserializeOwned, Serialize};

/// Serialize data to bytes using the specified format (serde-based formats only).
fn serialize_serde<T: Serialize>(data: &T, format: SerializationFormat) -> Result<Vec<u8>, SGError> {
    match format {
        SerializationFormat::Json | SerializationFormat::JsonLz4 => {
            serde_json::to_vec(data).map_err(|_| SGError::SerializationFailed)
        }
        SerializationFormat::Bitcode | SerializationFormat::BitcodeLz4 => {
            bitcode::serialize(data).map_err(|_| SGError::SerializationFailed)
        }
        #[cfg(feature = "rkyv")]
        SerializationFormat::Rkyv | SerializationFormat::RkyvLz4 => {
            Err(SGError::SerializationFailed) // Use serialize_rkyv instead
        }
    }
}

/// Deserialize data from bytes using the specified format (serde-based formats only).
fn deserialize_serde<T: DeserializeOwned>(data: &[u8], format: SerializationFormat) -> Result<T, SGError> {
    match format {
        SerializationFormat::Json | SerializationFormat::JsonLz4 => {
            serde_json::from_slice(data).map_err(|_| SGError::DeserializationFailed)
        }
        SerializationFormat::Bitcode | SerializationFormat::BitcodeLz4 => {
            bitcode::deserialize(data).map_err(|_| SGError::DeserializationFailed)
        }
        #[cfg(feature = "rkyv")]
        SerializationFormat::Rkyv | SerializationFormat::RkyvLz4 => {
            Err(SGError::DeserializationFailed) // Use deserialize_rkyv instead
        }
    }
}

// Re-export rkyv traits for convenience when feature is enabled
#[cfg(feature = "rkyv")]
pub use rkyv::{Archive, Deserialize as RkyvDeserialize, Serialize as RkyvSerialize};

/// High-level serializer type used by rkyv's `to_bytes` function.
#[cfg(feature = "rkyv")]
pub type HighSerializer<'a, E> = rkyv::rancor::Strategy<
    rkyv::ser::Serializer<rkyv::util::AlignedVec, rkyv::ser::allocator::ArenaHandle<'a>, rkyv::ser::sharing::Share>,
    E
>;

/// High-level deserializer type used by rkyv's `from_bytes` function.
#[cfg(feature = "rkyv")]
pub type HighDeserializer<E> = rkyv::rancor::Strategy<rkyv::de::Pool, E>;

/// High-level validator type used by rkyv's `from_bytes` function.
#[cfg(feature = "rkyv")]
pub type HighValidator<'a, E> = rkyv::rancor::Strategy<
    rkyv::validation::Validator<rkyv::validation::archive::ArchiveValidator<'a>, rkyv::validation::shared::SharedValidator>,
    E
>;

/// Serialize data to bytes using rkyv.
/// 
/// The type must derive `rkyv::Archive` and `rkyv::Serialize`.
#[cfg(feature = "rkyv")]
pub fn serialize_rkyv<T>(data: &T) -> Result<Vec<u8>, SGError> 
where 
    T: for<'a> rkyv::Serialize<HighSerializer<'a, rkyv::rancor::Error>>,
{
    rkyv::to_bytes::<rkyv::rancor::Error>(data)
        .map(|v| v.to_vec())
        .map_err(|_| SGError::SerializationFailed)
}

/// Serialize data to bytes using rkyv with LZ4 compression.
#[cfg(feature = "rkyv")]
pub fn serialize_rkyv_lz4<T>(data: &T) -> Result<Vec<u8>, SGError> 
where 
    T: for<'a> rkyv::Serialize<HighSerializer<'a, rkyv::rancor::Error>>,
{
    let bytes = serialize_rkyv(data)?;
    Ok(lz4_flex::compress_prepend_size(&bytes))
}

/// Deserialize data from bytes using rkyv.
/// 
/// The type must derive `rkyv::Archive` and `rkyv::Deserialize`.
#[cfg(feature = "rkyv")]
pub fn deserialize_rkyv<T>(data: &[u8]) -> Result<T, SGError>
where
    T: rkyv::Archive,
    T::Archived: rkyv::Deserialize<T, HighDeserializer<rkyv::rancor::Error>> 
        + for<'a> rkyv::bytecheck::CheckBytes<HighValidator<'a, rkyv::rancor::Error>>,
{
    rkyv::from_bytes::<T, rkyv::rancor::Error>(data)
        .map_err(|_| SGError::DeserializationFailed)
}

/// Deserialize data from bytes using rkyv with LZ4 decompression.
#[cfg(feature = "rkyv")]
pub fn deserialize_rkyv_lz4<T>(data: &[u8]) -> Result<T, SGError>
where
    T: rkyv::Archive,
    T::Archived: rkyv::Deserialize<T, HighDeserializer<rkyv::rancor::Error>>
        + for<'a> rkyv::bytecheck::CheckBytes<HighValidator<'a, rkyv::rancor::Error>>,
{
    let decompressed = lz4_flex::decompress_size_prepended(data)
        .map_err(|_| SGError::LZ4DecompressionFailed)?;
    deserialize_rkyv(&decompressed)
}

/// Serialize data to bytes using the specified format.
/// Applies LZ4 compression if the format variant ends with Lz4.
/// Note: For rkyv formats, use serialize_rkyv or serialize_rkyv_lz4 directly.
pub fn serialize<T: Serialize>(data: &T, format: SerializationFormat) -> Result<Vec<u8>, SGError> {
    let bytes = serialize_serde(data, format)?;
    match format {
        SerializationFormat::JsonLz4 | SerializationFormat::BitcodeLz4 => {
            Ok(lz4_flex::compress_prepend_size(&bytes))
        }
        _ => Ok(bytes),
    }
}

/// Deserialize data from bytes using the specified format.
/// Applies LZ4 decompression if the format variant ends with Lz4.
/// Note: For rkyv formats, use deserialize_rkyv or deserialize_rkyv_lz4 directly.
pub fn deserialize<T: DeserializeOwned>(data: &[u8], format: SerializationFormat) -> Result<T, SGError> {
    match format {
        SerializationFormat::JsonLz4 | SerializationFormat::BitcodeLz4 => {
            let decompressed = lz4_flex::decompress_size_prepended(data)
                .map_err(|_| SGError::LZ4DecompressionFailed)?;
            deserialize_serde(&decompressed, format)
        }
        _ => deserialize_serde(data, format),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(serde::Serialize, serde::Deserialize, Debug, PartialEq)]
    struct TestData {
        values: Vec<f64>,
        name: String,
    }

    #[test]
    fn test_json_roundtrip() {
        let data = TestData {
            values: vec![1.0, 2.0, 3.0],
            name: "test".to_string(),
        };
        
        let bytes = serialize(&data, SerializationFormat::Json).unwrap();
        let result: TestData = deserialize(&bytes, SerializationFormat::Json).unwrap();
        assert_eq!(data, result);
    }

    #[test]
    fn test_json_lz4_roundtrip() {
        let data = TestData {
            values: vec![1.0, 2.0, 3.0],
            name: "test".to_string(),
        };
        
        let bytes = serialize(&data, SerializationFormat::JsonLz4).unwrap();
        let result: TestData = deserialize(&bytes, SerializationFormat::JsonLz4).unwrap();
        assert_eq!(data, result);
    }

    #[test]
    fn test_bitcode_roundtrip() {
        let data = TestData {
            values: vec![1.0, 2.0, 3.0],
            name: "test".to_string(),
        };
        
        let bytes = serialize(&data, SerializationFormat::Bitcode).unwrap();
        let result: TestData = deserialize(&bytes, SerializationFormat::Bitcode).unwrap();
        assert_eq!(data, result);
    }

    #[test]
    fn test_bitcode_lz4_roundtrip() {
        let data = TestData {
            values: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            name: "compressed_test".to_string(),
        };
        
        let bytes = serialize(&data, SerializationFormat::BitcodeLz4).unwrap();
        let result: TestData = deserialize(&bytes, SerializationFormat::BitcodeLz4).unwrap();
        assert_eq!(data, result);
    }

    #[cfg(feature = "rkyv")]
    mod rkyv_tests {
        use super::*;

        #[derive(Debug, PartialEq, rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)]
        #[rkyv(derive(Debug))]
        struct RkyvTestData {
            value: i32,
            count: u64,
        }

        #[test]
        fn test_rkyv_roundtrip() {
            let data = RkyvTestData {
                value: 42,
                count: 100,
            };
            
            let bytes = serialize_rkyv(&data).unwrap();
            let result: RkyvTestData = deserialize_rkyv(&bytes).unwrap();
            assert_eq!(data, result);
        }

        #[test]
        fn test_rkyv_lz4_roundtrip() {
            let data = RkyvTestData {
                value: 123,
                count: 456,
            };
            
            let bytes = serialize_rkyv_lz4(&data).unwrap();
            let result: RkyvTestData = deserialize_rkyv_lz4(&bytes).unwrap();
            assert_eq!(data, result);
        }

        #[test]
        fn test_sparse_grid_rkyv_roundtrip() {
             use crate::const_generic::grids::sparse_grid::SparseGridBase;
             let grid = SparseGridBase::<2, 1>::new();
             
             let bytes = serialize_rkyv(&grid).unwrap();
             let deserialized: SparseGridBase<2, 1> = deserialize_rkyv(&bytes).unwrap();
             
             assert_eq!(grid.storage.len(), deserialized.storage.len());
        }

        #[test]
        fn test_linear_grid_rkyv_roundtrip() {
            use crate::const_generic::grids::linear_grid::LinearGrid;
            let mut grid = LinearGrid::<2, 1>::new();
            grid.base_mut().storage.insert_point(crate::const_generic::storage::GridPoint::new([1, 1], [1, 1], true));

            let bytes = serialize_rkyv(&grid).unwrap();
            let deserialized: LinearGrid<2, 1> = deserialize_rkyv(&bytes).unwrap();

            assert_eq!(grid.base().storage.len(), deserialized.base().storage.len());
        }

        #[test]
        fn test_combination_grid_rkyv_roundtrip() {
            use crate::dynamic::grids::combination_grid::CombinationSparseGrid;
            use crate::basis::global::{GlobalBasis, GlobalBasisType};
            
            let grid = CombinationSparseGrid::new(2, vec![GlobalBasis::new(GlobalBasisType::Linear); 2]);
            
            let bytes = serialize_rkyv(&grid).unwrap();
            let deserialized: CombinationSparseGrid = deserialize_rkyv(&bytes).unwrap();
            
            assert_eq!(grid.ndim(), deserialized.ndim());
        }
    }
}
