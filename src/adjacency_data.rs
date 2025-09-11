use std::ops::{Deref, DerefMut};

use serde::{Deserialize, Deserializer, Serialize, Serializer};


#[derive(Default, Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
pub(crate) struct NodeAdjacencyData
{
    pub(crate) zero_index: usize,
    pub(crate) data: Vec<NodeAdjacency>,
}

impl Deref for NodeAdjacencyData
{
    type Target=Vec<NodeAdjacency>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl DerefMut for NodeAdjacencyData
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}


#[derive(Default, Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize))]
pub(crate) struct NodeAdjacency
{
    pub(crate) inner: NodeAdjacencyInner,
    pub(crate) left_zero: u32,
}


use bitfield_struct::bitfield;

#[bitfield(u128)]
#[derive(PartialEq, Eq)]
pub(crate) struct NodeAdjacencyInner {
    #[bits(25)]
    pub(crate) left: i64,
    #[bits(25)]
    pub(crate) right: i64,
    #[bits(25)]
    pub(crate) up: i64,
    #[bits(25)]
    pub(crate) down_left: i64,
    #[bits(25)]
    pub(crate) down_right: i64,
    #[bits(1)]
    pub(crate) has_left: bool,
    #[bits(1)]
    pub(crate) has_right: bool,
    #[bits(1)]
    pub(crate) has_parent: bool,
}
impl NodeAdjacencyInner
{ 
    #[inline]
    pub fn is_complete(&self) -> bool {
        self.has_left() && self.has_right() && self.has_left_child() && self.has_right_child() && self.has_parent()
    }
    
    #[inline]
    pub fn has_left_child(&self) -> bool {
        self.down_left() != 0
    }
    
    #[inline]
    pub fn has_right_child(&self) -> bool {
        self.down_right() != 0
    }
}

// Serialize Implementation
impl Serialize for NodeAdjacencyInner {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // Serialize as a single `u128` field
        serializer.serialize_u128(self.0)
    }
}

// Deserialize Implementation
impl<'de> Deserialize<'de> for NodeAdjacencyInner {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let raw = u128::deserialize(deserializer)?;
        Ok(Self(raw))
    }
}

// rkyv implementations for NodeAdjacencyInner (bitfield type)
#[cfg(feature = "rkyv")]
mod rkyv_impl {
    use super::NodeAdjacencyInner;
    
    impl rkyv::Archive for NodeAdjacencyInner {
        type Archived = rkyv::rend::u128_le;
        type Resolver = ();

        fn resolve(&self, _resolver: Self::Resolver, out: rkyv::Place<Self::Archived>) {
            out.write(rkyv::rend::u128_le::from_native(self.0));
        }
    }

    impl<S: rkyv::rancor::Fallible + ?Sized> rkyv::Serialize<S> for NodeAdjacencyInner {
        fn serialize(&self, _serializer: &mut S) -> Result<Self::Resolver, S::Error> {
            Ok(())
        }
    }

    impl<D: rkyv::rancor::Fallible + ?Sized> rkyv::Deserialize<NodeAdjacencyInner, D> for rkyv::rend::u128_le {
        fn deserialize(&self, _deserializer: &mut D) -> Result<NodeAdjacencyInner, D::Error> {
            Ok(NodeAdjacencyInner(self.to_native()))
        }
    }
}
