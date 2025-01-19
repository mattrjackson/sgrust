use std::ops::{Deref, DerefMut};

use serde::{Deserialize, Deserializer, Serialize, Serializer};


#[derive(Default, Clone, Debug, Serialize, Deserialize)]
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



use bitfield_struct::bitfield;

#[bitfield(u128)]
pub(crate) struct NodeAdjacency {
    #[bits(30)]
    pub(crate) left: i64,
    #[bits(30)]
    pub(crate) right: i64,
    #[bits(30)]
    pub(crate) up: i64,
    #[bits(30)]
    pub(crate) down: i64,
    #[bits(1)]
    pub(crate) has_left: bool,
    #[bits(1)]
    pub(crate) has_right: bool,
    #[bits(1)]
    pub(crate) has_parent: bool,
    #[bits(1)]
    pub(crate) has_left_child: bool,
    #[bits(1)]
    pub(crate) has_right_child: bool,
    #[bits(1)]
    pub(crate) has_child: bool,
    #[bits(2)]
    _other: u8
}
impl NodeAdjacency
{ 
    #[inline]
    pub fn is_complete(&self) -> bool {
        self.has_left() && self.has_right() && self.has_left_child() && self.has_right_child() && self.has_parent()
    }
}

// Serialize Implementation
impl Serialize for NodeAdjacency {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // Serialize as a single `u128` field
        serializer.serialize_u128(self.0)
    }
}

// Deserialize Implementation
impl<'de> Deserialize<'de> for NodeAdjacency {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let raw = u128::deserialize(deserializer)?;
        Ok(Self(raw))
    }
}
