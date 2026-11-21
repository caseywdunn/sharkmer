//! Target-local node lookup with owned, precharged hash-table buffers.

use super::*;
use ahash::RandomState;
use hashbrown::HashMap;

/// Node identities are never changed or reconstructed by table growth.
/// The table allocation drops before its reservation and owning budget.
pub(in crate::pcr) struct NodeLookup {
    map: HashMap<u64, NodeIndex, RandomState>,
    reservation: Option<MemoryReservation>,
    budget: Option<Arc<ResourceBudget>>,
}

impl NodeLookup {
    /// Construct an empty lookup under an optional shared target ceiling.
    pub(in crate::pcr) fn new(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self {
            map: HashMap::with_hasher(RandomState::new()),
            reservation: None,
            budget: budget.cloned(),
        }
    }

    /// Build lookup for exactly these stable IDs, inheriting graph ownership.
    pub(in crate::pcr) fn from_graph(graph: &AssemblyGraph) -> Result<Self> {
        let mut lookup = Self::new(graph.budget.as_ref());
        for node in graph.node_indices() {
            lookup.insert(graph[node].sub_kmer, node)?;
        }
        Ok(lookup)
    }

    /// Borrow the stable graph identity for an oriented node key.
    pub(in crate::pcr) fn get(&self, key: &u64) -> Option<&NodeIndex> {
        self.map.get(key)
    }

    // hashbrown 0.16.1: power-of-two bucket arrays, 7/8 load, one
    // control byte per bucket and at most 16 trailing control bytes. These
    // entries have size/alignment 16/8; layout tests protect the bound.
    fn allocation_bound(capacity: usize) -> Result<u64> {
        ensure!(capacity >= 14, "Node lookup reservation is too small");
        let entry_bytes = size_of::<(u64, NodeIndex)>();
        ensure!(entry_bytes == 16, "Node lookup entry layout changed");
        let bytes = capacity
            .checked_mul(8)
            .and_then(|n| n.checked_div(7))
            .and_then(usize::checked_next_power_of_two)
            .and_then(|n| n.checked_mul(entry_bytes + 1))
            .and_then(|n| n.checked_add(16))
            .filter(|bytes| *bytes <= isize::MAX as usize - 15)
            .context("Node lookup allocation size overflow")?;
        u64::try_from(bytes).context("Node lookup allocation exceeds supported range")
    }

    fn grow(&mut self) -> Result<()> {
        let Some(budget) = &self.budget else {
            return Ok(());
        };
        let capacity = self
            .map
            .capacity()
            .checked_mul(2)
            .context("Node lookup capacity overflow")?
            .max(14);
        // Keep the old charge alive while hashbrown allocates its replacement.
        let reservation =
            budget.reserve(MemoryClass::TargetGraphs, Self::allocation_bound(capacity)?)?;
        self.map
            .try_reserve(
                capacity
                    .checked_sub(self.map.len())
                    .context("Invalid node lookup capacity")?,
            )
            .map_err(|error| anyhow::anyhow!("Unable to reserve node lookup: {error:?}"))?;
        let actual = u64::try_from(self.map.allocation_size())?;
        // Attach the guard before a later invariant failure can return.
        self.reservation = Some(reservation);
        ensure!(
            self.map.capacity() >= capacity && actual <= Self::allocation_bound(capacity)?,
            "Node lookup allocation differs from its pinned bound"
        );
        if let Some(reservation) = self.reservation.as_mut() {
            reservation.resize(actual)?;
        }
        Ok(())
    }

    /// Insert an identity, reserving overlapping buffers before any growth.
    pub(in crate::pcr) fn insert(&mut self, key: u64, node: NodeIndex) -> Result<()> {
        if self.map.len() == self.map.capacity() && !self.map.contains_key(&key) {
            self.grow()?;
        }
        self.map.insert(key, node);
        Ok(())
    }
}

#[cfg(test)]
impl FromIterator<(u64, NodeIndex)> for NodeLookup {
    fn from_iter<T: IntoIterator<Item = (u64, NodeIndex)>>(iter: T) -> Self {
        Self {
            map: iter.into_iter().collect(),
            reservation: None,
            budget: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{ResourceExhausted, ResourceLimits};

    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn lookup_layout_bound_covers_every_growth_and_releases_owned_buffers() {
        let budget = budget(1 << 20);
        let mut lookup = NodeLookup::new(Some(&budget));
        for key in 0..4096 {
            lookup.insert(key, NodeIndex::new(key as usize)).unwrap();
            let actual = lookup.map.allocation_size() as u64;
            assert!(actual <= NodeLookup::allocation_bound(lookup.map.capacity()).unwrap());
            assert_eq!(budget.snapshot().reserved_memory_bytes, actual);
        }
        drop(lookup);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(NodeLookup::allocation_bound(usize::MAX).is_err());
    }

    #[test]
    fn lookup_resize_refusal_preserves_entries_and_can_retry_after_other_owner_drops() {
        let budget = budget(1000);
        let mut lookup = NodeLookup::new(Some(&budget));
        for key in 0..14 {
            lookup.insert(key, NodeIndex::new(key as usize)).unwrap();
        }
        let actual = lookup.map.allocation_size() as u64;
        let held = budget.reserve(MemoryClass::TargetGraphs, 300).unwrap();
        let before = budget.snapshot().reserved_memory_bytes;
        let error = lookup.insert(14, NodeIndex::new(14)).unwrap_err();
        assert!(error.downcast_ref::<ResourceExhausted>().is_some());
        assert_eq!(lookup.map.len(), 14);
        assert_eq!(lookup.map.allocation_size() as u64, actual);
        assert_eq!(budget.snapshot().reserved_memory_bytes, before);
        for key in 0..14 {
            assert_eq!(lookup.get(&key), Some(&NodeIndex::new(key as usize)));
        }
        drop(held);
        lookup.insert(14, NodeIndex::new(14)).unwrap();
        assert_eq!(lookup.get(&14), Some(&NodeIndex::new(14)));
        drop(lookup);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn lookup_inherits_graph_ceiling_and_stable_ids_including_holes() {
        let budget = budget(1 << 20);
        let mut graph = AssemblyGraph::new(Some(&budget)).unwrap();
        let a = graph
            .add_node(DBNode {
                sub_kmer: 1,
                is_start: true,
                is_end: false,
            })
            .unwrap();
        let b = graph
            .add_node(DBNode {
                sub_kmer: 2,
                is_start: false,
                is_end: false,
            })
            .unwrap();
        let c = graph
            .add_node(DBNode {
                sub_kmer: 3,
                is_start: false,
                is_end: true,
            })
            .unwrap();
        graph.remove_node(b);
        let lookup = NodeLookup::from_graph(&graph).unwrap();
        assert_eq!(lookup.get(&1), Some(&a));
        assert_eq!(lookup.get(&2), None);
        assert_eq!(lookup.get(&3), Some(&c));
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            graph.buffer_bytes().unwrap() + lookup.map.allocation_size() as u64
        );
        drop(graph);
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            lookup.map.allocation_size() as u64
        );
        drop(lookup);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }
}
