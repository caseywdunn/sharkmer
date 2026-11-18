//! Fallible ownership of target-local vector, queue and set backing buffers.
//! Element-owned allocations require their own guards; these types charge
//! capacity and resize overlap, not allocator rounding or whole-process RSS.

use super::*;
use ahash::RandomState;
use hashbrown::HashSet;
use std::collections::VecDeque;
use std::hash::Hash;

pub(in crate::pcr) use crate::kmer::budget::collections::BudgetVec;
use crate::kmer::budget::collections::{Charge, bytes};

/// A FIFO queue that retains its charge while wrapped or empty.
pub(in crate::pcr) struct BudgetQueue<T> {
    data: VecDeque<T>,
    charge: Charge,
}

impl<T> BudgetQueue<T> {
    /// Start empty, sharing an optional target ceiling.
    pub(in crate::pcr) fn new(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self {
            data: VecDeque::new(),
            charge: Charge::new(budget),
        }
    }

    /// Reserve old and replacement buffers before growing the queue.
    pub(in crate::pcr) fn push_back(&mut self, item: T) -> Result<()> {
        if self.data.len() == self.data.capacity() {
            let capacity = self
                .data
                .capacity()
                .checked_mul(2)
                .context("Workspace queue capacity overflow")?
                .max(4);
            let reservation = self.charge.reserve(bytes::<T>(capacity)?)?;
            self.data
                .try_reserve_exact(
                    capacity
                        .checked_sub(self.data.len())
                        .context("Invalid workspace queue capacity")?,
                )
                .context("Unable to allocate workspace queue")?;
            self.charge.replace(
                reservation,
                bytes::<T>(self.data.capacity())?,
                bytes::<T>(capacity)?,
            )?;
        }
        self.data.push_back(item);
        Ok(())
    }

    /// Remove the next value without freeing its reusable backing capacity.
    pub(in crate::pcr) fn pop_front(&mut self) -> Option<T> {
        self.data.pop_front()
    }
}

/// A pinned hashbrown set whose immutable API cannot grow its table implicitly.
pub(in crate::pcr) struct BudgetSet<T> {
    data: HashSet<T, RandomState>,
    charge: Charge,
}

impl<T: Eq + Hash> BudgetSet<T> {
    /// Start empty, sharing an optional target ceiling.
    pub(in crate::pcr) fn new(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self {
            data: HashSet::with_hasher(RandomState::new()),
            charge: Charge::new(budget),
        }
    }

    fn allocation_bound(capacity: usize) -> Result<u64> {
        ensure!(capacity >= 14, "Workspace set reservation is too small");
        let buckets = capacity
            .checked_mul(8)
            .and_then(|n| n.checked_div(7))
            .and_then(usize::checked_next_power_of_two)
            .context("Workspace set bucket overflow")?;
        // hashbrown 0.16.1 TableLayout: round data to max(entry alignment,
        // Group::WIDTH), then append bucket controls and a trailing group.
        // Group width is at most 16 on supported targets. Larger alignment is
        // included explicitly rather than assuming 16-byte entries.
        let alignment = std::mem::align_of::<T>().max(16);
        let data = size_of::<T>()
            .checked_mul(buckets)
            .and_then(|n| n.checked_add(alignment - 1))
            .map(|n| n & !(alignment - 1))
            .context("Workspace set data size overflow")?;
        let total = data
            .checked_add(buckets)
            .and_then(|n| n.checked_add(16))
            .filter(|n| *n <= isize::MAX as usize - (alignment - 1))
            .context("Workspace set allocation overflow")?;
        u64::try_from(total).context("Workspace set allocation exceeds supported range")
    }

    fn grow(&mut self) -> Result<()> {
        if !self.charge.is_enabled() {
            return Ok(());
        }
        let capacity = self
            .data
            .capacity()
            .checked_mul(2)
            .context("Workspace set capacity overflow")?
            .max(14);
        let reservation = self.charge.reserve(Self::allocation_bound(capacity)?)?;
        self.data
            .try_reserve(
                capacity
                    .checked_sub(self.data.len())
                    .context("Invalid workspace set capacity")?,
            )
            .map_err(|error| anyhow::anyhow!("Unable to allocate workspace set: {error:?}"))?;
        self.charge.replace(
            reservation,
            u64::try_from(self.data.allocation_size())?,
            Self::allocation_bound(capacity)?,
        )?;
        ensure!(
            self.data.capacity() >= capacity,
            "Workspace set capacity differs from its pinned bound"
        );
        Ok(())
    }

    /// Insert without uncharged growth; return whether this key was new.
    pub(in crate::pcr) fn insert(&mut self, key: T) -> Result<bool> {
        if self.data.len() == self.data.capacity() && !self.data.contains(&key) {
            self.grow()?;
        }
        Ok(self.data.insert(key))
    }

    /// Test exact membership without allocating.
    pub(in crate::pcr) fn contains(&self, key: &T) -> bool {
        self.data.contains(key)
    }

    /// Borrow keys without exposing a growing mutable table.
    pub(in crate::pcr) fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    /// Report the number of keys.
    pub(in crate::pcr) fn len(&self) -> usize {
        self.data.len()
    }

    /// Report whether no keys are present.
    pub(in crate::pcr) fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Clear values while retaining and charging reusable capacity.
    pub(in crate::pcr) fn clear(&mut self) {
        self.data.clear();
    }
}

impl<T: Eq + Hash> Default for BudgetSet<T> {
    fn default() -> Self {
        Self::new(None)
    }
}

impl<T: Eq + Hash + std::fmt::Debug> std::fmt::Debug for BudgetSet<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}

#[cfg(test)]
impl<T: Eq + Hash> From<ahash::AHashSet<T>> for BudgetSet<T> {
    fn from(source: ahash::AHashSet<T>) -> Self {
        Self {
            data: source.into_iter().collect(),
            charge: Charge::new(None),
        }
    }
}

/// Compact membership indexed by stable graph node IDs, including holes.
/// Only backing words grow; clear retains reusable capacity and its charge.
pub(in crate::pcr) struct BudgetNodeSet {
    words: BudgetVec<u64>,
    count: usize,
}

impl BudgetNodeSet {
    /// Start empty under the target's ceiling.
    pub(in crate::pcr) fn new(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self {
            words: BudgetVec::new(budget),
            count: 0,
        }
    }

    /// Insert an ID, charging any new word allocation before setting its bit.
    pub(in crate::pcr) fn insert(&mut self, node: NodeIndex) -> Result<bool> {
        let word = node.index() / 64;
        while self.words.len() <= word {
            self.words.push(0)?;
        }
        let mask = 1u64 << (node.index() % 64);
        if self.words[word] & mask != 0 {
            return Ok(false);
        }
        self.words[word] |= mask;
        self.count += 1;
        Ok(true)
    }

    /// Test an ID without allocating words for absent indices.
    pub(in crate::pcr) fn contains(&self, node: &NodeIndex) -> bool {
        self.words
            .get(node.index() / 64)
            .is_some_and(|word| word & (1u64 << (node.index() % 64)) != 0)
    }

    /// Count distinct IDs, independent of holes or duplicate insertions.
    pub(in crate::pcr) fn len(&self) -> usize {
        self.count
    }

    /// Forget IDs while keeping the charged word buffer for reuse.
    pub(in crate::pcr) fn clear(&mut self) {
        self.words.fill(0);
        self.count = 0;
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
    fn vector_refusal_retry_and_reuse_keep_correct_capacity_charges() {
        let budget = budget(100);
        let mut values = BudgetVec::new(Some(&budget));
        for value in 0..4u64 {
            values.push(value).unwrap();
        }
        let held = budget.reserve(MemoryClass::TargetGraphs, 8).unwrap();
        assert!(
            values
                .push(4)
                .unwrap_err()
                .downcast_ref::<ResourceExhausted>()
                .is_some()
        );
        assert_eq!(&*values, &[0, 1, 2, 3]);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 40);
        drop(held);
        values.push(4).unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 64);
        assert_eq!(budget.snapshot().peak_reserved_memory_bytes, 96);
        values.clear();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 64);
        drop(values);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(bytes::<u64>(usize::MAX).is_err());
        let mut zero = BudgetVec::new(Some(&budget));
        zero.push(()).unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn queue_wraparound_and_refused_growth_preserve_fifo_and_live_charge() {
        let budget = budget(100);
        let mut queue = BudgetQueue::new(Some(&budget));
        for value in 0..4u64 {
            queue.push_back(value).unwrap();
        }
        assert_eq!(queue.pop_front(), Some(0));
        assert_eq!(queue.pop_front(), Some(1));
        queue.push_back(4).unwrap();
        queue.push_back(5).unwrap();
        let held = budget.reserve(MemoryClass::TargetGraphs, 8).unwrap();
        assert!(
            queue
                .push_back(6)
                .unwrap_err()
                .downcast_ref::<ResourceExhausted>()
                .is_some()
        );
        assert_eq!(
            queue.data.iter().copied().collect::<Vec<_>>(),
            vec![2, 3, 4, 5]
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, 40);
        drop(held);
        queue.push_back(6).unwrap();
        for value in 2..7 {
            assert_eq!(queue.pop_front(), Some(value));
        }
        assert_eq!(queue.pop_front(), None);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 64);
        drop(queue);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn set_layout_bounds_cover_key_sizes_alignment_and_resize_overlap() {
        #[derive(PartialEq, Eq, Hash)]
        #[repr(align(64))]
        struct Aligned([u64; 8]);
        fn check<T: Eq + Hash>(values: impl IntoIterator<Item = T>) {
            let budget = budget(4 << 20);
            let mut set = BudgetSet::new(Some(&budget));
            for value in values {
                set.insert(value).unwrap();
                let actual = set.data.allocation_size() as u64;
                assert!(actual <= BudgetSet::<T>::allocation_bound(set.data.capacity()).unwrap());
                assert_eq!(budget.snapshot().reserved_memory_bytes, actual);
            }
            set.clear();
            assert!(set.is_empty());
            assert!(budget.snapshot().reserved_memory_bytes > 0);
            drop(set);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(BudgetSet::<T>::allocation_bound(usize::MAX).is_err());
        }
        check(0..=255u8);
        check(0..4096u32);
        check(0..4096u64);
        check((0..4096u64).map(|n| Aligned([n; 8])));
        let budget = budget(400);
        let mut set = BudgetSet::new(Some(&budget));
        for n in 0..14u64 {
            set.insert(n).unwrap();
        }
        let before = budget.snapshot().reserved_memory_bytes;
        assert!(
            set.insert(14)
                .unwrap_err()
                .downcast_ref::<ResourceExhausted>()
                .is_some()
        );
        assert_eq!(set.len(), 14);
        assert!((0..14).all(|n| set.contains(&n)));
        assert_eq!(budget.snapshot().reserved_memory_bytes, before);
        drop(set);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn indexed_node_flags_handle_holes_duplicates_refusal_and_reuse() {
        let budget = budget(120);
        let mut flags = BudgetNodeSet::new(Some(&budget));
        assert!(flags.insert(NodeIndex::new(130)).unwrap());
        assert!(!flags.insert(NodeIndex::new(130)).unwrap());
        assert!(!flags.contains(&NodeIndex::new(129)));
        let held = budget.reserve(MemoryClass::TargetGraphs, 40).unwrap();
        assert!(
            flags
                .insert(NodeIndex::new(257))
                .unwrap_err()
                .downcast_ref::<ResourceExhausted>()
                .is_some()
        );
        assert_eq!(flags.len(), 1);
        assert!(flags.contains(&NodeIndex::new(130)));
        assert!(!flags.contains(&NodeIndex::new(257)));
        drop(held);
        assert!(flags.insert(NodeIndex::new(257)).unwrap());
        assert_eq!(flags.len(), 2);
        flags.clear();
        assert_eq!(flags.len(), 0);
        assert!(!flags.contains(&NodeIndex::new(130)));
        assert!(flags.insert(NodeIndex::new(130)).unwrap());
        drop(flags);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn backing_reservation_outlives_element_destruction() {
        struct Witness {
            budget: Arc<ResourceBudget>,
            seen: Arc<std::sync::atomic::AtomicBool>,
        }
        impl Drop for Witness {
            fn drop(&mut self) {
                self.seen.store(
                    self.budget.snapshot().reserved_memory_bytes > 0,
                    std::sync::atomic::Ordering::SeqCst,
                );
            }
        }
        let budget = budget(4096);
        let seen = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let mut values = BudgetVec::new(Some(&budget));
        values
            .push(Witness {
                budget: Arc::clone(&budget),
                seen: Arc::clone(&seen),
            })
            .unwrap();
        drop(values);
        assert!(seen.load(std::sync::atomic::Ordering::SeqCst));
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }
}
