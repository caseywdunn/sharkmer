//! Exact v4 count-store interface and in-memory conformance candidate.
//!
//! The slot budget covers owned table buffers, including both sides of a
//! resize. Callers must separately reserve process/allocator overhead, input
//! buffers, graphs and evidence before using this as a production backend.

use super::budget::{MemoryClass, MemoryReservation, ResourceBudget, ResourceLimits};
use super::identity::{CanonicalKey, KmerSize, OrientedKey};
use anyhow::{Context, Result, bail, ensure};
use std::num::NonZeroU64;
use std::sync::Arc;

#[cfg(feature = "ahashmap")]
type TableHasher = ahash::RandomState;
#[cfg(feature = "fxhashmap")]
type TableHasher = std::hash::BuildHasherDefault<rustc_hash::FxHasher>;
#[cfg(feature = "ahashmap")]
const TABLE_HASHER_NAME: &str = "ahash";
#[cfg(feature = "fxhashmap")]
const TABLE_HASHER_NAME: &str = "fxhash";

/// Identity and cardinality of a finalized exact count store.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CountMetadata {
    pub format_version: u32,
    pub k: KmerSize,
    pub count_width_bits: u8,
    pub backend: &'static str,
    pub hash_algorithm: &'static str,
    pub canonicalization: &'static str,
    pub preprocessing: &'static str,
    pub overflow_policy: &'static str,
    pub selected_input_sha256: [u8; 32],
    pub distinct_keys: u64,
    pub occurrences: u64,
    pub max_lookup_batch: usize,
}

/// Shared read-only boundary for memory and disk-backed count stores.
///
/// Absence, a positive exact count, and an access failure are separate results.
/// Implementations verify the entire canonical key and reject mismatched k.
pub trait CountStore: Send + Sync {
    /// Return the immutable metadata bound at finalization.
    fn metadata(&self) -> &CountMetadata;

    /// Look up an exact canonical key, rejecting incompatible identities.
    fn lookup(&self, key: CanonicalKey) -> Result<Option<NonZeroU64>>;

    /// Resolve a directional graph key against canonical counts. Backends may
    /// avoid reverse-complement work on an exact forward hit.
    fn lookup_oriented(&self, key: OrientedKey) -> Result<Option<NonZeroU64>> {
        self.lookup(key.canonical())
    }

    /// Visit every stored key once. Stop and propagate a visitor/access error.
    fn visit(&self, visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>) -> Result<()>;

    /// Resolve a bounded batch in input order, preserving duplicate queries.
    ///
    /// Capacity/identity errors occur before changing output. An access error
    /// may leave partial output; callers must ignore output on any error.
    fn lookup_batch(&self, keys: &[CanonicalKey], output: &mut [Option<NonZeroU64>]) -> Result<()> {
        ensure!(
            keys.len() == output.len(),
            "Count lookup batch/output lengths differ"
        );
        ensure!(
            keys.len() <= self.metadata().max_lookup_batch,
            "Count lookup batch exceeds its reservation"
        );
        for key in keys {
            ensure!(
                key.size() == self.metadata().k,
                "Count lookup k differs from store k"
            );
        }
        for (key, destination) in keys.iter().zip(output.iter_mut()) {
            *destination = self.lookup(*key)?;
        }
        Ok(())
    }
}

/// A non-copying minimum-count view; missing and filtered keys return None.
pub struct CountView<'a> {
    store: &'a dyn CountStore,
    minimum: NonZeroU64,
}

impl<'a> CountView<'a> {
    /// Construct a threshold view without changing or scanning its store.
    pub fn new(store: &'a dyn CountStore, minimum: NonZeroU64) -> Self {
        Self { store, minimum }
    }

    /// Resolve one key under the threshold while preserving access errors.
    pub fn lookup(&self, key: CanonicalKey) -> Result<Option<NonZeroU64>> {
        Ok(self
            .store
            .lookup(key)?
            .filter(|count| *count >= self.minimum))
    }

    /// Visit only admitted entries, without collecting a second global table.
    pub fn visit(
        &self,
        visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>,
    ) -> Result<()> {
        self.store.visit(&mut |key, count| {
            if count >= self.minimum {
                visitor(key, count)?;
            }
            Ok(())
        })
    }
}

#[derive(Clone, Copy, Default)]
struct Slot {
    key: u64,
    // Zero means vacant. An all-A key is valid because vacancy uses count.
    count: u64,
}

struct Table {
    k: KmerSize,
    slots: Vec<Slot>,
    hasher: TableHasher,
    len: usize,
    occurrences: u64,
    slot_budget_bytes: usize,
    budget: Arc<ResourceBudget>,
    slot_reservation: Option<MemoryReservation>,
    peak_slot_bytes: usize,
    max_lookup_batch: usize,
}

impl Table {
    fn bytes(capacity: usize) -> Result<usize> {
        capacity
            .checked_mul(std::mem::size_of::<Slot>())
            .context("Count slot byte size overflow")
    }

    fn home(&self, bits: u64, capacity: usize) -> usize {
        std::hash::BuildHasher::hash_one(&self.hasher, bits) as usize & (capacity - 1)
    }

    fn find(&self, bits: u64) -> Result<usize> {
        let capacity = self.slots.len();
        ensure!(capacity.is_power_of_two(), "Invalid count table capacity");
        let mut index = self.home(bits, capacity);
        for _ in 0..capacity {
            let slot = self.slots[index];
            if slot.count == 0 || slot.key == bits {
                return Ok(index);
            }
            index = (index + 1) & (capacity - 1);
        }
        bail!("Count table has no vacant slot")
    }

    fn grow(&mut self) -> Result<()> {
        let capacity = if self.slots.is_empty() {
            16
        } else {
            self.slots
                .len()
                .checked_mul(2)
                .context("Count capacity overflow")?
        };
        let old_bytes = Self::bytes(self.slots.capacity())?;
        let peak = old_bytes
            .checked_add(Self::bytes(capacity)?)
            .context("Count resize byte size overflow")?;
        ensure!(
            peak <= self.slot_budget_bytes,
            "Count table resize needs {peak} slot bytes; budget is {}",
            self.slot_budget_bytes
        );
        let mut reservation = self.budget.reserve(
            MemoryClass::CountSlots,
            u64::try_from(Self::bytes(capacity)?)?,
        )?;
        let mut slots = Vec::new();
        slots
            .try_reserve_exact(capacity)
            .context("Unable to reserve count table slots")?;
        let actual_peak = old_bytes
            .checked_add(Self::bytes(slots.capacity())?)
            .context("Allocated count resize byte size overflow")?;
        ensure!(
            actual_peak <= self.slot_budget_bytes,
            "Allocated count slots exceed their reservation"
        );
        reservation.resize(u64::try_from(Self::bytes(slots.capacity())?)?)?;
        slots.resize(capacity, Slot::default());
        for slot in &self.slots {
            if slot.count == 0 {
                continue;
            }
            let mut index = self.home(slot.key, capacity);
            while slots[index].count != 0 {
                index = (index + 1) & (capacity - 1);
            }
            slots[index] = *slot;
        }
        self.peak_slot_bytes = self.peak_slot_bytes.max(actual_peak);
        self.slots = slots;
        self.slot_reservation = Some(reservation);
        Ok(())
    }

    fn check_key(&self, key: CanonicalKey) -> Result<()> {
        ensure!(key.size() == self.k, "Count key k differs from store k");
        Ok(())
    }
}

/// Mutable exact counter candidate. Finalization moves its allocation.
pub struct MemoryCountBuilder {
    table: Table,
}

impl MemoryCountBuilder {
    /// Reserve table buffers within a caller-provided slot budget.
    pub fn new(k: KmerSize, slot_budget_bytes: usize, max_lookup_batch: usize) -> Result<Self> {
        let bytes = u64::try_from(slot_budget_bytes)?;
        let mut classes = [0; 8];
        classes[MemoryClass::CountSlots as usize] = bytes;
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
        })?;
        Self::with_budget(k, budget, max_lookup_batch)
    }

    /// Share live memory accounting with other count/graph/evidence owners.
    pub fn with_budget(
        k: KmerSize,
        budget: Arc<ResourceBudget>,
        max_lookup_batch: usize,
    ) -> Result<Self> {
        ensure!(
            max_lookup_batch > 0,
            "Count lookup batch reservation must be positive"
        );
        let slot_budget_bytes = usize::try_from(budget.class_limit(MemoryClass::CountSlots))?;
        let mut table = Table {
            k,
            slots: Vec::new(),
            hasher: TableHasher::default(),
            len: 0,
            occurrences: 0,
            slot_budget_bytes,
            budget,
            slot_reservation: None,
            peak_slot_bytes: 0,
            max_lookup_batch,
        };
        table.grow()?;
        Ok(Self { table })
    }

    /// Add a positive exact count. Validation/overflow/capacity errors preserve
    /// previously accepted counts and totals; the caller must fail publication.
    pub fn add(&mut self, key: CanonicalKey, increment: NonZeroU64) -> Result<()> {
        self.table.check_key(key)?;
        let mut index = self.table.find(key.bits())?;
        let old = self.table.slots[index].count;
        let count = old
            .checked_add(increment.get())
            .context("Exact k-mer count overflow")?;
        let total = self
            .table
            .occurrences
            .checked_add(increment.get())
            .context("Exact k-mer occurrence total overflow")?;
        if old == 0 {
            let length = self
                .table
                .len
                .checked_add(1)
                .context("Distinct k-mer cardinality overflow")?;
            if length > self.table.slots.len() / 4 * 3 {
                self.table.grow()?;
                index = self.table.find(key.bits())?;
            }
            self.table.len = length;
        }
        self.table.slots[index] = Slot {
            key: key.bits(),
            count,
        };
        self.table.occurrences = total;
        Ok(())
    }

    /// Bind the counted input digest and move into a shared immutable store.
    pub fn finalize(self, selected_input_sha256: [u8; 32]) -> MemoryCountStore {
        let metadata = CountMetadata {
            format_version: 1,
            k: self.table.k,
            count_width_bits: 64,
            backend: "open_addressed_memory",
            hash_algorithm: TABLE_HASHER_NAME,
            canonicalization: "min_forward_reverse_complement",
            preprocessing: "uppercase_acgtn_n_breaks_windows",
            overflow_policy: "error_preserve_previous_counts",
            selected_input_sha256,
            distinct_keys: self.table.len as u64,
            occurrences: self.table.occurrences,
            max_lookup_batch: self.table.max_lookup_batch,
        };
        MemoryCountStore {
            table: self.table,
            metadata,
        }
    }
}

/// Immutable in-memory store implementing the same fallible API as disk stores.
pub struct MemoryCountStore {
    table: Table,
    metadata: CountMetadata,
}

impl MemoryCountStore {
    /// Number of owned hash slots (including vacant slots).
    pub fn capacity(&self) -> usize {
        self.table.slots.len()
    }

    /// Bytes currently owned by the finalized slot buffer, excluding overhead.
    pub fn slot_bytes(&self) -> usize {
        self.table.slots.capacity() * std::mem::size_of::<Slot>()
    }

    /// Peak owned slot bytes, including simultaneous old/new resize buffers.
    /// This excludes allocator/process overhead and is not an RSS measurement.
    pub fn peak_slot_bytes(&self) -> usize {
        self.table.peak_slot_bytes
    }
}

impl CountStore for MemoryCountStore {
    fn metadata(&self) -> &CountMetadata {
        &self.metadata
    }

    fn lookup(&self, key: CanonicalKey) -> Result<Option<NonZeroU64>> {
        self.table.check_key(key)?;
        Ok(NonZeroU64::new(
            self.table.slots[self.table.find(key.bits())?].count,
        ))
    }

    fn lookup_oriented(&self, key: OrientedKey) -> Result<Option<NonZeroU64>> {
        // A single canonical probe avoids an extra random miss for reverse
        // orientations; typed reverse complementation uses constant work.
        self.lookup(key.canonical())
    }

    fn visit(&self, visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>) -> Result<()> {
        for slot in &self.table.slots {
            if let Some(count) = NonZeroU64::new(slot.count) {
                // Insertion accepts only validated CanonicalKey identities;
                // visiting these immutable slots need not canonicalize again.
                visitor(CanonicalKey::from_stored(slot.key, self.metadata.k), count)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::identity::KmerWindows;
    use std::collections::BTreeMap;

    fn positive(n: u64) -> NonZeroU64 {
        NonZeroU64::new(n).unwrap()
    }
    fn key(bits: u64, k: usize) -> CanonicalKey {
        CanonicalKey::new(bits, KmerSize::new(k).unwrap()).unwrap()
    }

    #[test]
    fn live_graph_reservations_control_counter_resize_without_losing_counts() {
        let mut classes = [0; 8];
        classes[MemoryClass::CountSlots as usize] = 900;
        classes[MemoryClass::TargetGraphs as usize] = 500;
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 1000,
            process_overhead_bytes: 100,
            class_bytes: classes,
            temporary_disk_bytes: 0,
        })
        .unwrap();
        let mut builder =
            MemoryCountBuilder::with_budget(KmerSize::new(6).unwrap(), Arc::clone(&budget), 8)
                .unwrap();
        let graph = budget.reserve(MemoryClass::TargetGraphs, 500).unwrap();
        let mut keys = std::collections::BTreeSet::new();
        for bits in 0..4096 {
            keys.insert(key(bits, 6));
            if keys.len() == 13 {
                break;
            }
        }
        let keys = keys.into_iter().collect::<Vec<_>>();
        for key in &keys[..12] {
            builder.add(*key, positive(1)).unwrap();
        }
        assert!(builder.add(keys[12], positive(1)).is_err());
        assert_eq!(builder.table.occurrences, 12);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 856);
        drop(graph);
        builder.add(keys[12], positive(1)).unwrap();
        let store = builder.finalize([0; 32]);
        assert_eq!(store.lookup(keys[12]).unwrap(), Some(positive(1)));
        assert_eq!(store.metadata().occurrences, 13);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 612);
        assert_eq!(budget.snapshot().peak_reserved_memory_bytes, 868);
        drop(store);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 100);
    }

    #[test]
    fn empty_finalization_preserves_policy_and_has_no_false_keys() {
        let builder = MemoryCountBuilder::new(KmerSize::new(31).unwrap(), 4096, 4).unwrap();
        let store = builder.finalize([9; 32]);
        assert_eq!(store.lookup(key(0, 31)).unwrap(), None);
        assert_eq!(store.metadata().distinct_keys, 0);
        assert_eq!(store.metadata().occurrences, 0);
        assert_eq!(store.metadata().count_width_bits, 64);
        assert_eq!(store.metadata().hash_algorithm, TABLE_HASHER_NAME);
        assert_eq!(store.metadata().backend, "open_addressed_memory");
        assert_eq!(
            store.metadata().overflow_policy,
            "error_preserve_previous_counts"
        );
        assert_eq!(
            store.metadata().canonicalization,
            "min_forward_reverse_complement"
        );
        assert_eq!(
            store.metadata().preprocessing,
            "uppercase_acgtn_n_breaks_windows"
        );
        store.lookup_batch(&[], &mut []).unwrap();
        let mut called = false;
        store
            .visit(&mut |_, _| {
                called = true;
                Ok(())
            })
            .unwrap();
        assert!(!called);
    }

    #[test]
    fn counted_sequences_match_independent_scalar_counts_across_growth() {
        let k = KmerSize::new(5).unwrap();
        let mut builder = MemoryCountBuilder::new(k, 1 << 20, 16).unwrap();
        let mut expected = BTreeMap::<u64, u64>::new();
        // Enumerate literal 5-mers independently of the rolling extractor.
        for bits in 0..1024u64 {
            let bases = (0..5)
                .rev()
                .map(|i| b"ACGT"[((bits >> (2 * i)) & 3) as usize])
                .collect::<Vec<_>>();
            let reverse = bases.iter().rev().fold(0u64, |x, b| {
                (x << 2)
                    | match b {
                        b'A' => 3,
                        b'C' => 2,
                        b'G' => 1,
                        b'T' => 0,
                        _ => unreachable!(),
                    }
            });
            let canonical = bits.min(reverse);
            let increment = 1 + bits % 7;
            *expected.entry(canonical).or_default() += increment;
            for window in KmerWindows::new(&bases, 5).unwrap() {
                builder
                    .add(window.unwrap().canonical, positive(increment))
                    .unwrap();
            }
        }
        let store = builder.finalize([3; 32]);
        let mut observed = BTreeMap::new();
        store
            .visit(&mut |key, count| {
                assert!(observed.insert(key.bits(), count.get()).is_none());
                Ok(())
            })
            .unwrap();
        assert_eq!(observed, expected);
        assert_eq!(store.metadata().occurrences, expected.values().sum());
        assert_eq!(store.metadata().distinct_keys, expected.len() as u64);
        assert_eq!(store.metadata().selected_input_sha256, [3; 32]);
        assert!(store.peak_slot_bytes() <= 1 << 20);
    }

    #[test]
    fn hash_collisions_never_create_false_membership_and_zero_key_is_valid() {
        let mut builder = MemoryCountBuilder::new(KmerSize::new(6).unwrap(), 4096, 8).unwrap();
        let target = builder.table.home(0, 16);
        let mut colliders = Vec::new();
        for bits in 0..4096 {
            let candidate = key(bits, 6);
            if !colliders.contains(&candidate) && builder.table.home(candidate.bits(), 16) == target
            {
                colliders.push(candidate);
            }
            if colliders.len() == 4 {
                break;
            }
        }
        assert_eq!(colliders.len(), 4);
        assert_eq!(colliders[0], key(0, 6));
        for (i, k) in colliders[..3].iter().enumerate() {
            builder.add(*k, positive(i as u64 + 1)).unwrap();
        }
        let store = builder.finalize([0; 32]);
        for (i, k) in colliders[..3].iter().enumerate() {
            assert_eq!(store.lookup(*k).unwrap(), Some(positive(i as u64 + 1)));
        }
        assert_eq!(store.lookup(colliders[3]).unwrap(), None);
        for (i, key) in colliders[..3].iter().enumerate() {
            let oriented = OrientedKey::new(key.bits(), key.size()).unwrap();
            assert_eq!(
                store.lookup_oriented(oriented).unwrap(),
                Some(positive(i as u64 + 1))
            );
            assert_eq!(
                store
                    .lookup_oriented(oriented.reverse_complement())
                    .unwrap(),
                Some(positive(i as u64 + 1))
            );
        }
        let absent = OrientedKey::new(colliders[3].bits(), colliders[3].size()).unwrap();
        assert_eq!(store.lookup_oriented(absent).unwrap(), None);
        assert_eq!(
            store.lookup_oriented(absent.reverse_complement()).unwrap(),
            None
        );
        assert!(
            store
                .lookup_oriented(OrientedKey::new(0, KmerSize::new(5).unwrap()).unwrap())
                .is_err()
        );
        assert!(store.lookup(key(0, 5)).is_err());
    }

    #[test]
    fn exact_overflow_and_wrong_k_leave_counts_and_totals_intact() {
        let mut builder = MemoryCountBuilder::new(KmerSize::new(3).unwrap(), 4096, 4).unwrap();
        builder.add(key(0, 3), positive(u64::MAX)).unwrap();
        assert!(builder.add(key(0, 3), positive(1)).is_err());
        assert!(builder.add(key(1, 3), positive(1)).is_err());
        assert!(builder.add(key(0, 2), positive(1)).is_err());
        let store = builder.finalize([0; 32]);
        assert_eq!(store.lookup(key(0, 3)).unwrap(), Some(positive(u64::MAX)));
        assert_eq!(store.lookup(key(1, 3)).unwrap(), None);
        assert_eq!(store.metadata().occurrences, u64::MAX);
        assert_eq!(store.metadata().distinct_keys, 1);
    }

    #[test]
    fn resize_budget_includes_old_and_new_allocation_and_failed_growth_preserves_data() {
        // 32 slots would fit alone, but 16+32 during growth would exceed512.
        let mut builder = MemoryCountBuilder::new(KmerSize::new(6).unwrap(), 512, 8).unwrap();
        let mut accepted = Vec::new();
        for bits in 0..4096 {
            let k = key(bits, 6);
            if accepted.contains(&k) {
                continue;
            }
            if accepted.len() == 12 {
                assert!(builder.add(k, positive(1)).is_err());
                break;
            }
            builder.add(k, positive(1)).unwrap();
            accepted.push(k);
        }
        let pointer = builder.table.slots.as_ptr();
        let store = builder.finalize([0; 32]);
        assert_eq!(store.table.slots.as_ptr(), pointer);
        assert_eq!(store.metadata().distinct_keys, 12);
        assert_eq!(store.metadata().occurrences, 12);
        for k in accepted {
            assert_eq!(store.lookup(k).unwrap(), Some(positive(1)));
        }
        assert_eq!(store.peak_slot_bytes(), 256);
        assert!(MemoryCountBuilder::new(KmerSize::new(3).unwrap(), 255, 4).is_err());
        assert!(MemoryCountBuilder::new(KmerSize::new(3).unwrap(), 1024, 0).is_err());
    }

    #[test]
    fn bounded_batch_filter_and_visitor_errors_obey_contract() {
        let mut builder = MemoryCountBuilder::new(KmerSize::new(4).unwrap(), 4096, 3).unwrap();
        builder.add(key(0, 4), positive(2)).unwrap();
        builder.add(key(1, 4), positive(5)).unwrap();
        let store = builder.finalize([0; 32]);
        let mut output = [None; 3];
        store
            .lookup_batch(&[key(1, 4), key(0, 4), key(1, 4)], &mut output)
            .unwrap();
        assert_eq!(
            output,
            [Some(positive(5)), Some(positive(2)), Some(positive(5))]
        );
        let before = output;
        assert!(
            store
                .lookup_batch(&[key(0, 4), key(1, 3), key(1, 4)], &mut output)
                .is_err()
        );
        assert_eq!(output, before);
        assert!(store.lookup_batch(&[key(0, 4); 4], &mut [None; 4]).is_err());
        assert!(store.lookup_batch(&[key(0, 4)], &mut []).is_err());
        let view = CountView::new(&store, positive(3));
        assert_eq!(view.lookup(key(0, 4)).unwrap(), None);
        assert_eq!(view.lookup(key(1, 4)).unwrap(), Some(positive(5)));
        let mut visits = 0;
        assert!(
            view.visit(&mut |_, _| {
                visits += 1;
                bail!("stop")
            })
            .is_err()
        );
        assert_eq!(visits, 1);
        assert_eq!(store.metadata().occurrences, 7);
        assert!(store.visit(&mut |_, _| bail!("reader failure")).is_err());
    }
}
