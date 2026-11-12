//! Exact v4 count-store interface and in-memory conformance candidate.
//!
//! The slot budget covers owned table buffers, including both sides of a
//! resize. Callers must separately reserve process/allocator overhead, input
//! buffers, graphs and evidence before using this as a production backend.

use super::budget::{MemoryClass, MemoryReservation, ResourceBudget, ResourceLimits};
use super::identity::{CanonicalKey, KmerSize, OrientedKey};
use anyhow::{Context, Result, ensure};
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

type CountMap = hashbrown::HashMap<u64, u64, TableHasher>;

struct Table {
    k: KmerSize,
    map: CountMap,
    occurrences: u64,
    table_budget_bytes: usize,
    budget: Arc<ResourceBudget>,
    table_reservation: Option<MemoryReservation>,
    peak_table_bytes: usize,
    max_lookup_batch: usize,
}

impl Table {
    // hashbrown 0.16.1, (u64,u64): 16 bytes per bucket, one control byte
    // per bucket and at most 16 trailing control bytes. Normal growth starts
    // with 14 entries and doubles, giving power-of-two bucket arrays at 7/8
    // load. The pin and layout tests protect this pre-allocation bound.
    fn allocation_bound(capacity: usize) -> Result<usize> {
        ensure!(capacity >= 14, "Count table reservation is too small");
        capacity
            .checked_mul(8)
            .context("Count capacity overflow")?
            .checked_div(7)
            .context("Invalid count load factor")?
            .checked_next_power_of_two()
            .context("Count bucket overflow")?
            .checked_mul(17)
            .and_then(|bytes| bytes.checked_add(16))
            .filter(|bytes| *bytes <= isize::MAX as usize - 15)
            .context("Count table byte size overflow")
    }

    fn grow(&mut self) -> Result<()> {
        let capacity = if self.map.capacity() == 0 {
            14
        } else {
            self.map
                .capacity()
                .checked_mul(2)
                .context("Count capacity overflow")?
        };
        let bound = Self::allocation_bound(capacity)?;
        let old_bytes = self.map.allocation_size();
        let peak = old_bytes
            .checked_add(bound)
            .context("Count resize byte size overflow")?;
        ensure!(
            peak <= self.table_budget_bytes,
            "Count table resize needs {peak} table bytes; budget is {}",
            self.table_budget_bytes
        );
        let mut reservation = self
            .budget
            .reserve(MemoryClass::CountSlots, u64::try_from(bound)?)?;
        let mut replacement = CountMap::with_hasher(self.map.hasher().clone());
        replacement
            .try_reserve(capacity)
            .map_err(|error| anyhow::anyhow!("Unable to reserve exact count table: {error:?}"))?;
        let allocated = replacement.allocation_size();
        ensure!(
            allocated <= bound && replacement.capacity() >= capacity,
            "Count table allocation differs from its pinned reservation bound"
        );
        reservation.resize(u64::try_from(allocated)?)?;
        // Every fallible allocation/check precedes mutation of accepted counts.
        // This pre-sized destination has room for all old entries plus the new
        // key. Built-in integer hashing invokes no user code and cannot fail.
        for (key, count) in self.map.drain() {
            replacement.insert(key, count);
        }
        debug_assert_eq!(replacement.allocation_size(), allocated);
        self.peak_table_bytes = self.peak_table_bytes.max(old_bytes + allocated);
        self.map = replacement;
        self.table_reservation = Some(reservation);
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
    pub fn new(k: KmerSize, table_budget_bytes: usize, max_lookup_batch: usize) -> Result<Self> {
        let bytes = u64::try_from(table_budget_bytes)?;
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
        let table_budget_bytes = usize::try_from(budget.class_limit(MemoryClass::CountSlots))?;
        let mut table = Table {
            k,
            map: CountMap::with_hasher(TableHasher::default()),
            occurrences: 0,
            table_budget_bytes,
            budget,
            table_reservation: None,
            peak_table_bytes: 0,
            max_lookup_batch,
        };
        table.grow()?;
        Ok(Self { table })
    }

    /// Add a positive exact count. Validation/overflow/capacity errors preserve
    /// previously accepted counts and totals; the caller must fail publication.
    pub fn add(&mut self, key: CanonicalKey, increment: NonZeroU64) -> Result<()> {
        self.table.check_key(key)?;
        let total = self
            .table
            .occurrences
            .checked_add(increment.get())
            .context("Exact k-mer occurrence total overflow")?;
        // At capacity, an existing count must remain updatable without a new
        // allocation. Otherwise reserve before entry() could grow implicitly.
        if self.table.map.len() == self.table.map.capacity() {
            if let Some(count) = self.table.map.get_mut(&key.bits()) {
                *count = count
                    .checked_add(increment.get())
                    .context("Exact k-mer count overflow")?;
                self.table.occurrences = total;
                return Ok(());
            }
            self.table.grow()?;
        }
        match self.table.map.entry(key.bits()) {
            hashbrown::hash_map::Entry::Occupied(mut entry) => {
                let count = entry
                    .get()
                    .checked_add(increment.get())
                    .context("Exact k-mer count overflow")?;
                *entry.get_mut() = count;
            }
            hashbrown::hash_map::Entry::Vacant(entry) => {
                entry.insert(increment.get());
            }
        }
        self.table.occurrences = total;
        Ok(())
    }

    /// Bind the counted input digest and move into a shared immutable store.
    pub fn finalize(self, selected_input_sha256: [u8; 32]) -> MemoryCountStore {
        let metadata = CountMetadata {
            format_version: 1,
            k: self.table.k,
            count_width_bits: 64,
            backend: "hashbrown_memory",
            hash_algorithm: TABLE_HASHER_NAME,
            canonicalization: "min_forward_reverse_complement",
            preprocessing: "uppercase_acgtn_n_breaks_windows",
            overflow_policy: "error_preserve_previous_counts",
            selected_input_sha256,
            distinct_keys: self.table.map.len() as u64,
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
    /// Current entry capacity, excluding control bytes and vacant-bucket slack.
    pub fn capacity(&self) -> usize {
        self.table.map.capacity()
    }

    /// Owned table allocation bytes, including control bytes and bucket slack.
    pub fn table_bytes(&self) -> usize {
        self.table.map.allocation_size()
    }

    /// Peak owned slot bytes, including simultaneous old/new resize buffers.
    /// This excludes allocator/process overhead and is not an RSS measurement.
    pub fn peak_table_bytes(&self) -> usize {
        self.table.peak_table_bytes
    }
}

impl CountStore for MemoryCountStore {
    fn metadata(&self) -> &CountMetadata {
        &self.metadata
    }

    fn lookup(&self, key: CanonicalKey) -> Result<Option<NonZeroU64>> {
        self.table.check_key(key)?;
        Ok(self
            .table
            .map
            .get(&key.bits())
            .copied()
            .and_then(NonZeroU64::new))
    }

    fn lookup_oriented(&self, key: OrientedKey) -> Result<Option<NonZeroU64>> {
        // A single canonical probe avoids an extra random miss for reverse
        // orientations; typed reverse complementation uses constant work.
        self.lookup(key.canonical())
    }

    fn visit(&self, visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>) -> Result<()> {
        self.visit_entries(visitor)
    }
}

impl MemoryCountStore {
    /// Visit a concrete in-memory store with an inlineable visitor. The shared
    /// fallible trait boundary remains available for opaque and disk backends.
    pub fn visit_entries<F: FnMut(CanonicalKey, NonZeroU64) -> Result<()> + ?Sized>(
        &self,
        visitor: &mut F,
    ) -> Result<()> {
        for (&bits, &count) in &self.table.map {
            let count = NonZeroU64::new(count).context("Exact count table contains zero")?;
            visitor(CanonicalKey::from_stored(bits, self.metadata.k), count)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::identity::KmerWindows;
    use anyhow::bail;
    use std::collections::BTreeMap;

    fn positive(n: u64) -> NonZeroU64 {
        NonZeroU64::new(n).unwrap()
    }
    fn key(bits: u64, k: usize) -> CanonicalKey {
        CanonicalKey::new(bits, KmerSize::new(k).unwrap()).unwrap()
    }

    #[test]
    fn pinned_allocation_bounds_cover_control_bytes_and_all_growth_sizes() {
        let mut map = CountMap::with_hasher(TableHasher::default());
        for capacity in [14, 28, 56, 112, 224, 448, 896, 1792, 3584, 7168] {
            map.try_reserve(capacity).unwrap();
            assert!(map.capacity() >= capacity);
            assert!(map.allocation_size() <= Table::allocation_bound(capacity).unwrap());
            assert!(map.allocation_size() >= map.capacity() * 16);
            map.clear();
        }
        assert!(Table::allocation_bound(usize::MAX).is_err());
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
        let initial_bytes = builder.table.map.allocation_size() as u64;
        let next_bound = Table::allocation_bound(builder.table.map.capacity() * 2).unwrap() as u64;
        let graph = budget.reserve(MemoryClass::TargetGraphs, 500).unwrap();
        let mut keys = std::collections::BTreeSet::new();
        for bits in 0..4096 {
            keys.insert(key(bits, 6));
            if keys.len() == 15 {
                break;
            }
        }
        let keys = keys.into_iter().collect::<Vec<_>>();
        for key in &keys[..14] {
            builder.add(*key, positive(1)).unwrap();
        }
        assert!(builder.add(keys[14], positive(1)).is_err());
        assert_eq!(builder.table.occurrences, 14);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 600 + initial_bytes);
        drop(graph);
        builder.add(keys[14], positive(1)).unwrap();
        let store = builder.finalize([0; 32]);
        assert_eq!(store.lookup(keys[14]).unwrap(), Some(positive(1)));
        assert_eq!(store.metadata().occurrences, 15);
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            100 + store.table_bytes() as u64
        );
        assert_eq!(
            budget.snapshot().peak_reserved_memory_bytes,
            100 + initial_bytes + next_bound
        );
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
        assert_eq!(store.metadata().backend, "hashbrown_memory");
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
        assert!(store.peak_table_bytes() <= 1 << 20);
    }

    #[test]
    fn hash_collisions_never_create_false_membership_and_zero_key_is_valid() {
        let mut builder = MemoryCountBuilder::new(KmerSize::new(8).unwrap(), 4096, 8).unwrap();
        // Same initial bucket and seven-bit control fingerprint: full key
        // comparison must still distinguish every inserted and absent key.
        let home = |bits| {
            let hash = std::hash::BuildHasher::hash_one(builder.table.map.hasher(), bits);
            (hash & 15, hash >> 57)
        };
        let target = home(0u64);
        let mut colliders = Vec::new();
        for bits in 0..65536 {
            let candidate = key(bits, 8);
            if !colliders.contains(&candidate) && home(candidate.bits()) == target {
                colliders.push(candidate);
            }
            if colliders.len() == 4 {
                break;
            }
        }
        assert_eq!(colliders.len(), 4);
        assert_eq!(colliders[0], key(0, 8));
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
        // The next 560-byte table fits alone; the old+new allocation does not.
        let mut builder = MemoryCountBuilder::new(KmerSize::new(6).unwrap(), 600, 8).unwrap();
        let initial_bytes = builder.table.map.allocation_size();
        let mut accepted = Vec::new();
        for bits in 0..4096 {
            let k = key(bits, 6);
            if accepted.contains(&k) {
                continue;
            }
            if accepted.len() == 14 {
                assert!(builder.add(k, positive(1)).is_err());
                break;
            }
            builder.add(k, positive(1)).unwrap();
            accepted.push(k);
        }
        // A full table can update an existing key without reserving growth.
        builder.add(accepted[0], positive(1)).unwrap();
        let pointer = builder.table.map.get(&accepted[0].bits()).unwrap() as *const u64;
        let store = builder.finalize([0; 32]);
        assert_eq!(
            store.table.map.get(&accepted[0].bits()).unwrap() as *const u64,
            pointer
        );
        assert_eq!(store.metadata().distinct_keys, 14);
        assert_eq!(store.metadata().occurrences, 15);
        for (i, k) in accepted.into_iter().enumerate() {
            assert_eq!(
                store.lookup(k).unwrap(),
                Some(positive(if i == 0 { 2 } else { 1 }))
            );
        }
        assert_eq!(store.peak_table_bytes(), initial_bytes);
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
