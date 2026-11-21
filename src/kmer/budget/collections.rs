//! Owned target-local vector and hash-table allocations shared with assembly.
//! Capacity and resize overlap are charged before allocation. Element-owned
//! heap data and allocator/process overhead require separate ownership.

use super::{MemoryClass, MemoryReservation, ResourceBudget};
use anyhow::{Context, Result, ensure};
use std::ops::Deref;
use std::sync::Arc;

/// Backing-allocation guard retained after its container data drops.
pub(crate) struct Charge {
    reservation: Option<MemoryReservation>,
    budget: Option<Arc<ResourceBudget>>,
    class: MemoryClass,
}

impl Charge {
    /// Borrow the selected ledger without exposing accounting mutation.
    pub(crate) fn budget(&self) -> Option<&Arc<ResourceBudget>> {
        self.budget.as_ref()
    }

    /// Report the class inherited by copies and nested buffers.
    pub(crate) fn class(&self) -> MemoryClass {
        self.class
    }

    /// Report whether this owner has an active target ceiling.
    pub(crate) fn is_enabled(&self) -> bool {
        self.budget.is_some()
    }

    /// Start with no allocation and borrow the optional target ceiling.
    pub(crate) fn new(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self::in_class(budget, MemoryClass::TargetGraphs)
    }

    /// Start ownership in the named stage without reserving nonexistent storage.
    pub(crate) fn in_class(budget: Option<&Arc<ResourceBudget>>, class: MemoryClass) -> Self {
        Self {
            reservation: None,
            budget: budget.cloned(),
            class,
        }
    }

    /// Charge a replacement while the current allocation remains live.
    pub(crate) fn reserve(&self, bytes: u64) -> Result<Option<MemoryReservation>> {
        self.budget
            .as_ref()
            .map(|budget| budget.reserve(self.class, bytes))
            .transpose()
    }

    /// Attach the surviving guard before checking the observed allocation.
    pub(crate) fn replace(
        &mut self,
        reservation: Option<MemoryReservation>,
        actual: u64,
        bound: u64,
    ) -> Result<()> {
        self.reservation = reservation;
        if let Some(reservation) = self.reservation.as_mut() {
            ensure!(
                actual <= bound,
                "Workspace allocation exceeded its precharged bound"
            );
            reservation.resize(actual)?;
        }
        Ok(())
    }
}

/// Compute logical backing bytes with checked allocation-size arithmetic.
pub(crate) fn bytes<T>(capacity: usize) -> Result<u64> {
    let bytes = capacity
        .checked_mul(size_of::<T>())
        .filter(|n| *n <= isize::MAX as usize)
        .context("Workspace buffer size overflow")?;
    u64::try_from(bytes).context("Workspace buffer exceeds supported range")
}

/// A vector whose backing allocation drops before its ownership guard.
pub(crate) struct BudgetVec<T> {
    data: Vec<T>,
    charge: Charge,
}

impl<T> BudgetVec<T> {
    /// Start empty, sharing an optional target ceiling.
    pub(crate) fn new(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self::in_class(budget, MemoryClass::TargetGraphs)
    }

    /// Start an empty vector charged to a caller-selected stage.
    pub(crate) fn in_class(budget: Option<&Arc<ResourceBudget>>, class: MemoryClass) -> Self {
        Self {
            data: Vec::new(),
            charge: Charge::in_class(budget, class),
        }
    }

    /// Copy nested owners fallibly, retaining class and resize accounting.
    pub(crate) fn try_clone_with(&self, mut copy: impl FnMut(&T) -> Result<T>) -> Result<Self> {
        let mut result = Self::in_class(self.budget(), self.class());
        result.reserve_capacity(self.len())?;
        for value in self.iter() {
            result.push(copy(value)?)?;
        }
        Ok(result)
    }

    /// Collect values under the selected stage, charging before allocation.
    pub(crate) fn collect_in_class(
        iter: impl IntoIterator<Item = T>,
        budget: Option<&Arc<ResourceBudget>>,
        class: MemoryClass,
    ) -> Result<Self> {
        let mut result = Self::in_class(budget, class);
        let iter = iter.into_iter();
        result.reserve_capacity(iter.size_hint().0)?;
        for item in iter {
            result.push(item)?;
        }
        Ok(result)
    }

    /// Precharge exact indexed output before any parallel callback can run.
    /// Rayon fills the already allocated vector without exposing growing access.
    pub(crate) fn collect_indexed<I>(
        iter: I,
        budget: Option<&Arc<ResourceBudget>>,
        class: MemoryClass,
    ) -> Result<Self>
    where
        T: Send,
        I: rayon::iter::IndexedParallelIterator<Item = T>,
    {
        let length = iter.len();
        let mut result = Self::in_class(budget, class);
        result.reserve_capacity(length)?;
        iter.collect_into_vec(&mut result.data);
        ensure!(
            result.data.len() == length,
            "Indexed collection changed output length"
        );
        ensure!(
            result.data.capacity() == length || size_of::<T>() == 0,
            "Indexed collection changed its reserved capacity"
        );
        Ok(result)
    }

    /// Keep the stage identity when deriving nested or copied owners.
    pub(crate) fn class(&self) -> MemoryClass {
        self.charge.class
    }

    /// Drop trailing values without discarding reusable backing capacity.
    pub(crate) fn truncate(&mut self, length: usize) {
        self.data.truncate(length);
    }

    /// Test emptiness without exposing the underlying growing vector.
    pub(crate) fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Reserve known backing capacity before creating any values.
    pub(crate) fn with_capacity(
        capacity: usize,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<Self> {
        Self::with_capacity_in_class(capacity, budget, MemoryClass::TargetGraphs)
    }

    /// Admit exact backing capacity in the caller-selected memory class.
    pub(crate) fn with_capacity_in_class(
        capacity: usize,
        budget: Option<&Arc<ResourceBudget>>,
        class: MemoryClass,
    ) -> Result<Self> {
        let mut result = Self::in_class(budget, class);
        result.reserve_capacity(capacity)?;
        Ok(result)
    }

    /// Collect with fallible growth and retained ownership, without a raw copy.
    pub(crate) fn collect(
        iter: impl IntoIterator<Item = T>,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<Self> {
        let mut result = Self::new(budget);
        let iter = iter.into_iter();
        let lower = iter.size_hint().0;
        if lower > 0 {
            result.reserve_capacity(lower)?;
        }
        for item in iter {
            result.push(item)?;
        }
        Ok(result)
    }

    fn reserve_capacity(&mut self, capacity: usize) -> Result<()> {
        if capacity <= self.data.capacity() {
            return Ok(());
        }
        let reservation = self.charge.reserve(bytes::<T>(capacity)?)?;
        self.data
            .try_reserve_exact(
                capacity
                    .checked_sub(self.data.len())
                    .context("Invalid workspace capacity")?,
            )
            .context("Unable to allocate workspace vector")?;
        self.charge.replace(
            reservation,
            bytes::<T>(self.data.capacity())?,
            bytes::<T>(capacity)?,
        )
    }

    /// Reserve before push can grow; failed growth preserves existing entries.
    /// Admit the next element before an external operation creates its payload.
    pub(crate) fn reserve_one(&mut self) -> Result<()> {
        let needed = self
            .data
            .len()
            .checked_add(1)
            .context("Vector length overflow")?;
        self.reserve_capacity(needed)
    }

    pub(crate) fn push(&mut self, item: T) -> Result<()> {
        if self.data.len() == self.data.capacity() {
            self.reserve_capacity(
                self.data
                    .capacity()
                    .checked_mul(2)
                    .context("Workspace capacity overflow")?
                    .max(4),
            )?;
        }
        self.data.push(item);
        Ok(())
    }

    /// Append copied elements after admitting old/new backing overlap.
    pub(crate) fn extend_from_slice(&mut self, items: &[T]) -> Result<()>
    where
        T: Copy,
    {
        let needed = self
            .data
            .len()
            .checked_add(items.len())
            .context("Vector length overflow")?;
        if needed > self.data.capacity() {
            let capacity = self
                .data
                .capacity()
                .checked_mul(2)
                .context("Workspace capacity overflow")?
                .max(needed);
            self.reserve_capacity(capacity)?;
        }
        self.data.extend_from_slice(items);
        Ok(())
    }

    /// Report backing capacity without exposing a growing mutable allocation.
    pub(crate) fn capacity(&self) -> usize {
        self.data.capacity()
    }

    /// Borrow the ceiling for nested or replacement storage.
    pub(crate) fn budget(&self) -> Option<&Arc<ResourceBudget>> {
        self.charge.budget.as_ref()
    }

    /// Remove a value while retaining and charging reusable capacity.
    pub(crate) fn pop(&mut self) -> Option<T> {
        self.data.pop()
    }

    /// Clear values while retaining and charging reusable capacity.
    pub(crate) fn clear(&mut self) {
        self.data.clear();
    }
}

impl<T> Deref for BudgetVec<T> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        &self.data
    }
}

impl<T> std::ops::DerefMut for BudgetVec<T> {
    fn deref_mut(&mut self) -> &mut [T] {
        &mut self.data
    }
}

/// Consuming iteration keeps the backing reservation until its storage drops.
pub(crate) struct BudgetIntoIter<T> {
    iter: std::vec::IntoIter<T>,
    _charge: Charge,
}

impl<T> Iterator for BudgetIntoIter<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        self.iter.next()
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<T> IntoIterator for BudgetVec<T> {
    type Item = T;
    type IntoIter = BudgetIntoIter<T>;
    fn into_iter(self) -> Self::IntoIter {
        BudgetIntoIter {
            iter: self.data.into_iter(),
            _charge: self.charge,
        }
    }
}

impl<'a, T> IntoIterator for &'a BudgetVec<T> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

impl<T> AsRef<[T]> for BudgetVec<T> {
    fn as_ref(&self) -> &[T] {
        &self.data
    }
}

impl<T: Copy> BudgetVec<T> {
    /// Copy values into an independently reserved buffer under the same ceiling.
    pub(crate) fn try_clone(&self) -> Result<Self> {
        let mut result = Self::in_class(self.budget(), self.class());
        result.reserve_capacity(self.data.len())?;
        result.data.extend_from_slice(&self.data);
        Ok(result)
    }
}

impl<T: std::fmt::Debug> std::fmt::Debug for BudgetVec<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}

impl<T: PartialEq> PartialEq for BudgetVec<T> {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl<T> Default for BudgetVec<T> {
    fn default() -> Self {
        Self::in_class(None, MemoryClass::RetainedResults)
    }
}
#[cfg(test)]
impl<T: PartialEq> PartialEq<Vec<T>> for BudgetVec<T> {
    fn eq(&self, other: &Vec<T>) -> bool {
        self.data == *other
    }
}
impl<T: Eq> Eq for BudgetVec<T> {}
impl<T: std::hash::Hash> std::hash::Hash for BudgetVec<T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.data.hash(state);
    }
}
impl<T> std::borrow::Borrow<[T]> for BudgetVec<T> {
    fn borrow(&self) -> &[T] {
        &self.data
    }
}
#[cfg(test)]
impl<T: Clone> Clone for BudgetVec<T> {
    fn clone(&self) -> Self {
        self.try_clone_with(|value| Ok(value.clone())).unwrap()
    }
}
#[cfg(test)]
impl<T> From<Vec<T>> for BudgetVec<T> {
    fn from(values: Vec<T>) -> Self {
        Self::collect_in_class(values, None, MemoryClass::RetainedResults).unwrap()
    }
}

/// Fallible deep copying preserves resource ownership and typed refusal.
pub(crate) trait TryClone: Sized {
    fn try_clone(&self) -> Result<Self>;
}
macro_rules! copy_values {($($value:ty),*)=>{$(impl TryClone for $value {fn try_clone(&self)->Result<Self> {Ok(*self)}})*};}
copy_values!(usize, u64, u32, bool, &'static str);
impl<T: TryClone> TryClone for Option<T> {
    fn try_clone(&self) -> Result<Self> {
        self.as_ref().map(TryClone::try_clone).transpose()
    }
}
impl<T: TryClone> TryClone for BudgetVec<T> {
    fn try_clone(&self) -> Result<Self> {
        self.try_clone_with(TryClone::try_clone)
    }
}

/// Bound the pinned hashbrown 0.16.1 bucket/control layout before allocation.
pub(crate) fn hash_allocation_bound<T>(capacity: usize) -> Result<u64> {
    ensure!(capacity >= 14, "Owned table reservation is too small");
    let buckets = capacity
        .checked_mul(8)
        .and_then(|n| n.checked_div(7))
        .and_then(usize::checked_next_power_of_two)
        .context("Owned table bucket overflow")?;
    let alignment = std::mem::align_of::<T>().max(16);
    let data = size_of::<T>()
        .checked_mul(buckets)
        .and_then(|n| n.checked_add(alignment - 1))
        .map(|n| n & !(alignment - 1))
        .context("Owned table data size overflow")?;
    let total = data
        .checked_add(buckets)
        .and_then(|n| n.checked_add(16))
        .filter(|n| *n <= isize::MAX as usize - (alignment - 1))
        .context("Owned table allocation overflow")?;
    u64::try_from(total).context("Owned table allocation exceeds supported range")
}

/// A local exact map with fallible growth and no growing mutable escape.
pub(crate) struct BudgetMap<K, V> {
    data: hashbrown::HashMap<K, V, ahash::RandomState>,
    charge: Charge,
}

impl<K: Eq + std::hash::Hash, V> BudgetMap<K, V> {
    /// Start empty with optional ownership under the target ceiling.
    pub(crate) fn new(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self::in_class(budget, MemoryClass::TargetGraphs)
    }

    /// Start a map whose backing allocation belongs to the selected stage.
    pub(crate) fn in_class(budget: Option<&Arc<ResourceBudget>>, class: MemoryClass) -> Self {
        Self {
            data: hashbrown::HashMap::with_hasher(ahash::RandomState::new()),
            charge: Charge::in_class(budget, class),
        }
    }

    /// Borrow the ceiling for temporary buffers derived from this map.
    pub(crate) fn budget(&self) -> Option<&Arc<ResourceBudget>> {
        self.charge.budget.as_ref()
    }

    /// Reserve old/replacement allocations before inserting a new key.
    pub(crate) fn insert(&mut self, key: K, value: V) -> Result<()> {
        if self.charge.is_enabled()
            && self.data.len() == self.data.capacity()
            && !self.data.contains_key(&key)
        {
            let capacity = self
                .data
                .capacity()
                .checked_mul(2)
                .context("Owned table capacity overflow")?
                .max(14);
            let bound = hash_allocation_bound::<(K, V)>(capacity)?;
            let reservation = self.charge.reserve(bound)?;
            self.data
                .try_reserve(
                    capacity
                        .checked_sub(self.data.len())
                        .context("Invalid owned table capacity")?,
                )
                .map_err(|error| anyhow::anyhow!("Unable to allocate local table: {error:?}"))?;
            self.charge.replace(
                reservation,
                u64::try_from(self.data.allocation_size())?,
                bound,
            )?;
            ensure!(
                self.data.capacity() >= capacity,
                "Owned table capacity differs from its pinned bound"
            );
        }
        self.data.insert(key, value);
        Ok(())
    }

    /// Borrow a value without canonicalizing or allocating.
    pub(crate) fn get<Q: ?Sized + Eq + std::hash::Hash>(&self, key: &Q) -> Option<&V>
    where
        K: std::borrow::Borrow<Q>,
    {
        self.data.get(key)
    }

    /// Mutate an existing value without exposing a growing table operation.
    pub(crate) fn get_mut<Q: ?Sized + Eq + std::hash::Hash>(&mut self, key: &Q) -> Option<&mut V>
    where
        K: std::borrow::Borrow<Q>,
    {
        self.data.get_mut(key)
    }

    /// Borrow entries in unspecified order.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        self.data.iter()
    }

    /// Borrow values without exposing a growing table operation.
    pub(crate) fn values(&self) -> impl Iterator<Item = &V> {
        self.data.values()
    }

    /// Report local cardinality.
    pub(crate) fn len(&self) -> usize {
        self.data.len()
    }

    /// Report whether the local table is empty.
    pub(crate) fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

/// Consuming table iteration retains its backing charge through destruction.
pub(crate) struct BudgetMapIntoIter<K, V> {
    iter: hashbrown::hash_map::IntoIter<K, V>,
    _charge: Charge,
}

impl<K, V> Iterator for BudgetMapIntoIter<K, V> {
    type Item = (K, V);
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<K, V> IntoIterator for BudgetMap<K, V> {
    type Item = (K, V);
    type IntoIter = BudgetMapIntoIter<K, V>;
    fn into_iter(self) -> Self::IntoIter {
        BudgetMapIntoIter {
            iter: self.data.into_iter(),
            _charge: self.charge,
        }
    }
}

#[cfg(test)]
impl<K: Eq + std::hash::Hash + std::fmt::Debug, V: std::fmt::Debug> std::fmt::Debug
    for BudgetMap<K, V>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}

#[cfg(test)]
impl<K: Eq + std::hash::Hash, V: PartialEq> PartialEq for BudgetMap<K, V> {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl<T: serde::Serialize> serde::Serialize for BudgetVec<T> {
    fn serialize<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error> {
        serde::Serialize::serialize(self.data.as_slice(), serializer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iterator_budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(super::super::ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn copied_result_vectors_and_consuming_iteration_preserve_the_result_class() {
        let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
        classes[MemoryClass::RetainedResults as usize] = 1024;
        let budget = ResourceBudget::new(super::super::ResourceLimits {
            memory_bytes: 1024,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap();
        let mut values = BudgetVec::in_class(Some(&budget), MemoryClass::RetainedResults);
        for value in 0..4u64 {
            values.push(value).unwrap();
        }
        let copy = values.try_clone().unwrap();
        assert_eq!(copy.class(), MemoryClass::RetainedResults);
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::RetainedResults as usize],
            64
        );
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
        let mut iter = values.into_iter();
        while iter.next().is_some() {}
        assert_eq!(budget.snapshot().reserved_memory_bytes, 64);
        drop(iter);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 32);
        drop(copy);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn consuming_map_iteration_retains_table_charge_and_nested_values_independently() {
        let budget = iterator_budget(8192);
        let mut map = BudgetMap::new(Some(&budget));
        map.insert(1u64, BudgetVec::collect(0..4u64, Some(&budget)).unwrap())
            .unwrap();
        map.insert(2u64, BudgetVec::collect(4..8u64, Some(&budget)).unwrap())
            .unwrap();
        let original = budget.snapshot().reserved_memory_bytes;
        let mut iter = map.into_iter();
        let (_, first) = iter.next().unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, original);
        drop(iter);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 32);
        drop(first);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        let mut map = BudgetMap::new(Some(&budget));
        map.insert(1u64, 2u64).unwrap();
        let live = budget.snapshot().reserved_memory_bytes;
        let mut iter = map.into_iter();
        assert_eq!(iter.next(), Some((1, 2)));
        assert_eq!(iter.next(), None);
        assert_eq!(budget.snapshot().reserved_memory_bytes, live);
        drop(iter);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn fallible_copy_and_consuming_iteration_keep_backing_charges_until_storage_drops() {
        let budget = iterator_budget(200);
        let values = BudgetVec::collect(0..4u64, Some(&budget)).unwrap();
        let held = budget.reserve(MemoryClass::TargetGraphs, 140).unwrap();
        assert!(
            values
                .try_clone()
                .unwrap_err()
                .downcast_ref::<super::super::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(&*values, &[0, 1, 2, 3]);
        drop(held);
        let copy = values.try_clone().unwrap();
        let mut iter = values.into_iter();
        for value in 0..4 {
            assert_eq!(iter.next(), Some(value));
        }
        assert_eq!(iter.next(), None);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 64);
        drop(iter);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 32);
        drop(copy);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn nested_buffer_iteration_releases_outer_and_unconsumed_owners_independently() {
        let budget = iterator_budget(1024);
        let mut paths = BudgetVec::new(Some(&budget));
        paths
            .push(BudgetVec::collect(0..4u64, Some(&budget)).unwrap())
            .unwrap();
        paths
            .push(BudgetVec::collect(4..8u64, Some(&budget)).unwrap())
            .unwrap();
        let mut iter = paths.into_iter();
        let first = iter.next().unwrap();
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            2 * 32 + 4 * size_of::<BudgetVec<u64>>() as u64
        );
        drop(iter);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 32);
        assert_eq!(&*first, &[0, 1, 2, 3]);
        drop(first);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn local_map_layout_and_growth_keep_actual_allocation_owned() {
        let budget = ResourceBudget::new(super::super::ResourceLimits {
            memory_bytes: 1 << 20,
            process_overhead_bytes: 0,
            class_bytes: [1 << 20; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap();
        let mut map = BudgetMap::new(Some(&budget));
        for key in 0..4096u64 {
            map.insert(key, key + 1).unwrap();
            assert_eq!(map.get(&key), Some(&(key + 1)));
            let actual = map.data.allocation_size() as u64;
            assert_eq!(budget.snapshot().reserved_memory_bytes, actual);
            assert!(actual <= hash_allocation_bound::<(u64, u64)>(map.data.capacity()).unwrap());
        }
        let before = budget.snapshot().reserved_memory_bytes;
        map.insert(0, 100).unwrap();
        assert_eq!(map.get(&0), Some(&100));
        assert_eq!(budget.snapshot().reserved_memory_bytes, before);
        drop(map);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(hash_allocation_bound::<(u64, u64)>(usize::MAX).is_err());
    }
    #[test]
    fn indexed_parallel_collection_reserves_before_callbacks_and_holds_backing_through_drop() {
        use crate::kmer::budget::text::BudgetString;
        use rayon::prelude::*;
        use std::sync::atomic::{AtomicUsize, Ordering};
        struct Observed {
            index: usize,
            text: BudgetString,
            budget: Arc<ResourceBudget>,
            errors: Arc<AtomicUsize>,
            backing: u64,
        }
        impl Drop for Observed {
            fn drop(&mut self) {
                if self.budget.snapshot().reserved_memory_bytes < self.backing {
                    self.errors.fetch_add(1, Ordering::SeqCst);
                }
            }
        }
        let count = 64;
        let bytes = (count * size_of::<Observed>()) as u64;
        let callbacks = AtomicUsize::new(0);
        let mut classes = [0; super::super::MEMORY_CLASSES];
        classes[MemoryClass::RetainedResults as usize] = bytes - 1;
        let tiny = ResourceBudget::new(super::super::ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap();
        let refused = BudgetVec::collect_indexed(
            (0..count).into_par_iter().map(|_| {
                callbacks.fetch_add(1, Ordering::SeqCst);
                [0u8; 128]
            }),
            Some(&tiny),
            MemoryClass::RetainedResults,
        )
        .err()
        .unwrap();
        assert!(
            refused
                .downcast_ref::<super::super::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(callbacks.load(Ordering::SeqCst), 0);
        assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
        let budget = iterator_budget(1 << 20);
        let errors = Arc::new(AtomicUsize::new(0));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let values = pool
            .install(|| {
                BudgetVec::collect_indexed(
                    (0..count).into_par_iter().map(|index| {
                        assert!(
                            budget.snapshot().reserved_by_class
                                [MemoryClass::RetainedResults as usize]
                                >= bytes
                        );
                        Observed {
                            index,
                            text: BudgetString::format(
                                format_args!("row-{index}"),
                                Some(&budget),
                                MemoryClass::RetainedResults,
                            )
                            .unwrap(),
                            budget: Arc::clone(&budget),
                            errors: Arc::clone(&errors),
                            backing: bytes,
                        }
                    }),
                    Some(&budget),
                    MemoryClass::RetainedResults,
                )
            })
            .unwrap();
        assert_eq!(
            values.iter().map(|value| value.index).collect::<Vec<_>>(),
            (0..count).collect::<Vec<_>>()
        );
        assert_eq!(values[17].text, "row-17");
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
        drop(values);
        assert_eq!(errors.load(Ordering::SeqCst), 0);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn indexed_parallel_panic_releases_partial_elements_and_backing_reservation() {
        use crate::kmer::budget::text::BudgetString;
        use rayon::prelude::*;
        let budget = iterator_budget(1 << 20);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pool.install(|| {
                BudgetVec::collect_indexed(
                    (0..64usize).into_par_iter().map(|index| {
                        if index == 17 {
                            panic!("intentional indexed callback failure");
                        }
                        BudgetString::format(
                            format_args!("row-{index}"),
                            Some(&budget),
                            MemoryClass::RetainedResults,
                        )
                        .unwrap()
                    }),
                    Some(&budget),
                    MemoryClass::RetainedResults,
                )
                .unwrap()
            })
        }));
        assert!(outcome.is_err());
        assert!(budget.snapshot().valid);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }
}
