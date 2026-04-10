// kmer/histogram.rs
//! Histogram of kmer count frequencies.

use ahash::AHashMap;
use anyhow::{Context, Result};

use super::counting::KmerCounts;

/// A structure to hold a histogram of kmer counts.
/// For each count, the histogram stores the number of kmers (n_kmers) with that count.
#[derive(Debug, Clone)]
pub struct Histogram {
    histo: Vec<u64>, // for 0 through histo_max, the number of kmers with that count. 0 isn't actually expected to have a value.
    histo_large: AHashMap<u64, u64>, // for counts larger than histo_max, the number of kmers with that count
    histo_max: u64,
}

impl Histogram {
    pub fn new(histo_max: &u64) -> Histogram {
        let length = *histo_max as usize + 2;
        let histo: Vec<u64> = vec![0; length];
        let histo_large: AHashMap<u64, u64> = AHashMap::default();
        Histogram {
            histo,
            histo_large,
            histo_max: *histo_max,
        }
    }

    #[allow(dead_code)]
    pub fn ingest_kmer_counts(&mut self, kmer_counts: &KmerCounts) {
        for (_, count) in kmer_counts.iter() {
            let count = *count as u64;
            if count <= self.histo_max {
                self.histo[count as usize] += 1;
            } else {
                let counter = self.histo_large.entry(count).or_insert(0);
                *counter += 1;
            }
        }
    }

    /// Move a kmer from one count bin to another.
    /// Decrements the bin for `old_count` and increments the bin for `new_count`.
    ///
    /// Callers must pass a valid `old_count` that actually has at least one
    /// kmer in the corresponding bin, or `old_count == new_count` (no-op).
    /// This is asserted in debug builds; release builds tolerate the
    /// inconsistency via saturating subtraction rather than panicking,
    /// because a histogram miscount is strictly less harmful than a crash.
    pub(crate) fn move_count(&mut self, old_count: u64, new_count: u64) {
        if old_count == new_count {
            return;
        }
        // Decrement old bin
        if old_count > 0 {
            if old_count <= self.histo_max {
                let idx = old_count as usize;
                debug_assert!(
                    self.histo[idx] > 0,
                    "Histogram::move_count: bin for old_count={} is empty; caller invariant violated",
                    old_count
                );
                self.histo[idx] = self.histo[idx].saturating_sub(1);
            } else if let Some(n) = self.histo_large.get_mut(&old_count) {
                debug_assert!(*n > 0, "Histogram::move_count: large bin is zero");
                *n = n.saturating_sub(1);
                if *n == 0 {
                    self.histo_large.remove(&old_count);
                }
            } else {
                debug_assert!(
                    false,
                    "Histogram::move_count: no large bin for old_count={}",
                    old_count
                );
            }
        }
        // Increment new bin
        if new_count <= self.histo_max {
            self.histo[new_count as usize] += 1;
        } else {
            *self.histo_large.entry(new_count).or_insert(0) += 1;
        }
    }

    #[allow(dead_code)]
    pub fn from_kmer_counts(kmer_counts: &KmerCounts, histo_max: &u64) -> Histogram {
        let mut histo = Histogram::new(histo_max);
        histo.ingest_kmer_counts(kmer_counts);
        histo
    }

    #[allow(dead_code)]
    pub fn get(&self, count: &u64) -> u64 {
        if *count <= self.histo_max {
            self.histo[*count as usize]
        } else {
            *self.histo_large.get(count).unwrap_or(&0)
        }
    }

    pub fn get_n_kmers(&self) -> u64 {
        let small: u64 = self
            .histo
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &n)| n * i as u64)
            .sum();
        let large: u64 = self
            .histo_large
            .iter()
            .map(|(count, n_kmers)| count * n_kmers)
            .sum();
        small + large
    }

    pub fn get_n_unique_kmers(&self) -> u64 {
        let small: u64 = self.histo.iter().skip(1).sum();
        let large: u64 = self.histo_large.values().sum();
        small + large
    }

    pub fn get_vector(&self) -> Result<Vec<u64>> {
        let mut histo_vec = self.histo.clone();
        let last = histo_vec.last_mut().context("Histogram vector is empty")?;

        for n_kmers in self.histo_large.values() {
            *last += n_kmers;
        }

        Ok(histo_vec)
    }
}
