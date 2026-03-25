// kmer/histogram.rs
//! Histogram of kmer count frequencies.

use anyhow::{Context, Result};

use super::counting::KmerCounts;

/// A structure to hold a histogram of kmer counts.
/// For each count, the histogram stores the number of kmers (n_kmers) with that count.
#[derive(Debug, Clone)]
pub struct Histogram {
    histo: Vec<u64>, // for 0 through histo_max, the number of kmers with that count. 0 isn't actually expected to have a value.
    histo_large: rustc_hash::FxHashMap<u64, u64>, // for counts larger than histo_max, the number of kmers with that count
    histo_max: u64,
}

impl Histogram {
    pub fn new(histo_max: &u64) -> Histogram {
        let length = *histo_max as usize + 2;
        let histo: Vec<u64> = vec![0; length];
        let histo_large: rustc_hash::FxHashMap<u64, u64> = rustc_hash::FxHashMap::default();
        Histogram {
            histo,
            histo_large,
            histo_max: *histo_max,
        }
    }

    #[allow(dead_code)]
    pub fn ingest_kmer_counts(&mut self, kmer_counts: &KmerCounts) {
        for (_, count) in kmer_counts.iter() {
            if *count <= self.histo_max {
                self.histo[*count as usize] += 1;
            } else {
                let counter = self.histo_large.entry(*count).or_insert(0);
                *counter += 1;
            }
        }
    }

    /// Move a kmer from one count bin to another.
    /// Decrements the bin for `old_count` and increments the bin for `new_count`.
    pub(crate) fn move_count(&mut self, old_count: u64, new_count: u64) {
        // Decrement old bin
        if old_count > 0 {
            if old_count <= self.histo_max {
                self.histo[old_count as usize] -= 1;
            } else if let Some(n) = self.histo_large.get_mut(&old_count) {
                *n -= 1;
                if *n == 0 {
                    self.histo_large.remove(&old_count);
                }
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
        let mut sum: u64 = 0;
        for i in 1..self.histo.len() {
            sum += self.histo[i] * i as u64;
        }

        for (count, n_kmers) in self.histo_large.iter() {
            sum += count * n_kmers;
        }

        sum
    }

    pub fn get_n_unique_kmers(&self) -> u64 {
        let mut sum: u64 = 0;
        for i in 1..self.histo.len() {
            sum += self.histo[i];
        }

        for (_count, n_kmers) in self.histo_large.iter() {
            sum += n_kmers;
        }

        sum
    }

    pub fn get_vector(&self) -> Result<Vec<u64>> {
        let mut histo_vec = self.histo.clone();
        let last = histo_vec.last_mut().context("Histogram vector is empty")?;

        for (_, n_kmers) in self.histo_large.iter() {
            *last += n_kmers;
        }

        Ok(histo_vec)
    }
}
