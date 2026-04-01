// pcr/read_filter.rs — Primer-kmer-based read filtering for per-gene read selection
//
// Filters reads to only those containing primer-derived kmers,
// reducing the read set per gene before threading.

use std::collections::HashSet;

use crate::io::ReadRecord;
use crate::kmer::KmerCounts;
use crate::kmer::encoding::kmers_from_ascii;

/// Filters reads by presence of primer-derived kmers.
/// Designed for reuse in Phase 7 (read-backed runway) where
/// reads matching specific primer kmers are needed before the
/// full graph exists.
pub struct PrimerReadFilter {
    /// Union of all primer-derived canonical kmers
    primer_kmers: HashSet<u64>,
    k: usize,
}

impl PrimerReadFilter {
    /// Build a filter from the forward and reverse primer KmerCounts.
    pub fn from_primer_kmers(
        forward_primer_kmers: &KmerCounts,
        reverse_primer_kmers: &KmerCounts,
        k: usize,
    ) -> Self {
        let mut primer_kmers = HashSet::new();

        for (&kmer, _) in forward_primer_kmers.iter() {
            // KmerCounts stores canonical kmers
            primer_kmers.insert(kmer);
        }
        for (&kmer, _) in reverse_primer_kmers.iter() {
            primer_kmers.insert(kmer);
        }

        PrimerReadFilter { primer_kmers, k }
    }

    /// Check whether a read contains any primer kmer.
    pub fn matches(&self, sequence: &str) -> bool {
        let kmers = match kmers_from_ascii(sequence, self.k) {
            Ok(k) => k,
            Err(_) => return false,
        };
        kmers.iter().any(|kmer| self.primer_kmers.contains(kmer))
    }

    /// Filter a slice of reads to only those matching primer kmers.
    pub fn filter_reads<'a>(&self, reads: &'a [ReadRecord]) -> Vec<&'a ReadRecord> {
        reads.iter().filter(|r| self.matches(&r.sequence)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::Mate;

    #[test]
    fn test_empty_filter() {
        let fwd = KmerCounts::new(&3);
        let rev = KmerCounts::new(&3);
        let filter = PrimerReadFilter::from_primer_kmers(&fwd, &rev, 3);
        assert!(!filter.matches("ACGTACGT"));
    }
}
