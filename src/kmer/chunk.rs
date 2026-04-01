// kmer/chunk.rs
//! Chunk: groups reads for incremental kmer counting.

use anyhow::Result;

use super::counting::KmerCounts;
use super::encoding::{Read, count_valid_bases};

pub struct Chunk {
    kmer_counts: KmerCounts,
    n_reads: u64,
    n_bases: u64,
}

impl Chunk {
    pub fn new(k: &usize) -> Chunk {
        Chunk {
            kmer_counts: KmerCounts::new(k),
            n_reads: 0,
            n_bases: 0,
        }
    }

    #[allow(dead_code)]
    pub fn ingest_reads(&mut self, reads: &[Read]) -> Result<()> {
        self.kmer_counts.ingest_reads(reads)?;
        self.n_reads += reads.len() as u64;
        self.n_bases += reads.iter().map(|read| read.length as u64).sum::<u64>();
        Ok(())
    }

    /// Ingest a single FASTQ sequence directly from ASCII, bypassing Read encoding.
    pub fn ingest_seq(&mut self, seq: &str) -> Result<()> {
        self.kmer_counts.ingest_seq(seq)?;
        self.n_reads += 1;
        self.n_bases += count_valid_bases(seq);
        Ok(())
    }

    /// Ingest a sequence and check each kmer against an Oligo filter.
    /// Returns true if any kmer matched the filter.
    pub fn ingest_seq_with_filter(
        &mut self,
        seq: &str,
        filter: &crate::io::OligoFilter,
    ) -> Result<bool> {
        let matched = self.kmer_counts.ingest_seq_with_filter(seq, filter)?;
        self.n_reads += 1;
        self.n_bases += count_valid_bases(seq);
        Ok(matched)
    }

    pub fn get_kmer_counts(&self) -> &KmerCounts {
        &self.kmer_counts
    }

    pub fn get_n_kmers(&self) -> u64 {
        self.kmer_counts.get_n_kmers()
    }

    pub fn get_n_reads(&self) -> u64 {
        self.n_reads
    }

    pub fn get_n_bases(&self) -> u64 {
        self.n_bases
    }
}
