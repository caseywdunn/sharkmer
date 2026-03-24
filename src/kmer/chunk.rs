// kmer/chunk.rs
//! Chunk: groups reads for incremental kmer counting.

use anyhow::Result;

use super::counting::KmerCounts;
use super::encoding::Read;

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

    pub fn ingest_reads(&mut self, reads: &[Read]) -> Result<()> {
        self.kmer_counts.ingest_reads(reads)?;
        self.n_reads += reads.len() as u64;
        self.n_bases += reads.iter().map(|read| read.length as u64).sum::<u64>();
        Ok(())
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
