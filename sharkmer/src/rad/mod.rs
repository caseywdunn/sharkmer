use rustc_hash::FxHashMap;
use bio::io::fasta;

use crate::COLOR_FAIL;
use crate::COLOR_NOTE;
use crate::COLOR_SUCCESS;
use crate::COLOR_WARNING;

pub struct RADParams {
    pub cut1: String,
    pub cut2: String,
	pub min_length: usize,
    pub max_length: usize,
    pub name: String,
}

pub fn do_rad(
	kmer_counts_map: &FxHashMap<u64, u64>,
    k: &usize,
    sample_name: &str,
    verbosity: usize,
    params: &RADParams,
) -> Vec<bio::io::fasta::Record> {
	let records: Vec<fasta::Record> = Vec::new();

	records
}