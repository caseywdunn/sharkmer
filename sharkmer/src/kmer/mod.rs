// kmer/mod.rs
//! This module provides kmer functions.

#[cfg(feature = "nohashmap")]
use std::collections::HashMap;

#[cfg(feature = "nohashmap")]
use std::hash::{BuildHasherDefault, Hasher};

#[cfg(feature = "intmap")]
use intmap::{IntMap, Entry};

#[cfg(feature = "fxhashmap")]
use rustc_hash::FxHashMap;




#[cfg(feature = "nohashmap")]
type NoHashHashMap<K, V> = HashMap<K, V, BuildHasherDefault<NoHashHasher>>;

#[cfg(feature = "nohashmap")]
struct NoHashHasher(u64); // Add a field to store the hash

#[cfg(feature = "nohashmap")]
impl Hasher for NoHashHasher {
    fn finish(&self) -> u64 {
        self.0 // Return the stored hash
    }

    fn write(&mut self, _: &[u8]) {
        // This method is required but won't be used
    }

    fn write_u64(&mut self, i: u64) {
        self.0 = i; // Store the key as the hash
    }
}

#[cfg(feature = "nohashmap")]
impl Default for NoHashHasher {
    fn default() -> Self {
        NoHashHasher(0) // Initialize the hash with 0
    }
}

/// A structure to hold a read in two bit encoding of bases
/// * `00` represents `A`
/// * `01` represents `C`
/// * `10` represents `G`
/// * `11` represents `T`
///
/// If the length is not divisible by 4, the least significant bits of the last byte will be ignored.
#[derive(Debug)]
pub struct Read {
    pub sequence: Vec<u8>,
    pub length: usize, // Number of bases in the read
}

impl Read {
    pub fn new(sequence: Vec<u8>, length: usize) -> Read {
        Read { sequence, length }
    }

    /// This function encodes a DNA sequence into a series of 8-bit integers,
    /// where each integer represents 4 consecutive bases from the sequence.
    /// The encoding uses a 2-bit representation for each base:
    ///
    /// * `A` -> `00`
    /// * `C` -> `01`
    /// * `G` -> `10`
    /// * `T` -> `11`
    ///
    /// Any other character will cause the function to panic.
    pub fn from_str(seq: &str) -> Read {
        let mut ints: Vec<u8> = Vec::with_capacity(seq.len() / 4 + 1);
        let mut frame: u8 = 0; // The current 8-bit integer being constructed
        let mut length: usize = 0; // length of the current subread
        for c in seq.chars() {
            length += 1;
            let base = match c {
                'A' => 0, // 00
                'C' => 1, // 01
                'G' => 2, // 10
                'T' => 3, // 11
                _ => panic!(
                    "Invalid character {} in sequence {}. Only ACGT allowed.",
                    c, seq
                ),
            };

            frame = (frame << 2) | base;
            if (length) % 4 == 0 {
                ints.push(frame);
                frame = 0; // Reset frame after pushing to the vector
            }
        }

        // Check if there is anything left in the frame, and rotate it into place if so
        let modulo = length % 4;
        if modulo != 0 {
            let shift = 2 * (4 - modulo);
            frame <<= shift;
            ints.push(frame);
        }

        // Create and return the read
        Read::new(ints, length)
    }

    pub fn validate(&self) -> bool {
        let mut valid = true;
        if self.length % 4 == 0 {
            if self.sequence.len() != self.length / 4 {
                valid = false;
            }
        } else {
            if self.sequence.len() != (self.length / 4 + 1) {
                valid = false;
            }
        }

        valid
    }

    /// Get canonical kmers from a Read.
    ///
    /// Given a sequence that's been encoded as a series of 8-bit integers (where each integer
    /// represents 4 bases), this function extracts canonical kmers from it. A canonical kmer
    /// is the lexicographically smaller of a kmer and its reverse complement.
    ///
    /// The encoding uses a 2-bit representation for each base:
    ///
    /// * `A` -> `00`
    /// * `C` -> `01`
    /// * `G` -> `10`
    /// * `T` -> `11`
    ///
    /// # Arguments
    ///
    /// * `ints` - A reference to the vector of 8-bit integers which encode a DNA sequence.
    /// * `k` - A reference to the length of the desired kmers (i.e., the number of bases in each kmer).
    ///
    /// # Returns
    ///
    /// A vector of 64-bit integers, where each integer represents a canonical kmer extracted
    /// from the encoded sequence.
    ///
    /// # Example
    ///
    /// ```rust
    /// let encoded_sequence: Vec<u8> = /* some encoded sequence */;
    /// let kmer_length: usize = /* desired kmer length */;
    /// let canonical_kmers = get_kmers(&encoded_sequence, &kmer_length);
    /// ```
    pub fn get_kmers(&self, k: &usize) -> Vec<u64> {
        let ints = &self.sequence;
        let mut kmers: Vec<u64> = Vec::with_capacity((ints.len() * 4 / k) + 1);
        let mut frame: u64 = 0; // read the bits for each base into the least significant end of this integer
        let mut revframe: u64 = 0; // read the bits for complement into the least significant end of this integer
        let mut n_valid = 0; // number of valid bases in the frame
        let mask: u64 = (1 << (2 * k)) - 1;

        // No kmers if the read is shorter than k
        if self.length < *k {
            return kmers;
        }

        // Iterate over the bases
        for (_i, &int) in ints.iter().enumerate() {
            // Iterate over the bases in the integer
            for j in 0..4 {
                // Get the base from the left side of the integer,
                // move it to the least two significant bits and mask it
                let base = ((int >> ((3 - j) * 2)) & 3) as u64;

                frame = (frame << 2) | base;
                revframe = (revframe >> 2) | ((3 - base) << (2 * (k - 1)));
                n_valid += 1;
                if n_valid >= *k {
                    let forward = frame & mask;
                    let reverse = revframe & mask;

                    // assert_eq!(forward, revcomp_kmer(reverse, k));

                    if forward < reverse {
                        kmers.push(forward);
                    } else {
                        kmers.push(reverse);
                    }
                }
            }
        }

        // If the read length is not divisible by 4, the last byte will have some extra bases
        // that result in extra kmers we need to remove
        let modulo = self.length % 4;
        if modulo != 0 {
            let n_extra_bases = 4 - modulo;
            kmers.truncate(kmers.len() - n_extra_bases);
        }

        if kmers.len() != self.length - *k + 1 {
            panic!(
                "Number of kmers ({}) does not match expected ({}) for read of length {}",
                kmers.len(),
                self.length - *k + 1,
                self.length
            );
        }

        kmers
    }
}

impl PartialEq for Read {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence && self.length == other.length
    }
}

#[cfg(feature = "intmap")]
type MapType = IntMap<u64>;

#[cfg(feature = "fxhashmap")]
type MapType = rustc_hash::FxHashMap<u64, u64>;

#[cfg(feature = "nohashmap")]
type MapType = NoHashHashMap<u64, u64>;

pub struct KmerCounts {
    kmers: MapType,
    k: usize,
}

#[cfg(feature = "intmap")]
impl KmerCounts {
    pub fn new(k: &usize) -> KmerCounts {
        KmerCounts {
            kmers: IntMap::new(),
            k: *k,
        }
    }

    pub fn ingest_reads(&mut self, reads: &[Read]) {
        for read in reads {
            for kmer in read.get_kmers(&self.k) {
                let counter = match self.kmers.entry(kmer) {
                    Entry::Occupied(entry) => entry.into_mut(),
                    Entry::Vacant(entry) => entry.insert(0),
                };
                *counter += 1;
            }
        }
    }

    /// Adds a kmer and its count to the KmerCounts object.
    pub fn insert(&mut self, kmer: &u64, count: &u64) {
        let counter = match self.kmers.entry(*kmer) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => entry.insert(0),
        };
        *counter += count;
    }

    /// Adds the counts from another KmerCounts object to this one.
    pub fn extend(&mut self, other: &KmerCounts) {
        if self.k != other.k {
            panic!("Cannot extend KmerCounts with different k");
        }

        for (kmer, count) in other.iter() {
            let counter = match self.kmers.entry(*kmer) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => entry.insert(0),
            };
            *counter += count;
        }
    }

    pub fn get_canonical(&self, kmer: &u64) -> u64 {
        let revcomp = revcomp_kmer(kmer, &self.k);
        let canonical = if *kmer < revcomp { *kmer } else { revcomp };
        *self.kmers.get(canonical).unwrap_or(&0)
    }

    pub fn contains(&self, kmer: &u64) -> bool {
        self.kmers.contains_key(*kmer)
    }

    pub fn remove_low_count_kmers(&mut self, min_count: &u64) {
        self.kmers.retain(|_, count| count >= min_count);
    }
}

#[cfg(feature = "fxhashmap")]
impl KmerCounts {
    pub fn new(k: &usize) -> KmerCounts {
        KmerCounts {
            kmers: FxHashMap::default(),
            k: *k,
        }
    }

    pub fn ingest_reads(&mut self, reads: &[Read]) {
        for read in reads {
            for kmer in read.get_kmers(&self.k) {
                let count = self.kmers.entry(kmer).or_insert(0);
                *count += 1;
            }
        }
    }

    /// Adds a kmer and its count to the KmerCounts object.
    pub fn insert(&mut self, kmer: &u64, count: &u64) {
        let old_count = self.kmers.entry(*kmer).or_insert(0);
        *old_count += count;
    }

    /// Adds the counts from another KmerCounts object to this one.
    pub fn extend(&mut self, other: &KmerCounts) {
        if self.k != other.k {
            panic!("Cannot extend KmerCounts with different k");
        }

        for (kmer, other_count) in other.iter() {
            let count = self.kmers.entry(*kmer).or_insert(0);
            *count += other_count;
        }
    }

    pub fn get_canonical(&self, kmer: &u64) -> u64 {
        let revcomp = revcomp_kmer(kmer, &self.k);
        let canonical = if *kmer < revcomp { *kmer } else { revcomp };
        *self.kmers.get(&canonical).unwrap_or(&0)
    }

    pub fn contains(&self, kmer: &u64) -> bool {
        self.kmers.contains_key(kmer)
    }

    pub fn remove_low_count_kmers(&mut self, min_count: &u64) {
        let mut min = *min_count;
        let min_ref = &mut min;
        self.kmers.retain(|_, count| count >= min_ref);
    }
}

#[cfg(feature = "nohashmap")]
impl KmerCounts {
    pub fn new(k: &usize) -> KmerCounts {
        KmerCounts {
            kmers: NoHashHashMap::<u64, u64>::default(),
            k: *k,
        }
    }

    pub fn ingest_reads(&mut self, reads: &[Read]) {
        for read in reads {
            for kmer in read.get_kmers(&self.k) {
                let count = self.kmers.entry(kmer).or_insert(0);
                *count += 1;
            }
        }
    }

    /// Adds a kmer and its count to the KmerCounts object.
    pub fn insert(&mut self, kmer: &u64, count: &u64) {
        let old_count = self.kmers.entry(*kmer).or_insert(0);
        *old_count += count;
    }

    /// Adds the counts from another KmerCounts object to this one.
    pub fn extend(&mut self, other: &KmerCounts) {
        if self.k != other.k {
            panic!("Cannot extend KmerCounts with different k");
        }

        for (kmer, other_count) in other.iter() {
            let count = self.kmers.entry(*kmer).or_insert(0);
            *count += other_count;
        }
    }

    pub fn get_canonical(&self, kmer: &u64) -> u64 {
        let revcomp = revcomp_kmer(kmer, &self.k);
        let canonical = if *kmer < revcomp { *kmer } else { revcomp };
        *self.kmers.get(&canonical).unwrap_or(&0)
    }

    pub fn contains(&self, kmer: &u64) -> bool {
        self.kmers.contains_key(kmer)
    }

    pub fn remove_low_count_kmers(&mut self, min_count: &u64) {
        let mut min = *min_count;
        let min_ref = &mut min;
        self.kmers.retain(|_, count| count >= min_ref);
    }
}

// Generic methods that don't depend on the hash map implementation
impl KmerCounts {
    pub fn iter(&self) -> impl Iterator<Item = (&u64, &u64)> {
        self.kmers.iter()
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    pub fn get_k(&self) -> usize {
        self.k
    }

    pub fn get_n_kmers(&self) -> u64 {
        self.kmers.iter().map(|(_, count)| count).sum()
    }

    pub fn get_n_unique_kmers(&self) -> u64 {
        self.kmers.len() as u64
    }

    pub fn counts(&self) -> Vec<u64> {
        self.kmers.values().map(|&count| count).collect()
    }

    pub fn get_count(&self, kmer: &u64) -> u64 {
        *self.kmers.get(kmer).unwrap_or(&0)
    }

    pub fn kmers(&self) -> Vec<u64> {
        self.kmers.keys().map(|&kmer| kmer).collect()
    }

    pub fn get_max_count(&self) -> u64 {
        *self.kmers.values().max().unwrap_or(&0)
    }

    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }

    // Expand the KmerCounts object to include the reverse complements of all the kmers, not just canonical kmers
    pub fn add_reverse_complements(&mut self) {
        let mut new_kmers = KmerCounts::new(&self.k);
        for (kmer, count) in self.iter() {
            new_kmers.insert(&crate::kmer::revcomp_kmer(kmer, &self.k), count);
        }

        // TODO - will increase counts if kmers already present, should
        // check if they are there before extending
        self.extend(&new_kmers);
    }

    // Create a new KmerCounts object that has kmers with at least min_count and 
    // includes the reverse complements of all the kmers, not just canonical kmers
    pub fn get_pcr_kmers(&self, min_count: &u64) -> KmerCounts {
        let mut pcr_kmers = KmerCounts::new(&self.k);
        for (kmer, count) in self.iter() {
            if count >= min_count {
                pcr_kmers.insert(kmer, count);
                pcr_kmers.insert(&crate::kmer::revcomp_kmer(kmer, &self.k), count);
            }
        }
        pcr_kmers
    }
}

impl Clone for KmerCounts {
    fn clone(&self) -> Self {
        KmerCounts {
            kmers: self.kmers.clone(),
            k: self.k,
        }
    }
}

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

    pub fn ingest_reads(&mut self, reads: &[Read]) {
        self.kmer_counts.ingest_reads(reads);
        self.n_reads += reads.len() as u64;
        self.n_bases += reads.iter().map(|read| read.length as u64).sum::<u64>();
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

/// A structure to hold a histogram of kmer counts.
/// For each count, the histogram stores the number of kmers (n_kmers) with that count.
#[derive(Debug)]
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

    pub fn from_kmer_counts(kmer_counts: &KmerCounts, histo_max: &u64) -> Histogram {
        let mut histo = Histogram::new(histo_max);
        histo.ingest_kmer_counts(kmer_counts);
        histo
    }

    pub fn _get(&self, count: &u64) -> u64 {
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

    pub fn get_vector(&self) -> Vec<u64> {
        let mut histo_vec = self.histo.clone();
        let last = histo_vec.last_mut().unwrap();

        for (_, n_kmers) in self.histo_large.iter() {
            *last += n_kmers;
        }

        histo_vec
    }
}

/// Returns the reverse complement of a specified kmer.
///
/// The function computes the reverse complement of a kmer represented as a 64-bit unsigned integer.
/// Each base in the kmer is encoded using 2 bits, with the following scheme:
///
/// * `00` represents `A`
/// * `01` represents `C`
/// * `10` represents `G`
/// * `11` represents `T`
///
/// # Arguments
///
/// * `kmer` - A reference to the 64-bit unsigned integer representation of the kmer.
/// * `k` - A reference to the length of the kmer (i.e., the number of bases).
///
/// # Returns
///
/// A 64-bit unsigned integer representing the reverse complement of the given kmer.
///
/// # Example
///
/// ```rust
/// let kmer: u64 = /* some kmer encoding */;
/// let k: usize = /* kmer length */;
/// let reverse_complement = revcomp_kmer(&kmer, &k);
/// ```
pub fn revcomp_kmer(kmer: &u64, k: &usize) -> u64 {
    let mut revcomp = 0;
    for i in 0..*k {
        let base = (kmer >> (2 * i)) & 3;
        revcomp = (revcomp << 2) | (3 - base);
    }
    revcomp
}

/// Converts a DNA sequence into a vector of Reads (really subreads).
///
///
/// Any sequence containing the base `N` is split into multiple subsequences at that point,
/// and each sub-sequence is encoded separately. The result is a vector of Reads,
/// with each Read holding the integer representation of a sub-sequence.
///
/// # Arguments
///
/// * `seq` - A reference to the DNA sequence string to be encoded.
///
/// # Returns
///
/// A vector of Reads, where each inner vector represents an encoded
/// DNA sub-sequence.
///
/// # Panics
///
/// The function will exit prematurely if it encounters a base other than `A`, `C`, `G`, `T`, or `N`.
///
/// # Example
///
/// ```rust
/// let sequence = "ACGTNAGCT";
/// let reads = seq_to_reads(&sequence);
/// ```
pub fn seq_to_reads(seq: &str) -> Vec<Read> {
    let mut reads: Vec<Read> = Vec::new();

    // Split seq on N
    for subseq in seq.split('N') {
        if !subseq.is_empty() {
            let read = Read::from_str(subseq);
            if !read.validate() {
                panic!("Invalid read: {:?} from subsequence {}", read, subseq);
            }
            reads.push(read);
        }
    }
    reads
}

pub fn kmer_to_seq(kmer: &u64, k: &usize) -> String {
    let mut seq = String::new();
    for i in 0..*k {
        let base = (kmer >> (2 * (k - i - 1))) & 3;
        let base = match base {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => panic!("Invalid base"),
        };
        seq.push(base);
    }
    seq
}

#[allow(dead_code)]
pub fn seq_to_kmer(seq: &str) -> u64 {
    let mut kmer = 0;
    for (_i, c) in seq.chars().enumerate() {
        let base = match c {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => panic!("Invalid base"),
        };
        kmer = (kmer << 2) | base;
    }
    kmer
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tests() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn test_read_validate() {
        let good_read = Read::new(vec![0b01101100, 0b00111001, 0b10100110, 0b00000000], 13);
        assert!(good_read.validate());

        let bad_read1 = Read::new(vec![0b01101100, 0b00111001, 0b10100110, 0b00000000], 10);
        assert!(!bad_read1.validate());

        let bad_read2 = Read::new(vec![0b01101100, 0b00111001, 0b10100110, 0b00000000], 17);
        assert!(!bad_read2.validate());
    }

    #[test]
    fn test_read_parsing() {
        let seq = "TANCACN";
        let reads = seq_to_reads(seq);
        for read in reads {
            assert!(read.validate());
        }

        let seq = "NTANCACNAGAAAATC";
        let reads = seq_to_reads(seq);
        for read in reads {
            assert!(read.validate());
        }

        let seq = "TATTAGCTCATCTANAACAATGAAAAATTGCATTGGCTNTAACTATGGATTTNTTAGAAATTAGTATTNATTTATCATTTTTAATTGGCATTATTNAACTCTTAAGAATAGATNGGAGTTCNCAATTAATTGAAGNTANCACNAGAAAATC";
        let reads = seq_to_reads(seq);
        for read in reads {
            assert!(read.validate());
        }
    }

    #[test]
    fn test_seq_to_reads() {
        // 'A' => 0, // 00
        // 'C' => 1, // 01
        // 'G' => 2, // 10
        // 'T' => 3, // 11
        let seq = "CGTAATGCGGCGA";
        let expected = vec![Read::new(
            vec![0b01101100, 0b00111001, 0b10100110, 0b00000000],
            13,
        )];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);

        let seq = "C";
        let expected = vec![Read::new(vec![0b01000000], 1)];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);

        let seq = "CGTAATGCGGCG";
        let expected = vec![Read::new(vec![0b01101100, 0b00111001, 0b10100110], 12)];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_from_str() {
        // 'A' => 0, // 00
        // 'C' => 1, // 01
        // 'G' => 2, // 10
        // 'T' => 3, // 11
        let seq = "CGTAATGCGGCGA";
        let expected = Read::new(vec![0b01101100, 0b00111001, 0b10100110, 0b00000000], 13);
        let actual = Read::from_str(seq);
        assert_eq!(actual, expected);

        let seq = "C";
        let expected = Read::new(vec![0b01000000], 1);
        let actual = Read::from_str(seq);
        assert_eq!(actual, expected);

        let seq = "CGTAATGCGGCG";
        let expected = Read::new(vec![0b01101100, 0b00111001, 0b10100110], 12);
        let actual = Read::from_str(seq);
        assert_eq!(actual, expected);

        let seq = "";
        let expected = Read::new(vec![], 0);
        let actual = Read::from_str(seq);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_seq_to_reads_n() {
        // 'A' => 0, // 00
        // 'C' => 1, // 01
        // 'G' => 2, // 10
        // 'T' => 3, // 11

        let seq = "NCGTAATGCGGCG";
        let expected = vec![Read::new(vec![0b01101100, 0b00111001, 0b10100110], 12)];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);

        let seq = "CGTANATGCGGCGA";
        let expected = vec![
            Read::new(vec![0b01101100], 4),
            Read::new(vec![0b00111001, 0b10100110, 0b00000000], 9),
        ];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);

        let seq = "NCGTANATGCGGCGA";
        let expected = vec![
            Read::new(vec![0b01101100], 4),
            Read::new(vec![0b00111001, 0b10100110, 0b00000000], 9),
        ];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);

        let seq = "NCGTANATGCGGCGANN";
        let expected = vec![
            Read::new(vec![0b01101100], 4),
            Read::new(vec![0b00111001, 0b10100110, 0b00000000], 9),
        ];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);

        let seq = "NNCGTANATGCGGCGA";
        let expected = vec![
            Read::new(vec![0b01101100], 4),
            Read::new(vec![0b00111001, 0b10100110, 0b00000000], 9),
        ];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_revcomp_kmer() {
        let kmer = 0b0010_0110;
        let actual = revcomp_kmer(&kmer, &3usize);

        // Check that the reverse complement of the reverse complement is the original
        assert_eq!(kmer, revcomp_kmer(&actual, &3usize));

        // Check against hard coded expected value
        let expected = 0b0001_1001;
        assert_eq!(actual, expected);

        let kmer = 0b0110_1100_0011_1001_1010_0110;
        let actual = revcomp_kmer(&kmer, &12usize);
        assert_eq!(kmer, revcomp_kmer(&actual, &12usize));

        let expected = 0b0110_0101_1001_0011_1100_0110;

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_get_kmers() {
        let ints = vec![0b01101100, 0b00111001, 0b10100110];
        let read = Read::new(ints.clone(), 12);
        // original	                // reverse complement
        // 0b01101100_00111001_10   0b01_10010011_11000110  >
        // 0b101100_00111001_1010   0b0101_10010011_110001  >
        // 0b1100_00111001_101001   0b100101_10010011_1100  >
        // 0b00_00111001_10100110   0b01100101_10010011_11  <
        //
        let expected = vec![
            0b01_1001_0011_1100_0110,
            0b01_0110_0100_1111_0001,
            0b10_0101_1001_0011_1100,
            0b00_0011_1001_1010_0110,
        ];
        let actual = read.get_kmers(&9usize);
        assert_eq!(actual, expected);

        // Test cases where read length not divisible by 4
        let read = Read::new(ints.clone(), 11);
        let expected = vec![
            0b01_1001_0011_1100_0110,
            0b01_0110_0100_1111_0001,
            0b10_0101_1001_0011_1100,
        ];
        let actual = read.get_kmers(&9usize);
        assert_eq!(actual, expected);

        // Test cases where read length not divisible by 4
        let read = Read::new(ints.clone(), 10);
        let expected = vec![0b01_1001_0011_1100_0110, 0b01_0110_0100_1111_0001];
        let actual = read.get_kmers(&9usize);
        assert_eq!(actual, expected);

        // Test cases where read length is k
        let read = Read::new(ints.clone(), 9);
        let expected = vec![0b01_1001_0011_1100_0110];
        let actual = read.get_kmers(&9usize);
        assert_eq!(actual, expected);

        // Test cases where read length is < k
        let ints_short = vec![0b01101100, 0b00111001];
        let read = Read::new(ints_short, 8);
        let expected = vec![];
        let actual = read.get_kmers(&9usize);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_kmer_to_seq() {
        let kmer = 0b1001_1000;
        let seq = kmer_to_seq(&kmer, &4usize);
        assert_eq!(seq, "GCGA");

        let kmer = 0b1001_1000_1001_1000;
        let seq = crate::kmer::kmer_to_seq(&kmer, &8usize);
        assert_eq!(seq, "GCGAGCGA");
    }

    #[test]
    fn test_histogram() {
        let mut kmers = KmerCounts::new(&11usize);
        kmers.insert(&1, &5);
        kmers.insert(&20, &5);
        kmers.insert(&2, &7);
        kmers.insert(&11, &11);
        kmers.insert(&12, &12);

        let histo_max = 10;
        let histo = Histogram::from_kmer_counts(&kmers, &histo_max);

        let histo_vec = histo.get_vector();

        assert_eq!(histo_vec.len(), histo_max as usize + 2);
        let expected_histo: Vec<u64> = vec![0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 2];
        assert_eq!(histo_vec, expected_histo);
    }
}
