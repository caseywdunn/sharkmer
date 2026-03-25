// kmer/mod.rs
//! This module provides kmer functions.

pub mod chunk;
pub mod counting;
pub mod encoding;
pub mod histogram;

// Re-export public items so that `crate::kmer::X` continues to work.
pub use chunk::Chunk;
pub use counting::{FilteredKmerCounts, KmerCounts};
#[allow(unused_imports)]
pub use encoding::{Read, kmer_to_seq, revcomp_kmer, seq_to_kmer, seq_to_reads};
pub use histogram::Histogram;

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
        let reads = seq_to_reads(seq).unwrap();
        for read in reads {
            assert!(read.validate());
        }

        let seq = "NTANCACNAGAAAATC";
        let reads = seq_to_reads(seq).unwrap();
        for read in reads {
            assert!(read.validate());
        }

        let seq = "TATTAGCTCATCTANAACAATGAAAAATTGCATTGGCTNTAACTATGGATTTNTTAGAAATTAGTATTNATTTATCATTTTTAATTGGCATTATTNAACTCTTAAGAATAGATNGGAGTTCNCAATTAATTGAAGNTANCACNAGAAAATC";
        let reads = seq_to_reads(seq).unwrap();
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
        let actual = seq_to_reads(seq).unwrap();
        assert_eq!(actual, expected);

        let seq = "C";
        let expected = vec![Read::new(vec![0b01000000], 1)];
        let actual = seq_to_reads(seq).unwrap();
        assert_eq!(actual, expected);

        let seq = "CGTAATGCGGCG";
        let expected = vec![Read::new(vec![0b01101100, 0b00111001, 0b10100110], 12)];
        let actual = seq_to_reads(seq).unwrap();
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
        let actual = Read::from_str(seq).unwrap();
        assert_eq!(actual, expected);

        let seq = "C";
        let expected = Read::new(vec![0b01000000], 1);
        let actual = Read::from_str(seq).unwrap();
        assert_eq!(actual, expected);

        let seq = "CGTAATGCGGCG";
        let expected = Read::new(vec![0b01101100, 0b00111001, 0b10100110], 12);
        let actual = Read::from_str(seq).unwrap();
        assert_eq!(actual, expected);

        let seq = "";
        let expected = Read::new(vec![], 0);
        let actual = Read::from_str(seq).unwrap();
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
        let actual = seq_to_reads(seq).unwrap();
        assert_eq!(actual, expected);

        let seq = "CGTANATGCGGCGA";
        let expected = vec![
            Read::new(vec![0b01101100], 4),
            Read::new(vec![0b00111001, 0b10100110, 0b00000000], 9),
        ];
        let actual = seq_to_reads(seq).unwrap();
        assert_eq!(actual, expected);

        let seq = "NCGTANATGCGGCGA";
        let expected = vec![
            Read::new(vec![0b01101100], 4),
            Read::new(vec![0b00111001, 0b10100110, 0b00000000], 9),
        ];
        let actual = seq_to_reads(seq).unwrap();
        assert_eq!(actual, expected);

        let seq = "NCGTANATGCGGCGANN";
        let expected = vec![
            Read::new(vec![0b01101100], 4),
            Read::new(vec![0b00111001, 0b10100110, 0b00000000], 9),
        ];
        let actual = seq_to_reads(seq).unwrap();
        assert_eq!(actual, expected);

        let seq = "NNCGTANATGCGGCGA";
        let expected = vec![
            Read::new(vec![0b01101100], 4),
            Read::new(vec![0b00111001, 0b10100110, 0b00000000], 9),
        ];
        let actual = seq_to_reads(seq).unwrap();
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
        let actual = read.get_kmers(&9usize).unwrap();
        assert_eq!(actual, expected);

        // Test cases where read length not divisible by 4
        let read = Read::new(ints.clone(), 11);
        let expected = vec![
            0b01_1001_0011_1100_0110,
            0b01_0110_0100_1111_0001,
            0b10_0101_1001_0011_1100,
        ];
        let actual = read.get_kmers(&9usize).unwrap();
        assert_eq!(actual, expected);

        // Test cases where read length not divisible by 4
        let read = Read::new(ints.clone(), 10);
        let expected = vec![0b01_1001_0011_1100_0110, 0b01_0110_0100_1111_0001];
        let actual = read.get_kmers(&9usize).unwrap();
        assert_eq!(actual, expected);

        // Test cases where read length is k
        let read = Read::new(ints.clone(), 9);
        let expected = vec![0b01_1001_0011_1100_0110];
        let actual = read.get_kmers(&9usize).unwrap();
        assert_eq!(actual, expected);

        // Test cases where read length is < k
        let ints_short = vec![0b01101100, 0b00111001];
        let read = Read::new(ints_short, 8);
        let expected = vec![];
        let actual = read.get_kmers(&9usize).unwrap();
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

        let histo_vec = histo.get_vector().unwrap();

        assert_eq!(histo_vec.len(), histo_max as usize + 2);
        let expected_histo: Vec<u64> = vec![0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 2];
        assert_eq!(histo_vec, expected_histo);
    }
}
