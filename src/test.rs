use super::*;

fn revcomp_kmer(kmer: u64, k: u8) -> u64 {
    let mut revcomp = 0;
    for i in 0..k {
        let base = (kmer >> (2 * i)) & 3;
        revcomp = (revcomp << 2) | (3 - base);
    }
    revcomp
}

#[test]
fn test_tests() {
    assert_eq!(2 + 2, 4);
}

#[test]
fn test_seq_to_ints() {
    // 'A' => 0, // 00
    // 'C' => 1, // 01
    // 'G' => 2, // 10
    // 'T' => 3, // 11
    let seq = "CGTAATGCGGCGA";
    let expected = vec![vec![0b01101100, 0b00111001, 0b10100110]];
    let actual = seq_to_ints(seq);
    assert_eq!(actual, expected);
}

#[test]
fn test_seq_to_ints_n() {
    // 'A' => 0, // 00
    // 'C' => 1, // 01
    // 'G' => 2, // 10
    // 'T' => 3, // 11
    let seq = "CGTANATGCGGCGA";
    let expected = vec![vec![0b01101100], vec![0b00111001, 0b10100110]];
    let actual = seq_to_ints(seq);
    assert_eq!(actual, expected);
}

#[test]
fn test_revcomp_kmer() {
    let kmer = 0b0010_0110;
    let actual = revcomp_kmer(kmer, 3);

    // Check that the reverse complement of the reverse complement is the original
    assert_eq!(kmer, revcomp_kmer(actual, 3));

    // Check against hard coded expected value
    let expected = 0b0001_1001;
    assert_eq!(actual, expected);

    let kmer = 0b0110_1100_0011_1001_1010_0110;
    let actual = revcomp_kmer(kmer, 12);
    assert_eq!(kmer, revcomp_kmer(actual, 12));

    let expected = 0b0110_0101_1001_0011_1100_0110;

    assert_eq!(actual, expected);
}

#[test]
fn test_ints_to_kmers() {
    let ints = vec![0b01101100, 0b00111001, 0b10100110];
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
    let actual = ints_to_kmers(&ints, 9);
    assert_eq!(actual, expected);
}
