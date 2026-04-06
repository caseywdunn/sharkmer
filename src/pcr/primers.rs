// pcr/primers.rs — primer preprocessing functions

use anyhow::{Result, bail, ensure};
use log::debug;
use std::collections::HashSet;

use crate::kmer::{FilteredKmerCounts, KmerCounts};

use super::{Oligo, PCRParams, PrimerDirection};

pub(super) fn is_valid_nucleotide(c: char) -> bool {
    match c {
        'A' => true,
        'C' => true,
        'G' => true,
        'T' => true,
        'R' => true, // A or G
        'Y' => true, // C or T
        'S' => true, // C or G
        'W' => true, // A or T
        'K' => true, // G or T
        'M' => true, // A or C
        'B' => true, // C or G or T
        'D' => true, // A or G or T
        'H' => true, // A or C or T
        'V' => true, // A or C or G
        'N' => true, // A or C or G or T
        _ => false,
    }
}

// Given an oligo sequence, return a Oligo struct representing it
pub fn string_to_oligo(seq: &str) -> Result<Oligo> {
    ensure!(
        seq.len() <= 32,
        "Oligo sequence length {} exceeds maximum of 32 bases",
        seq.len()
    );
    let mut kmer: u64 = 0;
    let mut length: usize = 0;
    for c in seq.chars() {
        let base = match c {
            'A' => 0, // 00
            'C' => 1, // 01
            'G' => 2, // 10
            'T' => 3, // 11
            _ => bail!("Invalid nucleotide {} in {}", c, seq),
        };

        kmer = (kmer << 2) | base as u64;
        length += 1;
    }
    Ok(Oligo { length, kmer })
}

// Given a primer that may include ambiguous nucleotides, return a set
// of sequences that include all possible resolutions of the ambiguity. If
// there are no ambiguous nucleotides, the set contains only the original
// sequence.
pub(super) fn resolve_primer(primer: &str) -> HashSet<String> {
    let mut sequences: HashSet<String> = HashSet::new();

    for nuc in primer.chars() {
        let possible_nucs = match nuc {
            'R' => vec!['A', 'G'],
            'Y' => vec!['C', 'T'],
            'S' => vec!['G', 'C'],
            'W' => vec!['A', 'T'],
            'K' => vec!['G', 'T'],
            'M' => vec!['A', 'C'],
            'B' => vec!['C', 'G', 'T'],
            'D' => vec!['A', 'G', 'T'],
            'H' => vec!['A', 'C', 'T'],
            'V' => vec!['A', 'C', 'G'],
            'N' => vec!['A', 'C', 'G', 'T'],
            _ => vec![nuc], // Use the nucleotide directly
        };

        if sequences.is_empty() {
            for possible_nuc in possible_nucs {
                sequences.insert(possible_nuc.to_string());
            }
        } else {
            let mut new_sequences = HashSet::new();
            for seq in sequences.iter() {
                for &possible_nuc in &possible_nucs {
                    let new_seq = format!("{}{}", seq, possible_nuc);
                    new_sequences.insert(new_seq);
                }
            }
            sequences = new_sequences;
        }
    }

    sequences
}

/// Given a set of sequences, return a set of all sequences that differ
/// from each original sequence at up to r positions. Includes the original
/// sequences.
pub(super) fn permute_sequences(sequences: HashSet<String>, r: &usize) -> HashSet<String> {
    let mut permutations = HashSet::new();

    for seq in &sequences {
        for positions in combinations(seq.len(), *r) {
            generate_recursive_permutations(seq, &positions, 0, &mut permutations);
        }
    }

    permutations
}

pub(super) fn combinations(n: usize, r: usize) -> Vec<Vec<usize>> {
    if r > n {
        return vec![];
    }

    if r == 0 {
        return vec![Vec::new()];
    }

    if n == r {
        return vec![(0..r).collect()];
    }

    let without_last = combinations(n - 1, r);
    let mut with_last = combinations(n - 1, r - 1);
    for item in &mut with_last {
        item.push(n - 1);
    }

    without_last.into_iter().chain(with_last).collect()
}

pub(super) fn generate_recursive_permutations(
    seq: &str,
    positions: &Vec<usize>,
    current: usize,
    unique_sequences: &mut HashSet<String>,
) {
    let nucleotides = ['A', 'T', 'C', 'G'];

    if current == positions.len() {
        unique_sequences.insert(seq.to_string());
        return;
    }

    let pos = positions[current];
    for &nucleotide in nucleotides.iter() {
        let mut new_seq = seq.chars().collect::<Vec<_>>();
        new_seq[pos] = nucleotide;
        generate_recursive_permutations(
            &new_seq.iter().collect::<String>(),
            positions,
            current + 1,
            unique_sequences,
        );
    }
}

// Find the kmers that contain the oligos.
// Always match the oligo at the START of the kmer.
fn find_oligos_in_kmers(
    oligos: &[Oligo],
    kmers: &FilteredKmerCounts,
    min_count: &u32,
) -> KmerCounts {
    assert!(
        !oligos.is_empty(),
        "find_oligos_in_kmers called with no oligos"
    );

    // Assume all oligos have the same length
    let oligo_length = oligos[0].length;
    let k = kmers.get_k();
    assert!(
        k <= 31,
        "k={} exceeds 31; 2-bit encoding requires k <= 31 for u64",
        k
    );
    assert!(
        oligo_length > 0 && oligo_length <= k,
        "oligo length {} out of range for k={} (must be 1..=k)",
        oligo_length,
        k
    );

    // Create hash set of oligos, shifted to the start (high-order bits) of the kmer
    let oligo_set: HashSet<u64> = oligos
        .iter()
        .map(|oligo| oligo.kmer << (2 * (k - oligo_length)))
        .collect();

    // Create mask covering the high-order bits where the oligo sits
    let mut mask: u64 = 0;
    for _i in 0..(2 * oligo_length) {
        mask = (mask << 1) | 1;
    }
    mask <<= 2 * k - 2 * oligo_length;

    let mut kmers_match: KmerCounts = KmerCounts::new(&k);

    // RC mask covers low bits (reverse complement of a start-matching kmer has primer in low bits)
    let rc_mask = (1u64 << (2 * oligo_length)) - 1;

    // Build the set of RC oligos (reverse complement of each oligo, unshifted)
    let rc_oligo_set: HashSet<u64> = oligos
        .iter()
        .map(|oligo| crate::kmer::revcomp_kmer(&oligo.kmer, &oligo_length))
        .collect();

    for (kmer, count) in kmers.iter() {
        if count >= min_count {
            if oligo_set.contains(&(kmer & mask)) {
                kmers_match.insert(kmer, count);
            } else if rc_oligo_set.contains(&(kmer & rc_mask)) {
                // The primer matches the reverse complement orientation of this kmer,
                // so insert the revcomp form for correct graph construction
                let rc_kmer = crate::kmer::revcomp_kmer(kmer, &k);
                kmers_match.insert(&rc_kmer, count);
            }
        }
    }

    kmers_match
}

// Individual steps in the PCR process

/// Ingest pcr parameters, and return a set of all possible variants of the
/// forward or reverse primer. The variants are generated by resolving
/// ambiguous nucleotides, and then permuting the sequences to include
/// up to `mismatches` mismatches.
pub(super) fn preprocess_primer(
    params: &PCRParams,
    dir: PrimerDirection,
    k: &usize,
) -> Result<HashSet<String>> {
    let mut primer = params.forward_seq.clone();
    if dir == PrimerDirection::Reverse {
        primer = params.reverse_seq.clone();
    }

    let mut trim = params.trim;
    if trim > *k {
        gene_warn!(
            params.gene_name,
            "Trim length ({}) exceeds k ({}), adjusting trim to {}",
            trim,
            k,
            k
        );
        trim = *k;
    }

    // Check if either is longer than trim, if so retain only the last trim nucleotides
    if primer.len() > trim {
        primer = primer[primer.len() - trim..].to_string();
        gene_info!(
            params.gene_name,
            "Trimming the primer to {} so that it is within the trim length of {}.",
            primer,
            trim
        );
    }

    // Expand ambiguous nucleotides
    let mut primer_variants = resolve_primer(&primer);

    const MAX_RESOLVED_VARIANTS: usize = 10_000;
    if primer_variants.len() > MAX_RESOLVED_VARIANTS {
        bail!(
            "Primer {} has too many ambiguous bases: {} resolved variants exceeds limit of {}. \
             Reduce ambiguity or use a more specific primer.",
            primer,
            primer_variants.len(),
            MAX_RESOLVED_VARIANTS
        );
    }

    // Clamp mismatches to the trimmed primer length
    let mismatches = params.mismatches.min(primer.len());

    // Get all possible variants of the primers
    primer_variants = permute_sequences(primer_variants, &mismatches);

    debug!(
        "  There are {} variants of the primer",
        primer_variants.len()
    );

    Ok(primer_variants)
}

/// Given a set of primer variants, return a set of kmers from the data that contain the primers
pub(super) fn get_kmers_from_primers(
    primer_variants: &HashSet<String>,
    kmer_counts: &FilteredKmerCounts,
    min_count: &u32,
) -> Result<KmerCounts> {
    // Get the kmers that contain the primers
    let mut oligos: Vec<Oligo> = Vec::new();
    for variant in primer_variants.iter() {
        oligos.push(string_to_oligo(variant)?);
    }

    Ok(find_oligos_in_kmers(&oligos, kmer_counts, min_count))
}

/// Given a set of kmers that contain primers, filter them to retain only those with the highest counts.
/// If there are more than `max_primer_kmers`, retain only that many with the highest counts.
/// Returns a hash map of the kmers and their counts.
pub(super) fn filter_primer_kmers(matches: KmerCounts, max_primer_kmers: usize) -> KmerCounts {
    if matches.is_empty() {
        return matches;
    }

    let mut counts: Vec<u32> = matches.counts();

    counts.sort();
    counts.reverse();

    // set equal to last element
    let mut top_count_cutoff = counts.last().expect("counts is non-empty (checked above)");

    if counts.len() > max_primer_kmers {
        top_count_cutoff = &counts[max_primer_kmers - 1];
    }

    let mut matches_keep: KmerCounts = KmerCounts::new(&matches.get_k());
    for (kmer, count) in matches.iter() {
        let mut keep: bool = false;
        if count >= top_count_cutoff {
            // Add the kmer to the hash map
            matches_keep.insert(kmer, count);
            keep = true;
        }
        debug!(
            "    {}, count {}, keep {}",
            crate::kmer::kmer_to_seq(kmer, &matches.get_k()),
            count,
            keep,
        );
    }

    // Replace matches with matches_keep
    matches_keep
}

pub(super) fn get_primer_kmers(
    params: &PCRParams,
    kmer_counts: &FilteredKmerCounts,
) -> Result<(KmerCounts, KmerCounts)> {
    // Preprocess the primers to get all variants to be considered
    let forward_variants =
        preprocess_primer(params, PrimerDirection::Forward, &kmer_counts.get_k())?;
    let reverse_variants =
        preprocess_primer(params, PrimerDirection::Reverse, &kmer_counts.get_k())?;

    // Get the kmers that contain the primers
    gene_info!(
        params.gene_name,
        "Searching kmers that contain the forward primer variants"
    );
    let mut forward_primer_kmers =
        get_kmers_from_primers(&forward_variants, kmer_counts, &params.min_count)?;
    forward_primer_kmers = filter_primer_kmers(forward_primer_kmers, params.max_primer_kmers);

    gene_info!(
        params.gene_name,
        "Searching kmers that contain the reverse primer variants"
    );
    let mut reverse_primer_kmers =
        get_kmers_from_primers(&reverse_variants, kmer_counts, &params.min_count)?;
    reverse_primer_kmers = filter_primer_kmers(reverse_primer_kmers, params.max_primer_kmers);

    Ok((forward_primer_kmers, reverse_primer_kmers))
}

/// Pre-encoded primer Oligos for a single gene, used for read retention
/// during Pass 1 before the kmer table is available.
#[allow(dead_code)]
pub struct PrimerOligoSet {
    /// Gene name for attribution
    pub gene_name: String,
    /// 2-bit encoded forward primer Oligos (all variants after trim/ambiguity/mismatch)
    pub forward_oligos: Vec<Oligo>,
    /// 2-bit encoded reverse primer Oligos
    pub reverse_oligos: Vec<Oligo>,
    /// Oligo length (all Oligos in a set have the same length per direction)
    pub forward_oligo_length: usize,
    pub reverse_oligo_length: usize,
}

/// Pre-encode primer Oligos for all genes before read ingestion.
/// This runs the text-only portion of primer preprocessing (trim,
/// ambiguity resolution, mismatch permutation) and encodes each
/// variant as a 2-bit Oligo. No kmer table access.
pub fn preprocess_primer_oligos(pcr_runs: &[PCRParams], k: usize) -> Result<Vec<PrimerOligoSet>> {
    let mut result = Vec::new();

    for params in pcr_runs {
        let forward_variants = preprocess_primer(params, PrimerDirection::Forward, &k)?;
        let reverse_variants = preprocess_primer(params, PrimerDirection::Reverse, &k)?;

        let forward_oligos: Vec<Oligo> = forward_variants
            .iter()
            .map(|v| string_to_oligo(v))
            .collect::<Result<Vec<_>>>()?;

        let reverse_oligos: Vec<Oligo> = reverse_variants
            .iter()
            .map(|v| string_to_oligo(v))
            .collect::<Result<Vec<_>>>()?;

        let forward_oligo_length = forward_oligos.first().map_or(0, |o| o.length);
        let reverse_oligo_length = reverse_oligos.first().map_or(0, |o| o.length);

        result.push(PrimerOligoSet {
            gene_name: params.gene_name.clone(),
            forward_oligos,
            reverse_oligos,
            forward_oligo_length,
            reverse_oligo_length,
        });
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- string_to_oligo ---

    #[test]
    fn test_string_to_oligo_single_base() {
        // A=00, C=01, G=10, T=11
        let a = string_to_oligo("A").unwrap();
        assert_eq!(a.kmer, 0b00);
        assert_eq!(a.length, 1);

        let t = string_to_oligo("T").unwrap();
        assert_eq!(t.kmer, 0b11);
        assert_eq!(t.length, 1);
    }

    #[test]
    fn test_string_to_oligo_encoding() {
        // ACGT => 00_01_10_11
        let oligo = string_to_oligo("ACGT").unwrap();
        assert_eq!(oligo.kmer, 0b00011011);
        assert_eq!(oligo.length, 4);
    }

    #[test]
    fn test_string_to_oligo_invalid_base() {
        assert!(string_to_oligo("ACNGT").is_err());
        assert!(string_to_oligo("X").is_err());
    }

    #[test]
    fn test_string_to_oligo_empty() {
        let oligo = string_to_oligo("").unwrap();
        assert_eq!(oligo.kmer, 0);
        assert_eq!(oligo.length, 0);
    }

    // --- resolve_primer ---

    #[test]
    fn test_resolve_primer_no_ambiguity() {
        let result = resolve_primer("ACGT");
        assert_eq!(result.len(), 1);
        assert!(result.contains("ACGT"));
    }

    #[test]
    fn test_resolve_primer_single_ambiguity() {
        // R = A or G
        let result = resolve_primer("AR");
        assert_eq!(result.len(), 2);
        assert!(result.contains("AA"));
        assert!(result.contains("AG"));
    }

    #[test]
    fn test_resolve_primer_multiple_ambiguities() {
        // R = A|G, Y = C|T => 2*2 = 4 variants
        let result = resolve_primer("RY");
        assert_eq!(result.len(), 4);
        assert!(result.contains("AC"));
        assert!(result.contains("AT"));
        assert!(result.contains("GC"));
        assert!(result.contains("GT"));
    }

    #[test]
    fn test_resolve_primer_n_gives_four() {
        let result = resolve_primer("N");
        assert_eq!(result.len(), 4);
        for base in &["A", "C", "G", "T"] {
            assert!(result.contains(*base));
        }
    }

    // --- combinations ---

    #[test]
    fn test_combinations_basic() {
        assert_eq!(combinations(4, 2).len(), 6); // C(4,2) = 6
        assert_eq!(combinations(5, 0).len(), 1); // C(n,0) = 1
        assert_eq!(combinations(3, 3).len(), 1); // C(n,n) = 1
    }

    #[test]
    fn test_combinations_r_exceeds_n() {
        assert!(combinations(2, 5).is_empty());
    }

    // --- permute_sequences ---

    #[test]
    fn test_permute_zero_mismatches() {
        let mut seqs = HashSet::new();
        seqs.insert("ACG".to_string());
        let result = permute_sequences(seqs, &0);
        assert_eq!(result.len(), 1);
        assert!(result.contains("ACG"));
    }

    #[test]
    fn test_permute_one_mismatch() {
        let mut seqs = HashSet::new();
        seqs.insert("AC".to_string());
        let result = permute_sequences(seqs, &1);
        // Position 0: AC, TC, CC, GC (4); Position 1: AA, AT, AC, AG (4)
        // Union = 7 unique (AC counted once)
        assert_eq!(result.len(), 7);
        assert!(result.contains("AC")); // original
        assert!(result.contains("TC")); // pos 0 mutated
        assert!(result.contains("AG")); // pos 1 mutated
    }

    // --- find_oligos_in_kmers ---

    /// Build a KmerCounts table from a short sequence for testing primer matching.
    fn build_kmer_counts(seq: &str, k: usize) -> KmerCounts {
        let mut kmer_counts = KmerCounts::new(&k);
        let reads = crate::kmer::seq_to_reads(seq).unwrap();
        kmer_counts.ingest_reads(&reads).unwrap();
        kmer_counts
    }

    #[test]
    fn test_find_oligos_exact_match_forward() {
        // Sequence: ACGTACGT, k=5
        // Oligo "ACG" (length 3) should match kmers starting with ACG
        let k = 5;
        let kmer_counts = build_kmer_counts("ACGTACGT", k);
        let filtered = kmer_counts.filtered_view(1);

        let oligo = string_to_oligo("ACG").unwrap();
        let result = find_oligos_in_kmers(&[oligo], &filtered, &1);

        // Kmers in ACGTACGT: ACGTA, CGTAC, GTACG, TACGT
        // Forward match (starts with ACG): ACGTA
        // RC match: reverse complement of a kmer ending with CGT (=ACG rc)
        assert!(
            result.len() > 0,
            "Should find at least one kmer matching oligo ACG"
        );
    }

    #[test]
    fn test_find_oligos_no_match() {
        // Sequence of all A's — oligo "GGG" should not match
        let k = 5;
        let kmer_counts = build_kmer_counts("AAAAAAAAAA", k);
        let filtered = kmer_counts.filtered_view(1);

        let oligo = string_to_oligo("GGG").unwrap();
        let result = find_oligos_in_kmers(&[oligo], &filtered, &1);

        assert_eq!(result.len(), 0, "Should find no matching kmers");
    }

    #[test]
    fn test_find_oligos_min_count_filter() {
        // Use a non-palindromic sequence so canonical kmer counts stay at 1.
        // AACCCAACC has k=5 kmers: AACCC, ACCCA, CCCAA, CCAAC, CAACC
        // None of these are their own reverse complement, so each count=1.
        let k = 5;
        let kmer_counts = build_kmer_counts("AACCCAACC", k);
        let filtered = kmer_counts.filtered_view(1);

        let oligo = string_to_oligo("AAC").unwrap();
        let result = find_oligos_in_kmers(&[oligo], &filtered, &2);

        assert_eq!(
            result.len(),
            0,
            "Single-copy kmers should not pass min_count=2"
        );
    }

    #[test]
    fn test_find_oligos_rc_match() {
        // Oligo "AAA" (length 3). The RC of "AAA" is "TTT".
        // Sequence "TTTTTTT" with k=5: all kmers are TTTTT.
        // RC of TTTTT = AAAAA. The high bits of AAAAA start with AAA → match.
        // The matched kmer should be stored as the RC (AAAAA).
        let k = 5;
        let kmer_counts = build_kmer_counts("TTTTTTT", k);
        let filtered = kmer_counts.filtered_view(1);

        let oligo = string_to_oligo("AAA").unwrap();
        let result = find_oligos_in_kmers(&[oligo], &filtered, &1);

        assert!(
            result.len() > 0,
            "Should find RC match: oligo AAA should match kmers in all-T sequence"
        );
        // The matched kmer should be stored in the forward-oligo orientation (AAAAA)
        let matched_kmer = result.iter().next().unwrap().0;
        let matched_seq = crate::kmer::kmer_to_seq(matched_kmer, &k);
        assert_eq!(matched_seq, "AAAAA");
    }

    #[test]
    fn test_find_oligos_oligo_equals_k() {
        // Edge case: oligo_length == k (shift = 0, mask covers all bits)
        let k = 5;
        let kmer_counts = build_kmer_counts("ACGTACGT", k);
        let filtered = kmer_counts.filtered_view(1);

        let oligo = string_to_oligo("ACGTA").unwrap();
        assert_eq!(oligo.length, k);
        let result = find_oligos_in_kmers(&[oligo], &filtered, &1);

        assert!(
            result.len() > 0,
            "Full-length oligo (length == k) should match the exact kmer"
        );
    }

    // --- filter_primer_kmers ---

    #[test]
    fn test_filter_primer_kmers_empty() {
        let empty = KmerCounts::new(&5);
        let result = filter_primer_kmers(empty, 10);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_filter_primer_kmers_cap() {
        // Insert 5 kmers with different counts, cap at 3
        let k = 5;
        let mut counts = KmerCounts::new(&k);
        for (i, seq) in ["AAAAA", "AAAAC", "AAACG", "AACGT", "ACGTA"]
            .iter()
            .enumerate()
        {
            let kmer = crate::kmer::seq_to_kmer(seq).unwrap();
            let count = (i + 1) as u32; // counts: 1, 2, 3, 4, 5
            counts.insert(&kmer, &count);
        }
        let result = filter_primer_kmers(counts, 3);
        // Top 3 by count: 5, 4, 3 → 3 kmers retained
        assert_eq!(result.len(), 3);
    }

    // --- is_valid_nucleotide ---

    #[test]
    fn test_valid_nucleotides() {
        for c in "ACGTRYWSKMBDHVN".chars() {
            assert!(is_valid_nucleotide(c), "{} should be valid", c);
        }
        assert!(!is_valid_nucleotide('X'));
        assert!(!is_valid_nucleotide('a')); // lowercase not accepted
    }
}
