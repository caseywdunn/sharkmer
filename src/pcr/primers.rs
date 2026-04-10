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

/// Ingest pcr parameters, and return primer variants grouped by mismatch
/// level. Index 0 contains the exact (ambiguity-resolved) variants, index 1
/// the 1-mismatch variants, etc. Each level contains only the NEW variants
/// introduced at that mismatch count (not cumulative).
pub(super) fn preprocess_primer_by_mismatch(
    params: &PCRParams,
    dir: PrimerDirection,
    k: &usize,
) -> Result<Vec<HashSet<String>>> {
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
    let base_variants = resolve_primer(&primer);

    const MAX_RESOLVED_VARIANTS: usize = 10_000;
    if base_variants.len() > MAX_RESOLVED_VARIANTS {
        bail!(
            "Primer {} has too many ambiguous bases: {} resolved variants exceeds limit of {}. \
             Reduce ambiguity or use a more specific primer.",
            primer,
            base_variants.len(),
            MAX_RESOLVED_VARIANTS
        );
    }

    // Clamp mismatches to the trimmed primer length
    let mismatches = params.mismatches.min(primer.len());

    // Build variants level by level: 0 mismatches, 1 mismatch, ..., n mismatches.
    // Each level contains only the NEW variants not seen at lower levels.
    let mut levels: Vec<HashSet<String>> = Vec::with_capacity(mismatches + 1);
    let mut seen: HashSet<String> = HashSet::new();

    // Level 0: exact matches (ambiguity-resolved only)
    levels.push(base_variants.clone());
    seen.extend(base_variants.iter().cloned());

    // Levels 1..=mismatches: permutations at each distance
    for _m in 1..=mismatches {
        let all_up_to_m = permute_sequences(seen.clone(), &1);
        let new_at_m: HashSet<String> = all_up_to_m.difference(&seen).cloned().collect();
        seen.extend(new_at_m.iter().cloned());
        levels.push(new_at_m);
    }

    let total: usize = levels.iter().map(|l| l.len()).sum();
    debug!(
        "  There are {} variants of the primer across {} mismatch levels",
        total,
        levels.len()
    );

    Ok(levels)
}

/// Convenience wrapper: returns all primer variants as a flat set (all mismatch
/// levels combined). Used by `preprocess_primer_oligos()` for read retention,
/// where mismatch-level ordering is not needed.
pub(super) fn preprocess_primer(
    params: &PCRParams,
    dir: PrimerDirection,
    k: &usize,
) -> Result<HashSet<String>> {
    let levels = preprocess_primer_by_mismatch(params, dir, k)?;
    Ok(levels.into_iter().flatten().collect())
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

/// Hard-cap a KmerCounts to at most `max_primer_kmers` entries.
/// Sorts by count descending, then by kmer value ascending (deterministic
/// tiebreaker). Keeps exactly `min(len, max_primer_kmers)` entries.
#[cfg(test)]
pub(super) fn filter_primer_kmers(matches: KmerCounts, max_primer_kmers: usize) -> KmerCounts {
    if matches.len() <= max_primer_kmers {
        return matches;
    }

    // Collect (kmer, count) pairs and sort: count DESC, kmer ASC
    let mut entries: Vec<(u64, u32)> = matches.iter().map(|(k, c)| (*k, *c)).collect();
    entries.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    let mut result = KmerCounts::new(&matches.get_k());
    for &(kmer, count) in entries.iter().take(max_primer_kmers) {
        result.insert(&kmer, &count);
    }

    debug!(
        "    Hard-capped primer kmers from {} to {} (cutoff count >= {})",
        entries.len(),
        max_primer_kmers,
        entries[max_primer_kmers - 1].1
    );

    result
}

/// Discover primer kmers round by round, one mismatch level at a time.
/// Lower mismatch levels fill the cap first; higher levels only contribute
/// if there is remaining capacity. Within each round, kmers are prioritized
/// by count (descending) then kmer value (ascending, deterministic).
fn discover_primer_kmers_by_round(
    variant_levels: &[HashSet<String>],
    kmer_counts: &FilteredKmerCounts,
    min_count: &u32,
    max_primer_kmers: usize,
    gene_name: &str,
) -> Result<KmerCounts> {
    let mut result = KmerCounts::new(&kmer_counts.get_k());

    for (mismatch_level, variants) in variant_levels.iter().enumerate() {
        if result.len() >= max_primer_kmers {
            break;
        }
        if variants.is_empty() {
            continue;
        }

        // Find kmers matching this mismatch level's variants
        let round_kmers = get_kmers_from_primers(variants, kmer_counts, min_count)?;

        // Filter out kmers already found at lower mismatch levels
        let mut new_entries: Vec<(u64, u32)> = Vec::new();
        for (kmer, count) in round_kmers.iter() {
            if result.get(kmer).is_none() {
                new_entries.push((*kmer, *count));
            }
        }

        if new_entries.is_empty() {
            continue;
        }

        // Sort new entries: count DESC, kmer ASC (deterministic tiebreaker)
        new_entries.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

        // Take only as many as we have room for
        let remaining = max_primer_kmers - result.len();
        let take = new_entries.len().min(remaining);

        for &(kmer, count) in new_entries.iter().take(take) {
            result.insert(&kmer, &count);
        }

        gene_info!(
            gene_name,
            "Mismatch level {}: {} new primer kmers ({} total, cap {})",
            mismatch_level,
            take,
            result.len(),
            max_primer_kmers
        );

        if new_entries.len() > take {
            gene_info!(
                gene_name,
                "Mismatch level {}: dropped {} kmers at cap",
                mismatch_level,
                new_entries.len() - take
            );
        }
    }

    Ok(result)
}

pub(super) fn get_primer_kmers(
    params: &PCRParams,
    kmer_counts: &FilteredKmerCounts,
) -> Result<(KmerCounts, KmerCounts)> {
    // Preprocess the primers into mismatch-level groups
    let forward_levels =
        preprocess_primer_by_mismatch(params, PrimerDirection::Forward, &kmer_counts.get_k())?;
    let reverse_levels =
        preprocess_primer_by_mismatch(params, PrimerDirection::Reverse, &kmer_counts.get_k())?;

    // Discover primer kmers round by round, filling the cap from lowest mismatch first
    gene_info!(
        params.gene_name,
        "Searching kmers that contain the forward primer variants"
    );
    let forward_primer_kmers = discover_primer_kmers_by_round(
        &forward_levels,
        kmer_counts,
        &params.min_count,
        params.max_primer_kmers,
        &params.gene_name,
    )?;

    gene_info!(
        params.gene_name,
        "Searching kmers that contain the reverse primer variants"
    );
    let reverse_primer_kmers = discover_primer_kmers_by_round(
        &reverse_levels,
        kmer_counts,
        &params.min_count,
        params.max_primer_kmers,
        &params.gene_name,
    )?;

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
            !result.is_empty(),
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
            !result.is_empty(),
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
            !result.is_empty(),
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

    #[test]
    fn test_filter_primer_kmers_hard_cap_with_ties() {
        // Insert 5 kmers all with the SAME count, cap at 3.
        // Old implementation would keep all 5 (leaky cap).
        // New implementation must keep exactly 3, breaking ties by kmer value.
        let k = 5;
        let mut counts = KmerCounts::new(&k);
        for seq in ["AAAAA", "AAAAC", "AAACG", "AACGT", "ACGTA"] {
            let kmer = crate::kmer::seq_to_kmer(seq).unwrap();
            let count = 2u32; // all same count
            counts.insert(&kmer, &count);
        }
        let result = filter_primer_kmers(counts, 3);
        assert_eq!(
            result.len(),
            3,
            "Hard cap must enforce exactly 3 even when all counts are tied"
        );
    }

    #[test]
    fn test_filter_primer_kmers_under_cap() {
        // Insert 2 kmers, cap at 5 — should keep all
        let k = 5;
        let mut counts = KmerCounts::new(&k);
        for seq in ["AAAAA", "AAAAC"] {
            let kmer = crate::kmer::seq_to_kmer(seq).unwrap();
            counts.insert(&kmer, &3);
        }
        let result = filter_primer_kmers(counts, 5);
        assert_eq!(result.len(), 2, "Should keep all when under cap");
    }

    #[test]
    fn test_preprocess_primer_by_mismatch_levels() {
        // Test that mismatch levels are non-overlapping and cumulative union
        // equals the flat preprocess_primer output
        use super::super::*;
        use super::*;
        let params = PCRParams {
            gene_name: "test".to_string(),
            forward_seq: "ACGTACGT".to_string(),
            reverse_seq: "TGCATGCA".to_string(),
            min_length: 100,
            max_length: 500,
            min_count: 2,
            mismatches: 2,
            trim: 8,
            dedup_edit_threshold: DEFAULT_DEDUP_EDIT_THRESHOLD,
            max_primer_kmers: DEFAULT_MAX_NUM_PRIMER_KMERS,
            max_dfs_states: DEFAULT_MAX_DFS_STATES,
            max_paths_per_pair: DEFAULT_MAX_PATHS_PER_PAIR,
            max_node_visits: DEFAULT_MAX_NODE_VISITS,
            max_seed_nodes: DEFAULT_MAX_SEED_NODES,
            high_coverage_ratio: DEFAULT_HIGH_COVERAGE_RATIO,
            tip_coverage_fraction: DEFAULT_TIP_COVERAGE_FRACTION,
            stopping_criteria: StoppingCriteria::FirstProduct,
            min_component_budget: DEFAULT_MIN_COMPONENT_BUDGET,
            citation: String::new(),
            notes: String::new(),
            expected_length: None,
            source: String::new(),
        };
        let k = 8usize;
        let levels = preprocess_primer_by_mismatch(&params, PrimerDirection::Forward, &k).unwrap();
        assert_eq!(levels.len(), 3); // 0, 1, 2 mismatches

        // Level 0 should contain exact match
        assert!(
            levels[0].contains("ACGTACGT"),
            "Level 0 should contain the exact primer"
        );

        // Levels should be non-overlapping
        for i in 0..levels.len() {
            for j in (i + 1)..levels.len() {
                let overlap: HashSet<_> = levels[i].intersection(&levels[j]).collect();
                assert!(
                    overlap.is_empty(),
                    "Levels {} and {} should not overlap",
                    i,
                    j
                );
            }
        }

        // Union should match flat preprocess_primer
        let flat = preprocess_primer(&params, PrimerDirection::Forward, &k).unwrap();
        let union: HashSet<String> = levels.into_iter().flatten().collect();
        assert_eq!(union, flat, "Union of levels should match flat output");
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
