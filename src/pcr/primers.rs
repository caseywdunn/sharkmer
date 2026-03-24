// pcr/primers.rs — primer preprocessing functions

use anyhow::{bail, Result};
use log::debug;
use std::collections::HashSet;

use crate::kmer::{FilteredKmerCounts, KmerCounts};

use super::{Oligo, PCRParams, PrimerDirection, MAX_NUM_PRIMER_KMERS};

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
    let nucleotides = ["A", "T", "C", "G"];

    if current == positions.len() {
        unique_sequences.insert(seq.to_string());
        return;
    }

    let pos = positions[current];
    for &nucleotide in nucleotides.iter() {
        let mut new_seq = seq.chars().collect::<Vec<_>>();
        new_seq[pos] = nucleotide
            .chars()
            .next()
            .expect("nucleotide literal is non-empty");
        generate_recursive_permutations(
            &new_seq.iter().collect::<String>(),
            positions,
            current + 1,
            unique_sequences,
        );
    }
}

pub(super) fn reverse_complement(seq: &str) -> Result<String> {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => Ok('T'),
            'T' => Ok('A'),
            'G' => Ok('C'),
            'C' => Ok('G'),
            'Y' => Ok('R'), // C or T
            'R' => Ok('Y'), // A or G
            'S' => Ok('S'), // C or G
            'W' => Ok('W'), // A or T
            'K' => Ok('M'), // G or T
            'M' => Ok('K'), // A or C
            'B' => Ok('V'), // C or G or T
            'V' => Ok('B'), // A or C or G
            'D' => Ok('H'), // A or G or T
            'H' => Ok('D'), // A or C or T
            'N' => Ok('N'), // A or C or G or T
            _ => bail!("Invalid nucleotide: {}", c),
        })
        .collect()
}

// Find the kmers that contain the oligos.
// If direction is forward, match the oligo at the start of the kmer.
// If direction is reverse, match the oligo at the end of the kmer.
fn find_oligos_in_kmers(
    oligos: &[Oligo],
    kmers: &FilteredKmerCounts,
    dir: &PrimerDirection,
    min_count: &u64,
) -> KmerCounts {
    // Assume all oligos have the same length
    let oligo_length = oligos[0].length;

    // Create hash set of oligos
    let oligo_set: HashSet<u64> = match dir {
        PrimerDirection::Forward => {
            // Rotate the oligo kmers to the start of the kmer
            oligos
                .iter()
                .map(|oligo| oligo.kmer << (2 * (kmers.get_k() - oligo_length)))
                .collect()
        }
        PrimerDirection::Reverse => {
            // Just add the oligo kmers
            oligos.iter().map(|oligo| oligo.kmer).collect()
        }
    };

    // Create mask for kmers so they can be compared to the oligo subsequence
    // If dir is forward, mask is set to 1 starting at the first position of the kmer for 2*oligo_length
    // If dir is reverse, mask is set to 1 for the 2*oligo_length least significant bits
    let mut mask: u64 = 0;
    match dir {
        PrimerDirection::Forward => {
            for _i in 0..(2 * oligo_length) {
                mask = (mask << 1) | 1;
            }
            mask <<= 2 * kmers.get_k() - 2 * oligo_length;
        }
        PrimerDirection::Reverse => {
            for _i in 0..(2 * oligo_length) {
                mask = (mask << 1) | 1;
            }
        }
    }

    let mut kmers_match: KmerCounts = KmerCounts::new(&kmers.get_k());

    // Build the reverse complement mask by mirroring the mask bits
    let rc_mask = {
        let shift = 2 * kmers.get_k() - 2 * oligo_length;
        match dir {
            // Forward mask covers high bits; RC of a forward-matching kmer has primer in low bits
            PrimerDirection::Forward => (1u64 << (2 * oligo_length)) - 1,
            // Reverse mask covers low bits; RC of a reverse-matching kmer has primer in high bits
            PrimerDirection::Reverse => ((1u64 << (2 * oligo_length)) - 1) << shift,
        }
    };

    // Build the set of RC oligos (reverse complement of each oligo, shifted to match rc_mask)
    let rc_oligo_set: HashSet<u64> = oligos
        .iter()
        .map(|oligo| {
            let rc = crate::kmer::revcomp_kmer(&oligo.kmer, &oligo_length);
            match dir {
                PrimerDirection::Forward => rc,
                PrimerDirection::Reverse => {
                    let shift = 2 * kmers.get_k() - 2 * oligo_length;
                    rc << shift
                }
            }
        })
        .collect();

    for (kmer, count) in kmers.iter() {
        if count >= min_count {
            if oligo_set.contains(&(kmer & mask)) {
                kmers_match.insert(kmer, count);
            } else if rc_oligo_set.contains(&(kmer & rc_mask)) {
                // The primer matches the reverse complement orientation of this kmer,
                // so insert the revcomp form for correct graph construction
                let rc_kmer = crate::kmer::revcomp_kmer(kmer, &kmers.get_k());
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

    if dir == PrimerDirection::Reverse {
        // Replace the reverse variants with their reverse complements
        let mut primer_variants_revcomp = HashSet::new();
        for variant in primer_variants.iter() {
            primer_variants_revcomp.insert(reverse_complement(variant)?);
        }
        primer_variants = primer_variants_revcomp;
    }

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
    dir: PrimerDirection,
    min_count: &u64,
) -> Result<KmerCounts> {
    // Get the kmers that contain the primers
    let mut oligos: Vec<Oligo> = Vec::new();
    for variant in primer_variants.iter() {
        oligos.push(string_to_oligo(variant)?);
    }

    Ok(find_oligos_in_kmers(&oligos, kmer_counts, &dir, min_count))
}

/// Given a set of kmers that contain primers, filter them to retain only those with the highest counts
/// If there are more than MAX_NUM_PRIMER_KMERS, retain only the MAX_NUM_PRIMER_KMERS with the highest counts
/// If there are less than MAX_NUM_PRIMER_KMERS, retain all of them
/// Returns a hash map of the kmers and their counts
pub(super) fn filter_primer_kmers(matches: KmerCounts) -> KmerCounts {
    if matches.is_empty() {
        return matches;
    }

    let mut counts: Vec<u64> = matches.counts();

    counts.sort();
    counts.reverse();

    // set equal to last element
    let mut top_count_cutoff = counts.last().expect("counts is non-empty (checked above)");

    if counts.len() > MAX_NUM_PRIMER_KMERS {
        top_count_cutoff = &counts[MAX_NUM_PRIMER_KMERS - 1];
    }

    // If there are less than MAX_NUM_PRIMER_KMERS matches, this is the lowest count
    // If there are more than MAX_NUM_PRIMER_KMERS matches, get the count of the
    // MAX_NUM_PRIMER_KMERS highest count matches and use that as the cutoff

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
    let mut forward_primer_kmers = get_kmers_from_primers(
        &forward_variants,
        kmer_counts,
        PrimerDirection::Forward,
        &params.min_count,
    )?;
    forward_primer_kmers = filter_primer_kmers(forward_primer_kmers);

    gene_info!(
        params.gene_name,
        "Searching kmers that contain the reverse primer variants"
    );
    let mut reverse_primer_kmers = get_kmers_from_primers(
        &reverse_variants,
        kmer_counts,
        PrimerDirection::Reverse,
        &params.min_count,
    )?;
    reverse_primer_kmers = filter_primer_kmers(reverse_primer_kmers);

    Ok((forward_primer_kmers, reverse_primer_kmers))
}
