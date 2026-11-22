// pcr/primers.rs — primer preprocessing functions

use super::memory::buffers::{BudgetSet, BudgetVec};
use crate::kmer::budget::ResourceBudget;
use anyhow::{Context, Result, bail, ensure};
use log::debug;
#[cfg(test)]
use std::collections::HashSet;
use std::sync::Arc;

#[cfg(test)]
use crate::kmer::KmerCounts;
use crate::kmer::{FilteredKmerCounts, PrimerKmers};

#[cfg(test)]
use super::Oligo;
use super::{PCRParams, PrimerDirection};

const MAX_RESOLVED_VARIANTS: usize = 10_000;
const MAX_RETAINED_VARIANTS: usize = 1_000_000;
const MAX_INTERMEDIATE_VARIANTS: usize = 4_000_000;
const MAX_GENERATION_WORK: usize = 1_000_000;

struct PrimerExpansionPlan<'a> {
    primer: &'a str,
    mismatches: usize,
    #[cfg_attr(not(test), allow(dead_code))]
    retained_variants: usize,
    #[cfg_attr(not(test), allow(dead_code))]
    generation_work: usize,
}

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
#[cfg(test)]
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
#[cfg(test)]
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

fn ambiguity_width(nucleotide: char) -> Result<usize> {
    match nucleotide {
        'A' | 'C' | 'G' | 'T' => Ok(1),
        'R' | 'Y' | 'S' | 'W' | 'K' | 'M' => Ok(2),
        'B' | 'D' | 'H' | 'V' => Ok(3),
        'N' => Ok(4),
        _ => bail!("Invalid nucleotide {} in primer", nucleotide),
    }
}

fn checked_product(left: usize, right: usize, primer: &str, limit_name: &str) -> Result<usize> {
    left.checked_mul(right).ok_or_else(|| {
        anyhow::anyhow!(
            "Primer {} exceeds the {} limit before expansion",
            primer,
            limit_name
        )
    })
}

fn checked_sum(left: usize, right: usize, primer: &str, limit_name: &str) -> Result<usize> {
    left.checked_add(right).ok_or_else(|| {
        anyhow::anyhow!(
            "Primer {} exceeds the {} limit before expansion",
            primer,
            limit_name
        )
    })
}

fn sum_with_limit(values: &[usize], primer: &str, limit_name: &str, limit: usize) -> Result<usize> {
    let mut total = 0usize;
    for value in values {
        total = checked_sum(total, *value, primer, limit_name)?;
        if total > limit {
            bail!(
                "Primer {} requires {} retained variants, exceeding the limit of {}. Reduce ambiguity or mismatches.",
                primer,
                total,
                limit
            );
        }
    }
    Ok(total)
}

fn intermediate_peak(
    exact_variants: usize,
    previous_seen: usize,
    current_seen: usize,
    primer: &str,
) -> Result<usize> {
    let permutation_peak = checked_sum(
        checked_product(previous_seen, 3, primer, "intermediate-variant")?,
        current_seen,
        primer,
        "intermediate-variant",
    )?;
    let post_extend_peak = checked_product(current_seen, 3, primer, "intermediate-variant")?;
    checked_sum(
        exact_variants,
        permutation_peak.max(post_extend_peak),
        primer,
        "intermediate-variant",
    )
}

fn ensure_intermediate_limit(peak: usize, primer: &str) -> Result<()> {
    if peak > MAX_INTERMEDIATE_VARIANTS {
        bail!(
            "Primer {} requires {} intermediate variants, exceeding the limit of {}. Reduce ambiguity or mismatches.",
            primer,
            peak,
            MAX_INTERMEDIATE_VARIANTS
        );
    }
    Ok(())
}

fn effective_trim(params: &PCRParams, kmer_length: usize) -> Result<usize> {
    let max_trim = kmer_length.checked_sub(1).ok_or_else(|| {
        anyhow::anyhow!(
            "PCR primer processing requires k-mer length of at least 2, got {}",
            kmer_length
        )
    })?;
    let trim = params.trim.min(max_trim);
    ensure!(
        trim > 0,
        "Primer {} has effective trim length 0 for k-mer length {}. Set trim to at least 1 and use k-mer length at least 2.",
        params.gene_name,
        kmer_length
    );
    Ok(trim)
}

fn plan_primer_expansion(
    params: &PCRParams,
    dir: PrimerDirection,
    kmer_length: usize,
) -> Result<PrimerExpansionPlan<'_>> {
    ensure!(
        kmer_length < 32,
        "PCR primer processing requires k-mer length below 32, got {}",
        kmer_length
    );
    let trim = effective_trim(params, kmer_length)?;
    let source_primer = if dir == PrimerDirection::Reverse {
        &params.reverse_seq
    } else {
        &params.forward_seq
    };
    let primer = source_primer
        .get(source_primer.len().saturating_sub(trim)..)
        .context("Primer trimming crosses an invalid nucleotide boundary")?;
    ensure!(
        !primer.is_empty(),
        "Primer {} has no retained bases after trimming",
        params.gene_name
    );

    let mut exact_variants = 1usize;
    for nucleotide in primer.chars() {
        exact_variants = checked_product(
            exact_variants,
            ambiguity_width(nucleotide)?,
            primer,
            "resolved-variant",
        )?;
        if exact_variants > MAX_RESOLVED_VARIANTS {
            bail!(
                "Primer {} has too many ambiguous bases: more than {} resolved variants exceeds limit of {}. Reduce ambiguity or use a more specific primer.",
                primer,
                MAX_RESOLVED_VARIANTS,
                MAX_RESOLVED_VARIANTS
            );
        }
    }

    let mismatches = params.mismatches.min(primer.len());
    let mut distance_counts = [0usize; 32];
    distance_counts[0] = 1;
    for nucleotide in primer.chars() {
        let matching_bases = ambiguity_width(nucleotide)?;
        let mismatching_bases = 4 - matching_bases;
        let mut next_counts = [0usize; 32];
        for mismatch_count in 0..=mismatches {
            if distance_counts[mismatch_count] == 0 {
                continue;
            }
            let matching_count = checked_product(
                distance_counts[mismatch_count],
                matching_bases,
                primer,
                "retained-variant",
            )?;
            next_counts[mismatch_count] = checked_sum(
                next_counts[mismatch_count],
                matching_count,
                primer,
                "retained-variant",
            )?;
            if mismatch_count < mismatches {
                let mismatching_count = checked_product(
                    distance_counts[mismatch_count],
                    mismatching_bases,
                    primer,
                    "retained-variant",
                )?;
                next_counts[mismatch_count + 1] = checked_sum(
                    next_counts[mismatch_count + 1],
                    mismatching_count,
                    primer,
                    "retained-variant",
                )?;
            }
        }
        sum_with_limit(
            &next_counts[..=mismatches],
            primer,
            "retained-variant",
            MAX_RETAINED_VARIANTS,
        )?;
        distance_counts = next_counts;
    }

    let retained_variants = sum_with_limit(
        &distance_counts[..=mismatches],
        primer,
        "retained-variant",
        MAX_RETAINED_VARIANTS,
    )?;
    let mut previous_seen = distance_counts[0];
    let mut generation_work = 0usize;
    let initial_peak = checked_product(exact_variants, 3, primer, "intermediate-variant")?;
    ensure_intermediate_limit(initial_peak, primer)?;
    for mismatch_count in 1..=mismatches {
        let current_seen = sum_with_limit(
            &distance_counts[..=mismatch_count],
            primer,
            "retained-variant",
            MAX_RETAINED_VARIANTS,
        )?;
        let work = checked_product(
            checked_product(previous_seen, primer.len(), primer, "generation-work")?,
            4,
            primer,
            "generation-work",
        )?;
        generation_work = checked_sum(generation_work, work, primer, "generation-work")?;
        if generation_work > MAX_GENERATION_WORK {
            bail!(
                "Primer {} requires {} generated mismatch candidates through round {}, exceeding the limit of {}. Reduce ambiguity or mismatches.",
                primer,
                generation_work,
                mismatch_count,
                MAX_GENERATION_WORK
            );
        }
        let peak = intermediate_peak(exact_variants, previous_seen, current_seen, primer)?;
        ensure_intermediate_limit(peak, primer)?;
        previous_seen = current_seen;
    }

    Ok(PrimerExpansionPlan {
        primer,
        mismatches,
        retained_variants,
        generation_work,
    })
}

pub(super) fn validate_primer_expansion(params: &PCRParams, kmer_length: usize) -> Result<()> {
    plan_primer_expansion(params, PrimerDirection::Forward, kmer_length)?;
    plan_primer_expansion(params, PrimerDirection::Reverse, kmer_length)?;
    Ok(())
}

/// Given a set of sequences, return a set of all sequences that differ
/// from each original sequence at up to r positions. Includes the original
/// sequences.
#[cfg(test)]
pub(super) fn permute_sequences(sequences: HashSet<String>, r: &usize) -> HashSet<String> {
    let mut permutations = HashSet::new();

    for seq in &sequences {
        for positions in combinations(seq.len(), *r) {
            generate_recursive_permutations(seq, &positions, 0, &mut permutations);
        }
    }

    permutations
}

#[cfg(test)]
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

#[cfg(test)]
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
#[cfg(test)]
fn find_oligos_in_kmers(
    oligos: &[Oligo],
    kmers: &FilteredKmerCounts,
    min_count: &u32,
) -> Result<PrimerKmers> {
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
        oligo_length > 0 && oligo_length < k,
        "oligo length {} out of range for k={} (must be 1..k-1); trim must be < k",
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

    let mut kmers_match = PrimerKmers::new(&k)?;

    // RC mask covers low bits (reverse complement of a start-matching kmer has primer in low bits)
    let rc_mask = (1u64 << (2 * oligo_length)) - 1;

    // Build the set of RC oligos (reverse complement of each oligo, unshifted)
    let rc_oligo_set: HashSet<u64> = oligos
        .iter()
        .map(|oligo| crate::kmer::revcomp_kmer(&oligo.kmer, &oligo_length))
        .collect();

    // The released primer scan used all underlying entries and its target's
    // min_count; the global graph-view threshold must not change seed choice.
    kmers.visit_unfiltered(&mut |kmer, count| {
        if count >= u64::from(*min_count) {
            if oligo_set.contains(&(kmer & mask)) {
                kmers_match.insert(&kmer, &count)?;
            } else if rc_oligo_set.contains(&(kmer & rc_mask)) {
                // The primer matches the reverse complement orientation of this kmer,
                // so insert the revcomp form for correct graph construction
                let rc_kmer = crate::kmer::revcomp_kmer(&kmer, &k);
                kmers_match.insert(&rc_kmer, &count)?;
            }
        }
        Ok(())
    })?;

    Ok(kmers_match)
}

// Individual steps in the PCR process

/// Ingest pcr parameters, and return primer variants grouped by mismatch
/// level. Index 0 contains the exact (ambiguity-resolved) variants, index 1
/// the 1-mismatch variants, etc. Each level contains only the NEW variants
/// introduced at that mismatch count (not cumulative).
#[cfg(test)]
pub(super) fn preprocess_primer_by_mismatch(
    params: &PCRParams,
    dir: PrimerDirection,
    k: &usize,
) -> Result<Vec<HashSet<String>>> {
    let source_primer = if dir == PrimerDirection::Reverse {
        &params.reverse_seq
    } else {
        &params.forward_seq
    };
    let plan = plan_primer_expansion(params, dir, *k)?;
    // The seed node sub_kmer is the PREFIX of the primer k-mer (kmer >> 2),
    // which drops the last base of the k-mer.  When the trimmed primer is
    // exactly k bases it IS the k-mer, so its 3' terminal base is dropped
    // from the seed node and the assembled amplicon.  Clamping to k-1
    // ensures all trimmed primer bases are captured in the seed node.
    if params.trim >= *k {
        gene_warn!(
            params.gene_name,
            "Trim length ({}) must be less than k ({}); adjusting trim to k-1 = {}",
            params.trim,
            k,
            k - 1
        );
    }
    if source_primer.len() > plan.primer.len() {
        gene_info!(
            params.gene_name,
            "Trimming the primer to {} so that it is within the trim length of {}.",
            plan.primer,
            plan.primer.len()
        );
    }
    let primer = plan.primer;

    let base_variants = resolve_primer(primer);
    let mismatches = plan.mismatches;

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

/// Convenience wrapper: returns all primer variants as a flat set (all
/// mismatch levels combined). Used by tests.
#[cfg(test)]
pub(super) fn preprocess_primer(
    params: &PCRParams,
    dir: PrimerDirection,
    k: &usize,
) -> Result<HashSet<String>> {
    let levels = preprocess_primer_by_mismatch(params, dir, k)?;
    Ok(levels.into_iter().flatten().collect())
}

/// Given a set of primer variants, return a set of kmers from the data that contain the primers
#[cfg(test)]
pub(super) fn get_kmers_from_primers(
    primer_variants: &HashSet<String>,
    kmer_counts: &FilteredKmerCounts,
    min_count: &u32,
) -> Result<PrimerKmers> {
    // Get the kmers that contain the primers
    let mut oligos: Vec<Oligo> = Vec::new();
    for variant in primer_variants.iter() {
        oligos.push(string_to_oligo(variant)?);
    }

    find_oligos_in_kmers(&oligos, kmer_counts, min_count)
}

/// Hard-cap directional primer candidates to at most `max_primer_kmers` entries.
/// Sorts by count descending, then by kmer value ascending (deterministic
/// tiebreaker). Keeps exactly `min(len, max_primer_kmers)` entries.
#[cfg(test)]
pub(super) fn filter_primer_kmers(
    matches: PrimerKmers,
    max_primer_kmers: usize,
) -> Result<PrimerKmers> {
    if matches.len() <= max_primer_kmers {
        return Ok(matches);
    }

    // Collect (kmer, count) pairs and sort: count DESC, kmer ASC
    let mut entries: Vec<(u64, u64)> = matches.iter().map(|(k, c)| (*k, *c)).collect();
    entries.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    let mut result = PrimerKmers::new(&matches.get_k())?;
    for &(kmer, count) in entries.iter().take(max_primer_kmers) {
        result.insert(&kmer, &count)?;
    }

    debug!(
        "    Hard-capped primer kmers from {} to {} (cutoff count >= {})",
        entries.len(),
        max_primer_kmers,
        entries[max_primer_kmers - 1].1
    );

    Ok(result)
}

/// Discover primer kmers round by round, one mismatch level at a time.
/// Lower mismatch levels fill the cap first; higher levels only contribute
/// if there is remaining capacity. Within each round, kmers are prioritized
/// by count (descending) then kmer value (ascending, deterministic).
#[cfg(test)]
fn discover_primer_kmers_by_round(
    variant_levels: &[HashSet<String>],
    kmer_counts: &FilteredKmerCounts,
    min_count: &u32,
    max_primer_kmers: usize,
    gene_name: &str,
) -> Result<PrimerKmers> {
    let mut result = PrimerKmers::new(&kmer_counts.get_k())?;

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
        let mut new_entries: Vec<(u64, u64)> = Vec::new();
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
            result.insert(&kmer, &count)?;
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

#[cfg(test)]
fn get_primer_kmers_reference(
    params: &PCRParams,
    kmer_counts: &FilteredKmerCounts,
) -> Result<(PrimerKmers, PrimerKmers)> {
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

struct EncodedPrimerLevels {
    length: usize,
    levels: BudgetVec<BudgetSet<u64>>,
}

fn allowed_bases(nucleotide: u8) -> Result<u8> {
    match nucleotide {
        b'A' => Ok(0b0001),
        b'C' => Ok(0b0010),
        b'G' => Ok(0b0100),
        b'T' => Ok(0b1000),
        b'R' => Ok(0b0101),
        b'Y' => Ok(0b1010),
        b'S' => Ok(0b0110),
        b'W' => Ok(0b1001),
        b'K' => Ok(0b1100),
        b'M' => Ok(0b0011),
        b'B' => Ok(0b1110),
        b'D' => Ok(0b1101),
        b'H' => Ok(0b1011),
        b'V' => Ok(0b0111),
        b'N' => Ok(0b1111),
        _ => bail!("Invalid nucleotide {} in primer", char::from(nucleotide)),
    }
}

/// Enumerate exactly the released mismatch levels using owned integer sets.
/// The preflight limits and all lower-mismatch priority rules are unchanged.
fn encoded_primer_levels(
    params: &PCRParams,
    direction: PrimerDirection,
    k: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<EncodedPrimerLevels> {
    let source = if direction == PrimerDirection::Reverse {
        &params.reverse_seq
    } else {
        &params.forward_seq
    };
    let plan = plan_primer_expansion(params, direction, k)?;
    if params.trim >= k {
        gene_warn!(
            params.gene_name,
            "Trim length ({}) must be less than k ({}); adjusting trim to k-1 = {}",
            params.trim,
            k,
            k - 1
        );
    }
    if source.len() > plan.primer.len() {
        gene_info!(
            params.gene_name,
            "Trimming the primer to {} so that it is within the trim length of {}.",
            plan.primer,
            plan.primer.len()
        );
    }
    let mut exact = BudgetSet::new(budget);
    exact.insert(0u64)?;
    for nucleotide in plan.primer.bytes() {
        let mask = allowed_bases(nucleotide)?;
        let mut next = BudgetSet::new(budget);
        for &prefix in exact.iter() {
            for base in 0..4u64 {
                if mask & (1 << base) != 0 {
                    next.insert((prefix << 2) | base)?;
                }
            }
        }
        exact = next;
    }
    let mut seen = BudgetSet::new(budget);
    for &variant in exact.iter() {
        seen.insert(variant)?;
    }
    let mut levels = BudgetVec::new(budget);
    levels.push(exact)?;
    for _ in 1..=plan.mismatches {
        let mut next = BudgetSet::new(budget);
        // Iterate the same cumulative set as the released one-step permutation,
        // without cloning strings or constructing intermediate character Vecs.
        for &variant in seen.iter() {
            for position in 0..plan.primer.len() {
                let shift = 2 * position;
                let prefix = variant & !(3u64 << shift);
                for base in 0..4u64 {
                    let candidate = prefix | (base << shift);
                    if !seen.contains(&candidate) {
                        next.insert(candidate)?;
                    }
                }
            }
        }
        for &variant in next.iter() {
            seen.insert(variant)?;
        }
        levels.push(next)?;
    }
    debug!(
        "  There are {} variants of the primer across {} mismatch levels",
        seen.len(),
        levels.len()
    );
    Ok(EncodedPrimerLevels {
        length: plan.primer.len(),
        levels,
    })
}

fn visit_encoded_oligos(
    variants: &BudgetSet<u64>,
    length: usize,
    counts: &FilteredKmerCounts,
    minimum: u32,
    budget: Option<&Arc<ResourceBudget>>,
    mut visit: impl FnMut(u64, u64) -> Result<()>,
) -> Result<()> {
    let k = counts.get_k();
    ensure!(
        length > 0 && length < k && k <= 31,
        "Invalid encoded primer width {length} for k={k}"
    );
    let mut reverse = BudgetSet::new(budget);
    for variant in variants.iter() {
        reverse.insert(crate::kmer::revcomp_kmer(variant, &length))?;
    }
    let shift = 2 * (k - length);
    let rc_mask = (1u64 << (2 * length)) - 1;
    // Preserve canonical-first orientation preference and the target's own
    // count cutoff; the global filtered-view threshold must not select seeds.
    counts.visit_unfiltered(&mut |key, count| {
        if count >= u64::from(minimum) {
            if variants.contains(&(key >> shift)) {
                visit(key, count)?;
            } else if reverse.contains(&(key & rc_mask)) {
                visit(crate::kmer::revcomp_kmer(&key, &k), count)?;
            }
        }
        Ok(())
    })?;
    Ok(())
}

// Total biological ranking retained from the released full-sort selection:
// abundance descending, then oriented k-mer value ascending.
fn primer_candidate_rank(left: &(u64, u64), right: &(u64, u64)) -> std::cmp::Ordering {
    right.1.cmp(&left.1).then(left.0.cmp(&right.0))
}

/// Keep a worst-first heap with at most the remaining candidate allowance.
/// Every visited canonical entry maps to one unique oriented key in this round.
fn retain_primer_candidate(
    heap: &mut BudgetVec<(u64, u64)>,
    candidate: (u64, u64),
    cap: usize,
) -> Result<()> {
    if cap == 0 {
        return Ok(());
    }
    if heap.len() < cap {
        heap.push(candidate)?;
        let mut child = heap.len() - 1;
        while child > 0 {
            let parent = (child - 1) / 2;
            if !primer_candidate_rank(&heap[child], &heap[parent]).is_gt() {
                break;
            }
            heap.swap(child, parent);
            child = parent;
        }
    } else if primer_candidate_rank(&candidate, &heap[0]).is_lt() {
        heap[0] = candidate;
        let mut parent = 0;
        loop {
            let left = 2 * parent + 1;
            if left >= heap.len() {
                break;
            }
            let right = left + 1;
            let child =
                if right < heap.len() && primer_candidate_rank(&heap[right], &heap[left]).is_gt() {
                    right
                } else {
                    left
                };
            if !primer_candidate_rank(&heap[child], &heap[parent]).is_gt() {
                break;
            }
            heap.swap(parent, child);
            parent = child;
        }
    }
    Ok(())
}

fn discover_encoded_rounds(
    levels: &EncodedPrimerLevels,
    counts: &FilteredKmerCounts,
    params: &PCRParams,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<PrimerKmers> {
    let mut result = PrimerKmers::with_budget(&counts.get_k(), budget)?;
    for (mismatch_level, variants) in levels.levels.iter().enumerate() {
        if result.len() >= params.max_primer_kmers {
            break;
        }
        if variants.is_empty() {
            continue;
        }
        let remaining = params.max_primer_kmers - result.len();
        let mut entries = BudgetVec::new(budget);
        let mut new_count = 0u64;
        visit_encoded_oligos(
            variants,
            levels.length,
            counts,
            params.min_count,
            budget,
            |key, count| {
                if result.get(&key).is_none() {
                    new_count = new_count
                        .checked_add(1)
                        .context("Primer candidate count overflow")?;
                    retain_primer_candidate(&mut entries, (key, count), remaining)?;
                }
                Ok(())
            },
        )?;
        if entries.is_empty() {
            continue;
        }
        // Keys are unique, so this total comparator needs no stable-sort heap.
        entries.sort_unstable_by(primer_candidate_rank);
        let take = entries.len();
        for &(key, count) in entries.iter().take(take) {
            result.insert(&key, &count)?;
        }
        gene_info!(
            params.gene_name,
            "Mismatch level {}: {} new primer kmers ({} total, cap {})",
            mismatch_level,
            take,
            result.len(),
            params.max_primer_kmers
        );
        if new_count > take as u64 {
            gene_info!(
                params.gene_name,
                "Mismatch level {}: dropped {} kmers at cap",
                mismatch_level,
                new_count - take as u64
            );
        }
    }
    Ok(result)
}

/// Discover the same directional candidates under the active target ceiling.
pub(super) fn get_primer_kmers(
    params: &PCRParams,
    counts: &FilteredKmerCounts,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<(PrimerKmers, PrimerKmers)> {
    let forward_levels =
        encoded_primer_levels(params, PrimerDirection::Forward, counts.get_k(), budget)?;
    let reverse_levels =
        encoded_primer_levels(params, PrimerDirection::Reverse, counts.get_k(), budget)?;
    gene_info!(
        params.gene_name,
        "Searching kmers that contain the forward primer variants"
    );
    let forward = discover_encoded_rounds(&forward_levels, counts, params, budget)?;
    drop(forward_levels);
    gene_info!(
        params.gene_name,
        "Searching kmers that contain the reverse primer variants"
    );
    let reverse = discover_encoded_rounds(&reverse_levels, counts, params, budget)?;
    Ok((forward, reverse))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn primer_budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    fn primer_params(primer: &str, mismatches: usize) -> PCRParams {
        crate::cli::parse_pcr_primers_string(&format!("name=primer,forward={primer},reverse={primer},trim={},mismatches={mismatches},min-length=0,max-length=100",primer.len())).unwrap()
    }

    #[test]
    fn encoded_mismatch_levels_match_legacy_string_expansion_for_all_iupac_pairs() {
        let budget = primer_budget(2 << 20);
        let alphabet = b"ACGTRYSWKMBDHVN";
        for &left in alphabet {
            for &right in alphabet {
                let primer = String::from_utf8(vec![left, right]).unwrap();
                for mismatches in 0..=2 {
                    let params = primer_params(&primer, mismatches);
                    let expected =
                        preprocess_primer_by_mismatch(&params, PrimerDirection::Forward, &5)
                            .unwrap();
                    let actual =
                        encoded_primer_levels(&params, PrimerDirection::Forward, 5, Some(&budget))
                            .unwrap();
                    assert_eq!(actual.levels.len(), expected.len());
                    for (actual, expected) in actual.levels.iter().zip(expected) {
                        let expected: std::collections::BTreeSet<_> = expected
                            .iter()
                            .map(|value| string_to_oligo(value).unwrap().kmer)
                            .collect();
                        assert_eq!(
                            actual
                                .iter()
                                .copied()
                                .collect::<std::collections::BTreeSet<_>>(),
                            expected,
                            "{primer}/{mismatches}"
                        );
                    }
                    drop(actual);
                    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
                }
            }
        }
        for primer in ["ACGTAC", "ARYNCA", "NNRYSW"] {
            let params = primer_params(primer, 2);
            let expected =
                preprocess_primer_by_mismatch(&params, PrimerDirection::Reverse, &9).unwrap();
            let actual =
                encoded_primer_levels(&params, PrimerDirection::Reverse, 9, Some(&budget)).unwrap();
            for (actual, expected) in actual.levels.iter().zip(expected) {
                assert_eq!(
                    actual
                        .iter()
                        .map(|bits| crate::kmer::kmer_to_seq(bits, &primer.len()))
                        .collect::<HashSet<_>>(),
                    expected
                );
            }
            drop(actual);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        }
    }

    #[test]
    fn encoded_discovery_preserves_round_priority_caps_orientation_and_target_threshold() {
        let mut counts = KmerCounts::new(&5);
        for key in 0..1024u64 {
            let sequence = crate::kmer::kmer_to_seq(&key, &5);
            for _ in 0..(key % 5 + 1) {
                counts.ingest_seq(&sequence).unwrap();
            }
        }
        let budget = primer_budget(2 << 20);
        for primer in ["AAC", "AAR", "NCA", "YGT", "NNN"] {
            for mismatches in 0..=2 {
                for cap in [1, 2, 7, 64, 1024] {
                    for minimum in [1, 8] {
                        let mut params = primer_params(primer, mismatches);
                        params.max_primer_kmers = cap;
                        params.min_count = minimum;
                        let filtered = counts.filtered_view(100);
                        let expected = get_primer_kmers_reference(&params, &filtered).unwrap();
                        let actual = get_primer_kmers(&params, &filtered, Some(&budget)).unwrap();
                        for (old, new) in [(&expected.0, &actual.0), (&expected.1, &actual.1)] {
                            assert_eq!(
                                old.iter()
                                    .map(|(&k, &c)| (k, c))
                                    .collect::<std::collections::BTreeMap<_, _>>(),
                                new.iter()
                                    .map(|(&k, &c)| (k, c))
                                    .collect::<std::collections::BTreeMap<_, _>>(),
                                "{primer}/{mismatches}/{cap}/{minimum}"
                            );
                        }
                        drop(actual);
                        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
                    }
                }
            }
        }
    }

    #[test]
    fn primer_variant_refusal_releases_temporary_buffers_without_altering_input_counts() {
        let mut counts = KmerCounts::new(&5);
        counts.ingest_seq("AACGT").unwrap();
        let mut params = primer_params("AAR", 2);
        params.min_count = 1;
        let tiny = primer_budget(1);
        let error = get_primer_kmers(&params, &counts.filtered_view(1), Some(&tiny))
            .err()
            .unwrap();
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
        let enough = primer_budget(1 << 20);
        let actual = get_primer_kmers(&params, &counts.filtered_view(1), Some(&enough)).unwrap();
        assert!(!actual.0.is_empty());
        drop(actual);
        assert_eq!(enough.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn bounded_candidate_heap_matches_full_sort_independent_of_visit_order() {
        let values: Vec<_> = (0..501u64).map(|key| (key, key % 11)).collect();
        let mut expected = values.clone();
        expected.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        let budget = primer_budget(16384);
        let mut random = 157u64;
        for cap in [0, 1, 2, 3, 7, 16, 100] {
            for _ in 0..5 {
                let mut shuffled = values.clone();
                for index in (1..shuffled.len()).rev() {
                    random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                    shuffled.swap(index, (random >> 32) as usize % (index + 1));
                }
                let mut heap = BudgetVec::new(Some(&budget));
                for entry in shuffled {
                    retain_primer_candidate(&mut heap, entry, cap).unwrap();
                    assert!(heap.len() <= cap);
                }
                heap.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
                assert_eq!(&*heap, &expected[..cap]);
                assert!(budget.snapshot().reserved_memory_bytes <= (cap.max(4) * 2 * 16) as u64);
                drop(heap);
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            }
        }
    }

    #[test]
    fn broad_primer_scan_retains_only_the_cap_under_a_small_buffer_budget() {
        let mut counts = KmerCounts::new(&5);
        for key in 0..1024u64 {
            counts
                .ingest_seq(&crate::kmer::kmer_to_seq(&key, &5))
                .unwrap();
        }
        let mut params = primer_params("N", 0);
        params.max_primer_kmers = 2;
        params.min_count = 1;
        let expected = get_primer_kmers_reference(&params, &counts.filtered_view(1)).unwrap();
        let budget = primer_budget(4096);
        let actual = get_primer_kmers(&params, &counts.filtered_view(1), Some(&budget)).unwrap();
        assert_eq!(actual.0.len(), 2);
        assert_eq!(actual.1.len(), 2);
        assert_eq!(
            actual
                .0
                .iter()
                .map(|(&k, &c)| (k, c))
                .collect::<std::collections::BTreeMap<_, _>>(),
            expected
                .0
                .iter()
                .map(|(&k, &c)| (k, c))
                .collect::<std::collections::BTreeMap<_, _>>()
        );
        assert!(budget.snapshot().peak_reserved_memory_bytes <= 4096);
        drop(actual);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

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
        let result = find_oligos_in_kmers(&[oligo], &filtered, &1).unwrap();

        // Kmers in ACGTACGT: ACGTA, CGTAC, GTACG, TACGT
        // Forward match (starts with ACG): ACGTA
        // RC match: reverse complement of a kmer ending with CGT (=ACG rc)
        assert!(
            !result.is_empty(),
            "Should find at least one kmer matching oligo ACG"
        );
    }

    #[test]
    fn primer_scan_preserves_target_threshold_when_graph_filter_is_stricter() {
        use crate::kmer::identity::{CanonicalKey, KmerSize};
        use crate::kmer::store::MemoryCountBuilder;
        use std::num::NonZeroU64;
        let counts = build_kmer_counts("AACCCAACC", 5);
        let size = KmerSize::new(5).unwrap();
        let mut builder = MemoryCountBuilder::new(size, 4096, 8).unwrap();
        for (&bits, &count) in counts.iter() {
            builder
                .add(
                    CanonicalKey::new(bits, size).unwrap(),
                    NonZeroU64::new(u64::from(count)).unwrap(),
                )
                .unwrap();
        }
        let store = builder.finalize([0; 32]);
        for view in [
            counts.filtered_view(2),
            FilteredKmerCounts::from_store(&store, 2),
        ] {
            let oligo = string_to_oligo("AAC").unwrap();
            let matches = find_oligos_in_kmers(&[oligo], &view, &1).unwrap();
            assert!(!matches.is_empty());
            for (&bits, &count) in matches.iter() {
                assert_eq!(count, 1);
                assert_eq!(view.get_canonical(&bits).unwrap(), None);
            }
        }
    }

    #[test]
    fn test_find_oligos_no_match() {
        // Sequence of all A's — oligo "GGG" should not match
        let k = 5;
        let kmer_counts = build_kmer_counts("AAAAAAAAAA", k);
        let filtered = kmer_counts.filtered_view(1);

        let oligo = string_to_oligo("GGG").unwrap();
        let result = find_oligos_in_kmers(&[oligo], &filtered, &1).unwrap();

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
        let result = find_oligos_in_kmers(&[oligo], &filtered, &2).unwrap();

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
        let result = find_oligos_in_kmers(&[oligo], &filtered, &1).unwrap();

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
    fn test_find_oligos_oligo_equals_k_minus_1() {
        // Maximum valid oligo length is k-1: the seed node is the PREFIX of
        // the primer k-mer (kmer >> 2), which drops the last base, so trim
        // is clamped to k-1 to ensure all primer bases reach the seed node.
        let k = 5;
        let kmer_counts = build_kmer_counts("ACGTACGT", k);
        let filtered = kmer_counts.filtered_view(1);

        // oligo length = k-1 = 4: matches any kmer starting with "ACGT"
        let oligo = string_to_oligo("ACGT").unwrap();
        assert_eq!(oligo.length, k - 1);
        let result = find_oligos_in_kmers(&[oligo], &filtered, &1).unwrap();

        assert!(
            !result.is_empty(),
            "Oligo of length k-1 should match kmers starting with the primer"
        );
    }

    // --- filter_primer_kmers ---

    #[test]
    fn test_filter_primer_kmers_empty() {
        let empty = PrimerKmers::new(&5).unwrap();
        let result = filter_primer_kmers(empty, 10).unwrap();
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_filter_primer_kmers_cap() {
        // Insert 5 kmers with different counts, cap at 3
        let k = 5;
        let mut counts = PrimerKmers::new(&k).unwrap();
        for (i, seq) in ["AAAAA", "AAAAC", "AAACG", "AACGT", "ACGTA"]
            .iter()
            .enumerate()
        {
            let kmer = crate::kmer::seq_to_kmer(seq).unwrap();
            let count = (i + 1) as u64; // counts: 1, 2, 3, 4, 5
            counts.insert(&kmer, &count).unwrap();
        }
        let result = filter_primer_kmers(counts, 3).unwrap();
        // Top 3 by count: 5, 4, 3 → 3 kmers retained
        assert_eq!(result.len(), 3);
    }

    #[test]
    fn test_filter_primer_kmers_hard_cap_with_ties() {
        // Insert 5 kmers all with the SAME count, cap at 3.
        // Old implementation would keep all 5 (leaky cap).
        // New implementation must keep exactly 3, breaking ties by kmer value.
        let k = 5;
        let mut counts = PrimerKmers::new(&k).unwrap();
        for seq in ["AAAAA", "AAAAC", "AAACG", "AACGT", "ACGTA"] {
            let kmer = crate::kmer::seq_to_kmer(seq).unwrap();
            let count = 2u64; // all same count
            counts.insert(&kmer, &count).unwrap();
        }
        let result = filter_primer_kmers(counts, 3).unwrap();
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
        let mut counts = PrimerKmers::new(&k).unwrap();
        for seq in ["AAAAA", "AAAAC"] {
            let kmer = crate::kmer::seq_to_kmer(seq).unwrap();
            counts.insert(&kmer, &3).unwrap();
        }
        let result = filter_primer_kmers(counts, 5).unwrap();
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
            gene: None,
            region: None,
            index: None,
            compartment: None,
            gene_type: None,
            copy_number: None,
            deprecated: false,
            deprecated_by: None,
            deprecated_reason: None,
            forward_seq: "ACGTACGT".to_string(),
            reverse_seq: "TGCATGCA".to_string(),
            min_length: 100,
            max_length: 500,
            min_count: 2,
            mismatches: 2,
            trim: 7, // k-1: maximum valid trim for k=8
            dedup_edit_threshold: DEFAULT_DEDUP_EDIT_THRESHOLD,
            max_primer_kmers: DEFAULT_MAX_NUM_PRIMER_KMERS,
            max_dfs_states: DEFAULT_MAX_DFS_STATES,
            max_paths_per_pair: DEFAULT_MAX_PATHS_PER_PAIR,
            max_node_visits: DEFAULT_MAX_NODE_VISITS,
            high_coverage_ratio: DEFAULT_HIGH_COVERAGE_RATIO,
            tip_coverage_fraction: DEFAULT_TIP_COVERAGE_FRACTION,
            diagnose_withheld_paths: false,
            emit_uncertain: false,
            citation: String::new(),
            notes: String::new(),
            expected_length: None,
            source: String::new(),
        };
        let k = 8usize;
        let levels = preprocess_primer_by_mismatch(&params, PrimerDirection::Forward, &k).unwrap();
        assert_eq!(levels.len(), 3); // 0, 1, 2 mismatches

        // Level 0 should contain exact match (primer trimmed to last 7 bp = k-1)
        assert!(
            levels[0].contains("CGTACGT"),
            "Level 0 should contain the trimmed primer (last k-1 bases)"
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

    #[test]
    fn test_expansion_plan_matches_normal_and_iupac_variants() {
        for primer in ["ACGTAC", "ARNYAC"] {
            let params = crate::cli::parse_pcr_primers_string(&format!(
                "name=budget,forward={primer},reverse=TGCAAA,trim=6,mismatches=2"
            ))
            .unwrap();
            let plan = plan_primer_expansion(&params, PrimerDirection::Forward, 7).unwrap();
            let actual = preprocess_primer(&params, PrimerDirection::Forward, &7).unwrap();

            assert_eq!(plan.retained_variants, actual.len());
            assert!(plan.generation_work > 0);
        }
    }

    #[test]
    fn test_expansion_plan_rejects_excessive_ambiguity_before_resolution() {
        let params = crate::cli::parse_pcr_primers_string(
            "name=budget,forward=NNNNNNNN,reverse=ACGTACGT,trim=8,mismatches=0",
        )
        .unwrap();
        let error = validate_primer_expansion(&params, 9).unwrap_err();

        assert!(error.to_string().contains("too many ambiguous bases"));
    }

    #[test]
    fn test_expansion_plan_rejects_generation_work_before_permutation() {
        let primer = format!("{}{}{}", "N".repeat(5), "B".repeat(2), "A".repeat(23));
        let params = crate::cli::parse_pcr_primers_string(&format!(
            "name=budget,forward={primer},reverse=ACGTACGT,trim=30,mismatches=1"
        ))
        .unwrap();
        let error = validate_primer_expansion(&params, 31).unwrap_err();

        assert!(error.to_string().contains("generated mismatch candidates"));
    }

    #[test]
    fn test_expansion_plan_rejects_zero_effective_trim_without_panicking() {
        let params = crate::cli::parse_pcr_primers_string(
            "name=budget,forward=ACGTACGT,reverse=TGCATGCA,trim=1,mismatches=0",
        )
        .unwrap();

        for kmer_length in [0, 1] {
            let error = validate_primer_expansion(&params, kmer_length).unwrap_err();
            assert!(error.to_string().contains("k-mer length"));
        }
    }

    #[test]
    fn test_intermediate_peak_accounts_for_post_extend_growth() {
        let peak = intermediate_peak(1, 10, 100, "ACGT").unwrap();
        assert_eq!(peak, 301);

        let error = ensure_intermediate_limit(
            intermediate_peak(1, 100_000, 1_500_000, "ACGT").unwrap(),
            "ACGT",
        )
        .unwrap_err();
        assert!(error.to_string().contains("intermediate variants"));
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
