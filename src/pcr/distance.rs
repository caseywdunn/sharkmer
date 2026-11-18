//! Reserve the pinned edit-distance library's scratch before it allocates.
//!
//! editdistancek 1.0.2::edit_distance_bounded owns exactly two Vec<isize> arrays
//! of 2*k+3 elements after rejecting excessive length differences. Neither
//! array grows; its SIMD/naive prefix matching only borrows the input bytes.
//! The guard stays live through both arrays' destruction, including unwinding.
use crate::kmer::budget::{MemoryClass, ResourceBudget};
use anyhow::{Context, Result, ensure};
use std::sync::Arc;

/// Compute the exact bounded-library heap payload, checking Vec layout limits.
fn scratch_bytes(k: usize) -> Result<u64> {
    let array = k
        .checked_mul(2)
        .and_then(|n| n.checked_add(3))
        .and_then(|n| n.checked_mul(size_of::<isize>()))
        .context("Edit-distance scratch size overflow")?;
    ensure!(
        array <= isize::MAX as usize,
        "Edit-distance scratch array exceeds the addressable allocation size"
    );
    u64::try_from(array)
        .ok()
        .and_then(|n| n.checked_mul(2))
        .context("Combined edit-distance scratch size overflow")
}

/// Preserve bio 1.6's bounded-distance adapter with an owned scratch allowance.
/// Its threshold is clamped to the longer input and its result fits in u32.
pub(super) fn bounded_levenshtein(
    alpha: &[u8],
    beta: &[u8],
    threshold: u32,
    budget: Option<&Arc<ResourceBudget>>,
    class: MemoryClass,
) -> Result<Option<u32>> {
    let k = (threshold as usize).min(alpha.len().max(beta.len()));
    // This is the pinned library's own allocation-free rejection, before either
    // Vec is created. Do not require scratch for a comparison it cannot use.
    if alpha.len().abs_diff(beta.len()) > k {
        return Ok(None);
    }
    let _scratch = budget
        .map(|budget| {
            budget
                .reserve(class, scratch_bytes(k)?)
                .context("Failed to reserve edit-distance scratch")
        })
        .transpose()?;
    editdistancek::edit_distance_bounded(alpha, beta, k)
        .map(|distance| {
            u32::try_from(distance).context("Edit distance exceeds clamped u32 threshold")
        })
        .transpose()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceLimits};

    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    fn full_distance(alpha: &[u8], beta: &[u8]) -> u32 {
        let mut row: Vec<u32> = (0..=beta.len() as u32).collect();
        for (i, a) in alpha.iter().enumerate() {
            let mut diagonal = row[0];
            row[0] = i as u32 + 1;
            for (j, b) in beta.iter().enumerate() {
                let above = row[j + 1];
                row[j + 1] = (row[j] + 1)
                    .min(above + 1)
                    .min(diagonal + u32::from(a != b));
                diagonal = above;
            }
        }
        row[beta.len()]
    }

    #[test]
    fn reserved_distance_matches_full_oracle_and_previous_adapter_for_all_short_strings() {
        let mut strings = Vec::new();
        for length in 0..=5 {
            for bits in 0..1usize << length {
                strings.push(
                    (0..length)
                        .map(|bit| b"AC"[(bits >> bit) & 1])
                        .collect::<Vec<_>>(),
                );
            }
        }
        let budget = budget(4096);
        for a in &strings {
            for b in &strings {
                let exact = full_distance(a, b);
                for threshold in [0, 1, 2, 3, 5, 8, u32::MAX] {
                    let actual = bounded_levenshtein(
                        a,
                        b,
                        threshold,
                        Some(&budget),
                        MemoryClass::RetainedResults,
                    )
                    .unwrap();
                    assert_eq!(
                        actual,
                        (exact <= threshold).then_some(exact),
                        "{a:?} {b:?} k={threshold}"
                    );
                    assert_eq!(
                        actual,
                        bio::alignment::distance::simd::bounded_levenshtein(a, b, threshold)
                    );
                    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
                }
            }
        }
    }

    #[test]
    fn reserved_distance_matches_full_oracle_for_simd_prefixes_and_multiple_edits() {
        let budget = budget(1 << 20);
        let mut state = 157u64;
        for length in [16, 17, 31, 32, 33, 64, 151, 256] {
            let a: Vec<_> = (0..length)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[((state >> 32) & 3) as usize]
                })
                .collect();
            for position in [0, length / 2, length - 1] {
                let mut changed = a.clone();
                changed[position] = b'N';
                let mut inserted = a.clone();
                inserted.splice(position..position, *b"ACG");
                let mut deleted = a.clone();
                deleted.drain(position..(position + 4).min(length));
                for b in [&changed, &inserted, &deleted] {
                    let exact = full_distance(&a, b);
                    for threshold in [0, 1, 2, 4, 8, u32::MAX] {
                        for (x, y) in [(&a, b), (b, &a)] {
                            let actual = bounded_levenshtein(
                                x,
                                y,
                                threshold,
                                Some(&budget),
                                MemoryClass::RetainedResults,
                            )
                            .unwrap();
                            assert_eq!(actual, (exact <= threshold).then_some(exact));
                            assert_eq!(
                                actual,
                                bio::alignment::distance::simd::bounded_levenshtein(
                                    x, y, threshold
                                )
                            );
                            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn reserved_distance_refuses_before_library_scratch_and_retries_after_owner_drop() {
        let required = scratch_bytes(4).unwrap();
        let budget = budget(required);
        let other = budget.reserve(MemoryClass::RetainedResults, 1).unwrap();
        let err = bounded_levenshtein(
            b"AAAA",
            b"TTTT",
            4,
            Some(&budget),
            MemoryClass::RetainedResults,
        )
        .unwrap_err();
        assert!(err.downcast_ref::<ResourceExhausted>().is_some());
        assert_eq!(budget.snapshot().reserved_memory_bytes, 1);
        assert_eq!(budget.snapshot().peak_reserved_memory_bytes, 1);
        // Length-based rejection must work even when scratch cannot be admitted.
        assert_eq!(
            bounded_levenshtein(b"", b"AAAA", 3, Some(&budget), MemoryClass::RetainedResults)
                .unwrap(),
            None
        );
        drop(other);
        assert_eq!(
            bounded_levenshtein(
                b"AAAA",
                b"TTTT",
                u32::MAX,
                Some(&budget),
                MemoryClass::RetainedResults
            )
            .unwrap(),
            Some(4)
        );
        assert_eq!(
            budget.snapshot().peak_reserved_by_class[MemoryClass::RetainedResults as usize],
            required
        );
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn distance_scratch_checks_arithmetic_and_individual_array_layout_limits() {
        assert_eq!(scratch_bytes(0).unwrap(), (6 * size_of::<isize>()) as u64);
        assert_eq!(scratch_bytes(17).unwrap(), (74 * size_of::<isize>()) as u64);
        assert!(scratch_bytes(usize::MAX).is_err());
        assert!(scratch_bytes(isize::MAX as usize / size_of::<isize>()).is_err());
    }

    #[test]
    fn reserved_distance_covers_isolated_allocator_peaks_and_releases_both_arrays() {
        const ENV: &str = "SHARKMER_DISTANCE_ALLOCATION_AUDIT_CHILD";
        if std::env::var_os(ENV).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "pcr::distance::tests::reserved_distance_covers_isolated_allocator_peaks_and_releases_both_arrays", "--nocapture", "--test-threads=1"])
                .env(ENV, "1")
                .status()
                .unwrap();
            assert!(status.success());
            return;
        }
        for length in [0, 1, 19, 151, 1024, 32768] {
            let a = vec![b'A'; length];
            let mut b = a.clone();
            if let Some(last) = b.last_mut() {
                *last = b'T';
            }
            for threshold in [0, 1, 8, u32::MAX] {
                let required = scratch_bytes((threshold as usize).min(length)).unwrap();
                let budget = budget(required);
                let baseline = crate::PEAK_ALLOC.current_usage();
                crate::PEAK_ALLOC.reset_peak_usage();
                let result = bounded_levenshtein(
                    &a,
                    &b,
                    threshold,
                    Some(&budget),
                    MemoryClass::RetainedResults,
                )
                .unwrap();
                let actual = crate::PEAK_ALLOC.peak_usage() - baseline;
                let remaining = crate::PEAK_ALLOC.current_usage();
                assert_eq!(
                    result,
                    if length == 0 {
                        Some(0)
                    } else {
                        (threshold > 0).then_some(1)
                    }
                );
                assert_eq!(remaining, baseline, "library scratch leaked");
                assert_eq!(actual as u64, required, "length={length} k={threshold}");
                assert_eq!(budget.snapshot().peak_reserved_memory_bytes, required);
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
                println!(
                    "distance_allocator length={length} threshold={threshold} peak_bytes={actual} reserved_bytes={required}"
                );
            }
        }
    }
}
