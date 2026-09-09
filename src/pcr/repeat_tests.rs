use super::*;
use crate::io::{Mate, ReadRecord};
use crate::kmer::{KmerCounts, seq_to_kmer};
use petgraph::stable_graph::StableDiGraph;
use rand::{Rng, SeedableRng};

struct RepeatFixture {
    target: String,
    left_flank: String,
    right_flank: String,
    params: PCRParams,
}

fn random_dna(random: &mut rand::rngs::StdRng, length: usize) -> String {
    (0..length)
        .map(|_| b"ACGT"[random.gen_range(0..4)] as char)
        .collect()
}

fn repeat_fixture(interior: &str) -> RepeatFixture {
    let mut random = rand::rngs::StdRng::seed_from_u64(132);
    let mut left_flank = random_dna(&mut random, 60);
    let mut right_flank = random_dna(&mut random, 60);
    let first_base = interior.as_bytes()[0] as char;
    let last_base = interior.as_bytes()[interior.len() - 1] as char;
    let left_boundary = if first_base == 'A' { "C" } else { "A" };
    let right_boundary = if last_base == 'T' { "C" } else { "T" };
    left_flank.replace_range(59..60, left_boundary);
    right_flank.replace_range(0..1, right_boundary);
    let target = format!("{}{}{}", left_flank, interior, right_flank);
    let reverse =
        String::from_utf8(bio::alphabets::dna::revcomp(&right_flank.as_bytes()[45..])).unwrap();
    let params = crate::cli::parse_pcr_primers_string(&format!(
        "name=repeat,forward={},reverse={},trim=15,mismatches=0,min-length=100,max-length=220,dedup-edit-threshold=0",
        &left_flank[..15], reverse
    ))
    .unwrap();
    RepeatFixture {
        target,
        left_flank,
        right_flank,
        params,
    }
}

fn repeat_counts(targets: &[&str]) -> KmerCounts {
    let kmer_length = 19;
    let mut counts = KmerCounts::new(&kmer_length);
    for target in targets {
        for _ in 0..4 {
            counts.ingest_seq(target).unwrap();
        }
    }
    counts
}

fn run_repeat_pcr(
    counts: &KmerCounts,
    params: &PCRParams,
    reads: Option<&[ReadRecord]>,
) -> PcrOutcome {
    do_pcr(
        &counts.filtered_view(2),
        "repeat",
        params,
        false,
        "/tmp/",
        reads,
        10_000,
    )
    .unwrap()
}

fn assert_repeat_uncertainty(outcome: &PcrOutcome) {
    assert!(
        outcome.records.is_empty(),
        "unexpected products: {:?}",
        outcome
            .records
            .iter()
            .map(|record| String::from_utf8_lossy(record.seq()))
            .collect::<Vec<_>>()
    );
    assert!(
        outcome
            .failure_reason
            .as_deref()
            .is_some_and(|reason| reason.contains("repeat length unresolved")),
        "unexpected outcome: {:?}",
        outcome.failure_reason
    );
}

#[test]
fn homopolymers_at_or_above_k_are_withheld_in_every_orientation() {
    for repeat_base in ['A', 'T', 'C', 'G'] {
        let supported = repeat_fixture(&repeat_base.to_string().repeat(18));
        let supported_counts = repeat_counts(&[&supported.target]);
        let supported_outcome = run_repeat_pcr(&supported_counts, &supported.params, None);
        assert_eq!(supported_outcome.records.len(), 1);
        assert_eq!(
            supported_outcome.records[0].seq(),
            supported.target.as_bytes()
        );

        for repeat_length in [19, 40] {
            let unresolved = repeat_fixture(&repeat_base.to_string().repeat(repeat_length));
            let unresolved_counts = repeat_counts(&[&unresolved.target]);
            let unresolved_outcome = run_repeat_pcr(&unresolved_counts, &unresolved.params, None);
            assert_repeat_uncertainty(&unresolved_outcome);
        }
    }
}

#[test]
fn forty_base_homopolymer_never_emits_collapsed_product() {
    let fixture = repeat_fixture(&"A".repeat(40));
    let counts = repeat_counts(&[&fixture.target]);

    let outcome = run_repeat_pcr(&counts, &fixture.params, None);

    assert_eq!(fixture.target.len(), 160);
    assert_repeat_uncertainty(&outcome);
}

#[test]
fn repeat_uncertainty_precedes_collapsed_length_rejection() {
    let mut fixture = repeat_fixture(&"A".repeat(40));
    fixture.params.min_length = 155;
    fixture.params.max_length = 165;
    let counts = repeat_counts(&[&fixture.target]);

    let outcome = run_repeat_pcr(&counts, &fixture.params, None);

    assert_repeat_uncertainty(&outcome);
}

#[test]
fn tandem_repeat_cycle_and_first_visit_shortcut_are_withheld() {
    let fixture = repeat_fixture(&"AC".repeat(20));
    let counts = repeat_counts(&[&fixture.target]);
    let outcome = run_repeat_pcr(&counts, &fixture.params, None);
    assert_repeat_uncertainty(&outcome);

    let mut graph = StableDiGraph::new();
    let start = graph.add_node(DBNode {
        sub_kmer: seq_to_kmer("AAA").unwrap(),
        is_start: true,
        is_end: false,
    });
    let cycle_first = graph.add_node(DBNode {
        sub_kmer: seq_to_kmer("AAC").unwrap(),
        is_start: false,
        is_end: false,
    });
    let cycle_second = graph.add_node(DBNode {
        sub_kmer: seq_to_kmer("ACA").unwrap(),
        is_start: false,
        is_end: false,
    });
    let end = graph.add_node(DBNode {
        sub_kmer: seq_to_kmer("CAA").unwrap(),
        is_start: false,
        is_end: true,
    });
    let start_edge = graph.add_edge(
        start,
        cycle_first,
        DBEdge {
            count: 4,
            coverage_ratio: 1.0,
        },
    );
    graph.add_edge(
        cycle_first,
        cycle_second,
        DBEdge {
            count: 4,
            coverage_ratio: 1.0,
        },
    );
    let cycle_edge = graph.add_edge(
        cycle_second,
        cycle_first,
        DBEdge {
            count: 4,
            coverage_ratio: 1.0,
        },
    );
    let end_edge = graph.add_edge(
        cycle_first,
        end,
        DBEdge {
            count: 4,
            coverage_ratio: 1.0,
        },
    );
    let unresolved = graph::unresolved_repeat_sub_kmers(&graph, ahash::AHashSet::new());
    graph.remove_edge(cycle_edge);
    let shortcut = vec![
        (start, None),
        (cycle_first, Some(start_edge)),
        (end, Some(end_edge)),
    ];
    let (resolved, unresolved_count) =
        paths::filter_unresolved_repeat_paths(&graph, vec![shortcut], &unresolved);

    assert!(resolved.is_empty());
    assert_eq!(unresolved_count, 1);
}

#[test]
fn repeat_branch_does_not_suppress_disjoint_clean_path() {
    let fixture = repeat_fixture(&"AC".repeat(20));
    let mut random = rand::rngs::StdRng::seed_from_u64(13_200);
    let clean_interior = random_dna(&mut random, 40);
    let clean_target = format!(
        "{}{}{}",
        fixture.left_flank, clean_interior, fixture.right_flank
    );
    let counts = repeat_counts(&[&fixture.target, &clean_target]);

    let spanning_reads = vec![ReadRecord {
        sequence: clean_target.clone(),
        index: 0,
        mate: Mate::Unpaired,
    }];
    let nonspanning_reads = vec![
        ReadRecord {
            sequence: fixture.left_flank[..40].to_string(),
            index: 0,
            mate: Mate::Unpaired,
        },
        ReadRecord {
            sequence: fixture.right_flank[20..].to_string(),
            index: 1,
            mate: Mate::Unpaired,
        },
    ];

    for reads in [
        None,
        Some(spanning_reads.as_slice()),
        Some(nonspanning_reads.as_slice()),
    ] {
        let outcome = run_repeat_pcr(&counts, &fixture.params, reads);
        assert_eq!(outcome.records.len(), 1);
        assert_eq!(outcome.records[0].seq(), clean_target.as_bytes());
    }
}

#[test]
fn repeat_evidence_without_connectivity_is_explicit() {
    let fixture = repeat_fixture(&"A".repeat(40));
    let mut random = rand::rngs::StdRng::seed_from_u64(132_001);
    let forward_arm = format!(
        "{}{}{}",
        fixture.left_flank,
        "A".repeat(40),
        random_dna(&mut random, 40)
    );
    let reverse_arm = format!("{}{}", random_dna(&mut random, 40), fixture.right_flank);
    let counts = repeat_counts(&[&forward_arm, &reverse_arm]);

    let outcome = run_repeat_pcr(&counts, &fixture.params, None);

    assert!(outcome.records.is_empty());
    assert_eq!(
        outcome.failure_reason.as_deref(),
        Some("repeat evidence encountered before start-to-end connectivity could be established")
    );
}

#[test]
fn read_threading_does_not_guess_repeat_copy_count() {
    let fixture = repeat_fixture(&"A".repeat(40));
    let counts = repeat_counts(&[&fixture.target]);
    let spanning_reads = vec![ReadRecord {
        sequence: fixture.target.clone(),
        index: 0,
        mate: Mate::Unpaired,
    }];
    let nonspanning_reads = vec![
        ReadRecord {
            sequence: fixture.left_flank[..40].to_string(),
            index: 0,
            mate: Mate::Unpaired,
        },
        ReadRecord {
            sequence: fixture.right_flank[20..].to_string(),
            index: 1,
            mate: Mate::Unpaired,
        },
    ];

    let spanning_outcome = run_repeat_pcr(&counts, &fixture.params, Some(&spanning_reads));
    let nonspanning_outcome = run_repeat_pcr(&counts, &fixture.params, Some(&nonspanning_reads));

    assert_repeat_uncertainty(&spanning_outcome);
    assert_repeat_uncertainty(&nonspanning_outcome);
}
