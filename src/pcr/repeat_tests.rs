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

fn path_search_counts() -> KmerCounts {
    let kmer_length = 19;
    let mut counts = KmerCounts::new(&kmer_length);
    counts.insert(&0, &10);
    counts
}

fn path_search_params() -> PCRParams {
    let mut params = repeat_fixture("AC").params;
    params.min_length = 0;
    params.max_length = 100;
    params.max_dfs_states = 10_000;
    params.max_paths_per_pair = 20;
    params.max_node_visits = 2;
    params
}

fn repeat_quota_graph(
    repeat_depth: usize,
    clean_route_count: usize,
) -> (
    StableDiGraph<DBNode, DBEdge>,
    ahash::AHashSet<u64>,
    Vec<petgraph::graph::NodeIndex>,
) {
    let mut graph = StableDiGraph::new();
    let mut next_sub_kmer = 1_u64;
    let start = graph.add_node(DBNode {
        sub_kmer: next_sub_kmer,
        is_start: true,
        is_end: false,
    });
    next_sub_kmer += 1;
    let repeat_root = graph.add_node(DBNode {
        sub_kmer: next_sub_kmer,
        is_start: false,
        is_end: false,
    });
    next_sub_kmer += 1;
    graph.add_edge(
        start,
        repeat_root,
        DBEdge {
            count: 100,
            coverage_ratio: 1.0,
        },
    );

    let mut frontier = vec![repeat_root];
    for level in 0..repeat_depth {
        let mut next_frontier = Vec::new();
        for parent in frontier {
            for _branch in 0..2 {
                let child = graph.add_node(DBNode {
                    sub_kmer: next_sub_kmer,
                    is_start: false,
                    is_end: level + 1 == repeat_depth,
                });
                next_sub_kmer += 1;
                graph.add_edge(
                    parent,
                    child,
                    DBEdge {
                        count: 100,
                        coverage_ratio: 1.0,
                    },
                );
                next_frontier.push(child);
            }
        }
        frontier = next_frontier;
    }

    let mut clean_ends = Vec::new();
    for route_index in 0..clean_route_count {
        let clean_middle = graph.add_node(DBNode {
            sub_kmer: next_sub_kmer,
            is_start: false,
            is_end: false,
        });
        next_sub_kmer += 1;
        let clean_end = graph.add_node(DBNode {
            sub_kmer: next_sub_kmer,
            is_start: false,
            is_end: true,
        });
        next_sub_kmer += 1;
        let edge_count = u32::try_from(clean_route_count - route_index).unwrap();
        graph.add_edge(
            start,
            clean_middle,
            DBEdge {
                count: edge_count,
                coverage_ratio: 1.0,
            },
        );
        graph.add_edge(
            clean_middle,
            clean_end,
            DBEdge {
                count: edge_count,
                coverage_ratio: 1.0,
            },
        );
        clean_ends.push(clean_end);
    }

    let unresolved = ahash::AHashSet::from_iter([graph[repeat_root].sub_kmer]);
    (graph, unresolved, clean_ends)
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
    graph.add_edge(
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
    graph.add_edge(
        cycle_first,
        end,
        DBEdge {
            count: 4,
            coverage_ratio: 1.0,
        },
    );
    let unresolved = graph::unresolved_repeat_sub_kmers(&graph, ahash::AHashSet::new());
    graph.remove_edge(cycle_edge);
    let counts = path_search_counts();
    let result = paths::search_assembly_paths(
        &graph,
        &counts.filtered_view(1),
        &path_search_params(),
        None,
        &unresolved,
    );

    assert!(result.paths.is_empty());
    assert_eq!(result.unresolved_repeat_path_count, 1);
}

#[test]
fn repeat_rejected_paths_do_not_consume_clean_path_quota() {
    let (graph, unresolved, clean_ends) = repeat_quota_graph(5, 1);
    let counts = path_search_counts();
    let result = paths::search_assembly_paths(
        &graph,
        &counts.filtered_view(1),
        &path_search_params(),
        None,
        &unresolved,
    );

    assert_eq!(result.unresolved_repeat_path_count, 32);
    assert_eq!(result.paths.len(), 1);
    assert_eq!(result.paths[0].last().unwrap().0, clean_ends[0]);
    assert!(!result.path_limit_reached);
    assert!(!result.dfs_limit_reached);
}

#[test]
fn all_repetitive_search_remains_bounded_by_dfs_budget() {
    let (graph, unresolved, _) = repeat_quota_graph(5, 0);
    let counts = path_search_counts();
    let mut params = path_search_params();
    params.max_dfs_states = 10;
    let result =
        paths::search_assembly_paths(&graph, &counts.filtered_view(1), &params, None, &unresolved);

    assert!(result.paths.is_empty());
    assert!(result.unresolved_repeat_path_count > 0);
    assert!(result.dfs_limit_reached);
    assert!(!result.path_limit_reached);
}

#[test]
fn clean_candidates_still_consume_path_quota() {
    let (graph, unresolved, clean_ends) = repeat_quota_graph(1, 3);
    let counts = path_search_counts();
    let mut params = path_search_params();
    params.max_paths_per_pair = 2;
    let result =
        paths::search_assembly_paths(&graph, &counts.filtered_view(1), &params, None, &unresolved);

    assert_eq!(result.paths.len(), 2);
    assert_eq!(result.unresolved_repeat_path_count, 2);
    assert!(result.path_limit_reached);
    assert!(
        result
            .paths
            .iter()
            .all(|path| clean_ends.contains(&path.last().unwrap().0))
    );
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
