use super::*;
use crate::kmer::KmerCounts;
use petgraph::stable_graph::StableDiGraph;
use rand::{Rng, SeedableRng};

fn synthetic_targets(lengths: &[usize]) -> (Vec<String>, PCRParams) {
    let mut random = rand::rngs::StdRng::seed_from_u64(131);
    let mut random_dna = |length| -> String {
        (0..length)
            .map(|_| b"ACGT"[random.gen_range(0..4)] as char)
            .collect()
    };
    let prefix = random_dna(25);
    let suffix = random_dna(25);
    let targets = lengths
        .iter()
        .map(|length| {
            format!(
                "{}{}{}",
                prefix,
                random_dna(length - prefix.len() - suffix.len()),
                suffix
            )
        })
        .collect();
    let reverse =
        String::from_utf8(bio::alphabets::dna::revcomp(&suffix.as_bytes()[10..])).unwrap();
    let params = crate::cli::parse_pcr_primers_string(&format!(
        "name=threshold,forward={},reverse={},trim=15,mismatches=0,min-length=170,max-length=210,dedup-edit-threshold=0",
        &prefix[..15], reverse
    ))
    .unwrap();
    (targets, params)
}

fn counts_for(targets: &[(&str, usize)]) -> KmerCounts {
    let kmer_length = 19;
    let mut counts = KmerCounts::new(&kmer_length);
    for (target, copies) in targets {
        for _ in 0..*copies {
            counts.ingest_seq(target).unwrap();
        }
    }
    counts
}

fn run_pcr(counts: &KmerCounts, params: &PCRParams) -> PcrOutcome {
    do_pcr(
        &counts.filtered_view(2),
        "threshold",
        params,
        false,
        "/tmp/",
        None,
        10_000,
    )
    .unwrap()
}

#[test]
fn lower_threshold_is_tried_after_connected_path_is_too_long() {
    let (targets, params) = synthetic_targets(&[400, 180]);
    let counts = counts_for(&[(&targets[0], 100), (&targets[1], 4)]);

    let outcome = run_pcr(&counts, &params);

    assert_eq!(
        outcome.records.len(),
        1,
        "unexpected outcome: {:?}",
        outcome.failure_reason
    );
    assert_eq!(outcome.records[0].seq(), targets[1].as_bytes());
}

#[test]
fn out_of_range_products_report_length_failure() {
    for target_length in [160, 400] {
        let (targets, params) = synthetic_targets(&[target_length]);
        let counts = counts_for(&[(&targets[0], 100)]);

        let outcome = run_pcr(&counts, &params);

        assert!(outcome.records.is_empty());
        assert_eq!(
            outcome.failure_reason.as_deref(),
            Some("connectivity found but no path satisfied the requested length range")
        );
    }
}

#[test]
fn exact_length_boundaries_are_valid() {
    for target_length in [170, 210] {
        let (targets, mut params) = synthetic_targets(&[target_length]);
        params.min_length = target_length;
        params.max_length = target_length;
        let counts = counts_for(&[(&targets[0], 10)]);

        let outcome = run_pcr(&counts, &params);

        assert_eq!(outcome.records.len(), 1);
        assert_eq!(outcome.records[0].seq(), targets[0].as_bytes());
    }
}

#[test]
fn standard_mode_stops_at_first_threshold_with_a_valid_product() {
    let (targets, params) = synthetic_targets(&[200, 180]);
    let counts = counts_for(&[(&targets[0], 100), (&targets[1], 4)]);

    let outcome = run_pcr(&counts, &params);

    assert_eq!(outcome.records.len(), 1);
    assert_eq!(outcome.records[0].seq(), targets[0].as_bytes());
}

#[test]
fn max_length_shorter_than_k_is_explicitly_unsupported() {
    let (targets, mut params) = synthetic_targets(&[180]);
    params.max_length = 18;
    let counts = counts_for(&[(&targets[0], 10)]);

    let outcome = run_pcr(&counts, &params);

    assert!(outcome.records.is_empty());
    assert_eq!(
        outcome.failure_reason.as_deref(),
        Some("max-length 18 is shorter than k-mer length 19")
    );
}

#[test]
fn disconnected_graph_after_pruning_has_explicit_failure() {
    let mut graph = StableDiGraph::new();
    graph.add_node(DBNode {
        sub_kmer: 0,
        is_start: true,
        is_end: false,
    });
    graph.add_node(DBNode {
        sub_kmer: 1,
        is_start: false,
        is_end: true,
    });
    let mut counts = KmerCounts::new(&3);
    counts.insert(&0, &10);
    let filtered = counts.filtered_view(1);
    let (_, mut params) = synthetic_targets(&[180]);
    params.min_length = 0;
    params.max_length = 100;

    let evaluation = evaluate_threshold_graph(
        graph,
        1,
        &filtered,
        "threshold",
        &params,
        false,
        "/tmp/",
        None,
    )
    .unwrap();

    assert!(evaluation.records.is_empty());
    assert_eq!(
        evaluation.failure_reason.as_deref(),
        Some("connectivity found but no start-to-end path remained after pruning")
    );
}

#[test]
fn search_limits_have_explicit_failures() {
    let (targets, params) = synthetic_targets(&[180]);
    let counts = counts_for(&[(&targets[0], 10)]);

    let mut dfs_limited = params.clone();
    dfs_limited.max_dfs_states = 0;
    let dfs_outcome = run_pcr(&counts, &dfs_limited);
    assert!(dfs_outcome.records.is_empty());
    assert_eq!(
        dfs_outcome.failure_reason.as_deref(),
        Some("DFS state limit reached before a valid amplicon was found")
    );

    let mut path_limited = params;
    path_limited.max_paths_per_pair = 0;
    let path_outcome = run_pcr(&counts, &path_limited);
    assert!(path_outcome.records.is_empty());
    assert_eq!(
        path_outcome.failure_reason.as_deref(),
        Some("path limit reached before a valid amplicon was found")
    );
}

#[test]
fn node_budget_failure_remains_explicit() {
    let (targets, params) = synthetic_targets(&[180]);
    let counts = counts_for(&[(&targets[0], 10)]);

    let outcome = do_pcr(
        &counts.filtered_view(2),
        "threshold",
        &params,
        false,
        "/tmp/",
        None,
        1,
    )
    .unwrap();

    assert!(outcome.records.is_empty());
    assert_eq!(
        outcome.failure_reason.as_deref(),
        Some("node budget exceeded")
    );
}
