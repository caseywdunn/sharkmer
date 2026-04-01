// pcr/seed_eval.rs — Bounded seed evaluation before full graph extension

use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;
use std::collections::HashMap;

use crate::kmer::FilteredKmerCounts;
use crate::kmer::encoding::kmers_from_ascii;

use super::graph::get_suffix_mask;
use super::{DBEdge, DBNode, PCRParams};

/// Maximum number of nodes to explore per seed during evaluation.
const MAX_NODES_PER_SEED: usize = 500;

/// Abandon a seed if branching ratio exceeds this (too branchy overall).
const MAX_BRANCHING_RATIO: f64 = 0.4;

/// Abandon a seed if it terminates early with fewer than this fraction of
/// expected linear nodes.
const MIN_EXTENSION_FRACTION: f64 = 0.1;

/// Abandon a seed if per-seed budget was exhausted AND branching ratio
/// exceeds this threshold.
const BUDGET_BRANCHING_THRESHOLD: f64 = 0.2;

/// Metrics collected during bounded extension of a single seed.
struct SeedMetrics {
    node_count: usize,
    branching_events: usize,
    budget_exhausted: bool,
    reached_opposite: bool,
    terminated: bool,
}

/// Direction a seed extends.
enum SeedDirection {
    Forward,
    Reverse,
}

/// Evaluate each seed node with a bounded local exploration. Seeds that
/// look off-target are marked terminal so they do not participate in
/// full graph extension. When retained reads are available, also checks
/// for read divergence (reads branching within k bases of the primer).
pub(super) fn evaluate_seeds(
    graph: &mut StableDiGraph<DBNode, DBEdge>,
    node_lookup: &HashMap<u64, NodeIndex>,
    kmer_counts: &FilteredKmerCounts,
    params: &PCRParams,
    min_count: u32,
    retained_reads: &[&str],
) {
    let k = kmer_counts.get_k();
    let expected_linear_nodes = params.max_length.saturating_sub(k);

    // Collect seed nodes to evaluate (cannot mutate graph while iterating)
    let seeds: Vec<(NodeIndex, bool, bool, u64)> = graph
        .node_indices()
        .filter(|&node| {
            let data = &graph[node];
            data.is_start || data.is_end
        })
        .map(|node| {
            let data = &graph[node];
            (node, data.is_start, data.is_end, data.sub_kmer)
        })
        .collect();

    for (node_idx, is_start, is_end, sub_kmer) in seeds {
        // Skip dual-role nodes unconditionally
        if is_start && is_end {
            gene_info!(
                params.gene_name,
                "Seed {} (dual): Kept — dual-role node (both start and end)",
                node_idx.index()
            );
            continue;
        }

        let direction = if is_start {
            SeedDirection::Forward
        } else {
            SeedDirection::Reverse
        };

        let metrics = bounded_extend(
            sub_kmer,
            &direction,
            kmer_counts,
            min_count,
            k,
            node_lookup,
            graph,
        );

        let direction_str = match direction {
            SeedDirection::Forward => "forward",
            SeedDirection::Reverse => "reverse",
        };

        // Keep unconditionally if opposite seed was reached
        if metrics.reached_opposite {
            gene_info!(
                params.gene_name,
                "Seed {} ({}): Kept — reached opposite-direction seed (early product recovery)",
                node_idx.index(),
                direction_str
            );
            continue;
        }

        let branching_ratio = if metrics.node_count > 0 {
            metrics.branching_events as f64 / metrics.node_count as f64
        } else {
            0.0
        };

        // Check abandonment criteria
        if metrics.budget_exhausted && branching_ratio > BUDGET_BRANCHING_THRESHOLD {
            gene_info!(
                params.gene_name,
                "Seed {} ({}): Abandoned — budget exhausted with high branching (ratio {:.2} > {:.2}, {} nodes)",
                node_idx.index(),
                direction_str,
                branching_ratio,
                BUDGET_BRANCHING_THRESHOLD,
                metrics.node_count
            );
            graph[node_idx].is_terminal = true;
            graph[node_idx].visited = true;
            continue;
        }

        if branching_ratio > MAX_BRANCHING_RATIO {
            gene_info!(
                params.gene_name,
                "Seed {} ({}): Abandoned — excessive branching (ratio {:.2} > {:.2}, {} nodes)",
                node_idx.index(),
                direction_str,
                branching_ratio,
                MAX_BRANCHING_RATIO,
                metrics.node_count
            );
            graph[node_idx].is_terminal = true;
            graph[node_idx].visited = true;
            continue;
        }

        let min_nodes = (expected_linear_nodes as f64 * MIN_EXTENSION_FRACTION) as usize;
        if metrics.node_count < min_nodes && metrics.terminated {
            gene_info!(
                params.gene_name,
                "Seed {} ({}): Abandoned — too small and terminated ({} nodes < {} minimum)",
                node_idx.index(),
                direction_str,
                metrics.node_count,
                min_nodes
            );
            graph[node_idx].is_terminal = true;
            graph[node_idx].visited = true;
            continue;
        }

        // Read divergence check: if retained reads are available,
        // thread them through the bounded subgraph and check for
        // divergence within k bases of the primer.
        if !retained_reads.is_empty() {
            let divergent = check_read_divergence(
                sub_kmer,
                &direction,
                kmer_counts,
                min_count,
                k,
                retained_reads,
            );
            if divergent {
                gene_info!(
                    params.gene_name,
                    "Seed {} ({}): Abandoned — read divergence detected (reads branch within k bases of primer)",
                    node_idx.index(),
                    direction_str
                );
                graph[node_idx].is_terminal = true;
                graph[node_idx].visited = true;
                continue;
            }
        }

        gene_info!(
            params.gene_name,
            "Seed {} ({}): Kept — passed evaluation ({} nodes, branching ratio {:.2})",
            node_idx.index(),
            direction_str,
            metrics.node_count,
            branching_ratio
        );
    }
}

/// Check if retained reads show divergence at a seed.
///
/// For each read, extract kmers and trace them through the local kmer
/// graph starting from the seed. Divergence means reads take different
/// paths within the first k bases after the primer — evidence of
/// off-target binding in repetitive regions.
///
/// Returns true if divergence is detected, false otherwise (including
/// when there are 0 or 1 matching reads — not enough to detect divergence).
fn check_read_divergence(
    seed_sub_kmer: u64,
    direction: &SeedDirection,
    _kmer_counts: &FilteredKmerCounts,
    _min_count: u32,
    k: usize,
    retained_reads: &[&str],
) -> bool {
    // We need at least 2 reads that extend from this seed to detect divergence
    let suffix_mask = get_suffix_mask(&k);

    // For each read, find the sequence of sub_kmers it follows from the seed.
    // We only care about the first k nodes (divergence within k bases).
    let max_trace_depth = k;
    let mut traces: Vec<Vec<u64>> = Vec::new();

    for read_seq in retained_reads {
        let read_kmers = match kmers_from_ascii(read_seq, k) {
            Ok(k) => k,
            Err(_) => continue,
        };

        // Check if any kmer in this read corresponds to the seed sub_kmer
        let mut trace = Vec::new();
        let mut found_seed = false;

        for &kmer in &read_kmers {
            // Check both orientations of this kmer against the seed
            let sub_kmer_prefix = kmer >> 2; // prefix (k-1 mer)
            let sub_kmer_suffix = kmer & suffix_mask; // suffix (k-1 mer)
            let rc_kmer = crate::kmer::revcomp_kmer(&kmer, &k);
            let rc_prefix = rc_kmer >> 2;
            let rc_suffix = rc_kmer & suffix_mask;

            if !found_seed {
                // Look for the seed sub_kmer as a prefix or suffix
                if sub_kmer_prefix == seed_sub_kmer
                    || sub_kmer_suffix == seed_sub_kmer
                    || rc_prefix == seed_sub_kmer
                    || rc_suffix == seed_sub_kmer
                {
                    found_seed = true;
                    // Start tracing from the next position
                    match direction {
                        SeedDirection::Forward => {
                            trace.push(sub_kmer_suffix);
                        }
                        SeedDirection::Reverse => {
                            trace.push(sub_kmer_prefix);
                        }
                    }
                }
            } else if trace.len() < max_trace_depth {
                // Continue tracing: record the next sub_kmer in the read's direction
                match direction {
                    SeedDirection::Forward => {
                        trace.push(sub_kmer_suffix);
                    }
                    SeedDirection::Reverse => {
                        trace.push(sub_kmer_prefix);
                    }
                }
            } else {
                break;
            }
        }

        if found_seed && !trace.is_empty() {
            traces.push(trace);
        }
    }

    // Need at least 2 traces to detect divergence
    if traces.len() < 2 {
        return false;
    }

    // Check for divergence: do traces share the same path for the first
    // few nodes? If any two traces diverge at the same position within
    // the first k nodes, that's divergence.
    let reference_trace = &traces[0];
    let mut n_divergent = 0;
    for trace in &traces[1..] {
        let min_len = reference_trace.len().min(trace.len());
        if min_len == 0 {
            continue;
        }
        // Check if they diverge at position 0 (immediately after seed)
        if trace[0] != reference_trace[0] {
            n_divergent += 1;
        }
    }

    // Divergence: majority of reads diverge immediately from the reference
    let total_comparisons = traces.len() - 1;
    if total_comparisons > 0 && n_divergent > total_comparisons / 2 {
        return true;
    }

    false
}

/// Perform a bounded extension from a single seed node, building a
/// throwaway local exploration for evaluation only. Checks the main
/// graph's node_lookup to detect convergence with opposite-direction seeds.
fn bounded_extend(
    seed_sub_kmer: u64,
    direction: &SeedDirection,
    kmer_counts: &FilteredKmerCounts,
    min_count: u32,
    k: usize,
    main_node_lookup: &HashMap<u64, NodeIndex>,
    main_graph: &StableDiGraph<DBNode, DBEdge>,
) -> SeedMetrics {
    let suffix_mask = get_suffix_mask(&k);
    let prefix_shift = 2 * (k - 1);

    // Track explored sub_kmers to avoid revisiting
    let mut explored: HashMap<u64, ()> = HashMap::new();
    let mut frontier: Vec<u64> = vec![seed_sub_kmer];
    explored.insert(seed_sub_kmer, ());

    let mut node_count: usize = 0;
    let mut branching_events: usize = 0;
    let mut reached_opposite = false;
    let mut budget_exhausted = false;
    let mut terminated = true; // assume terminated unless budget runs out

    while let Some(current_sub_kmer) = frontier.pop() {
        if node_count >= MAX_NODES_PER_SEED {
            budget_exhausted = true;
            terminated = false;
            break;
        }

        node_count += 1;

        let mut successors: Vec<u64> = Vec::new();

        match direction {
            SeedDirection::Forward => {
                // Forward: extend rightward using suffix-based lookup
                for base in 0u64..4 {
                    let candidate_kmer = (current_sub_kmer << 2) | base;
                    if let Some(count) = kmer_counts.get_canonical(&candidate_kmer) {
                        if count >= min_count {
                            let new_sub_kmer = candidate_kmer & suffix_mask;
                            successors.push(new_sub_kmer);

                            // Check if we reached an end seed in the main graph
                            if let Some(&main_node_idx) = main_node_lookup.get(&new_sub_kmer) {
                                if main_graph[main_node_idx].is_end {
                                    reached_opposite = true;
                                }
                            }
                        }
                    }
                }
            }
            SeedDirection::Reverse => {
                // Reverse: extend leftward using prefix-based lookup
                for base in 0u64..4 {
                    let candidate_kmer = (base << prefix_shift) | current_sub_kmer;
                    if let Some(count) = kmer_counts.get_canonical(&candidate_kmer) {
                        if count >= min_count {
                            let new_sub_kmer = candidate_kmer >> 2;
                            successors.push(new_sub_kmer);

                            // Check if we reached a start seed in the main graph
                            if let Some(&main_node_idx) = main_node_lookup.get(&new_sub_kmer) {
                                if main_graph[main_node_idx].is_start {
                                    reached_opposite = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Count branching events
        if successors.len() > 1 {
            branching_events += 1;
        }

        // Add new sub_kmers to frontier
        for new_sub_kmer in successors {
            if explored.contains_key(&new_sub_kmer) {
                continue;
            }
            explored.insert(new_sub_kmer, ());
            frontier.push(new_sub_kmer);
        }
    }

    SeedMetrics {
        node_count,
        branching_events,
        budget_exhausted,
        reached_opposite,
        terminated,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::KmerCounts;
    use crate::pcr::DEFAULT_DEDUP_EDIT_THRESHOLD;

    fn make_params() -> PCRParams {
        PCRParams {
            forward_seq: "AAAA".to_string(),
            reverse_seq: "TTTT".to_string(),
            min_length: 0,
            max_length: 2500,
            gene_name: "test".to_string(),
            min_count: 1,
            mismatches: 0,
            trim: 0,
            citation: "".to_string(),
            notes: "".to_string(),
            dedup_edit_threshold: DEFAULT_DEDUP_EDIT_THRESHOLD,
            source: "test".to_string(),
        }
    }

    #[test]
    fn test_linear_seed_kept() {
        // Create a linear chain: 24 bases, k=7 => 18 kmers forming a linear path.
        // Use a short max_length so that the chain is long enough relative to
        // expected_linear_nodes (max_length - k) * MIN_EXTENSION_FRACTION.
        let k = 7;
        let seq = "ACGTACGTACGTACGTACGTACGT"; // 24 bases => 18 kmers
        let mut kc = KmerCounts::new(&k);
        for _ in 0..10 {
            kc.ingest_seq(seq).unwrap();
        }
        let fkc = kc.filtered_view(1);

        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let mut node_lookup: HashMap<u64, NodeIndex> = HashMap::new();

        // First kmer's prefix (first k-1 bases) is the start node's sub_kmer
        let first_kmer = crate::kmer::seq_to_kmer(&seq[..k].to_string()).unwrap();
        let prefix = first_kmer >> 2;

        let node = graph.add_node(DBNode {
            sub_kmer: prefix,
            is_start: true,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        node_lookup.insert(prefix, node);

        // max_length=30 => expected_linear_nodes=23, min_nodes=2
        // The chain has ~18 nodes which is above 2, so it passes.
        let mut params = make_params();
        params.max_length = 30;

        evaluate_seeds(&mut graph, &node_lookup, &fkc, &params, 1, &[]);

        // The seed should NOT be marked terminal (linear chain, not branchy)
        assert!(
            !graph[node].is_terminal,
            "Linear seed should not be abandoned"
        );
    }

    #[test]
    fn test_branchy_seed_abandoned() {
        // Create kmer counts where from a seed sub_kmer, every step fans out
        // to all 4 bases, creating massive branching. Use k=15 to avoid
        // sub_kmer space collisions. With 5 levels of full 4-way branching
        // the reachable tree exceeds MAX_NODES_PER_SEED (500), so the budget
        // is exhausted with branching ratio above BUDGET_BRANCHING_THRESHOLD.
        let k: usize = 15;
        let suffix_mask: u64 = (1 << (2 * (k - 1))) - 1;

        let mut kc = KmerCounts::new(&k);

        // Use a non-trivial start sub_kmer to avoid self-loops
        // 14 bases: ACGTACGTACGTAC = 0b_00_01_10_11_00_01_10_11_00_01_10_11_00_01
        let start_sub_kmer: u64 = 0b_00_01_10_11_00_01_10_11_00_01_10_11_00_01;

        // For 5 levels of branching, add all possible extension kmers so
        // every node fans out to 4 successors.
        let mut frontier: Vec<u64> = vec![start_sub_kmer];
        for _level in 0..5 {
            let mut next_frontier: Vec<u64> = Vec::new();
            for &sub_kmer in &frontier {
                for base in 0u64..4 {
                    let candidate_kmer = (sub_kmer << 2) | base;
                    kc.insert(&candidate_kmer, &10);
                    let new_sub_kmer = candidate_kmer & suffix_mask;
                    next_frontier.push(new_sub_kmer);
                }
            }
            frontier = next_frontier;
        }

        let fkc = kc.filtered_view(1);

        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let mut node_lookup: HashMap<u64, NodeIndex> = HashMap::new();

        let node = graph.add_node(DBNode {
            sub_kmer: start_sub_kmer,
            is_start: true,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        node_lookup.insert(start_sub_kmer, node);

        let params = make_params();

        evaluate_seeds(&mut graph, &node_lookup, &fkc, &params, 1, &[]);

        // The seed should be marked terminal (budget exhausted with high branching)
        assert!(graph[node].is_terminal, "Branchy seed should be abandoned");
    }

    #[test]
    fn test_small_dead_end_seed_abandoned() {
        // Create a seed with only 1-2 reachable kmers, then it terminates.
        let k: usize = 7;

        let mut kc = KmerCounts::new(&k);

        // Start sub_kmer for the seed
        let start_sub_kmer: u64 = 0b_00_01_10_11_00_01; // ACGTAC in 2-bit

        // Add just one extension kmer
        let candidate_kmer = (start_sub_kmer << 2) | 0b11; // extends with T
        kc.insert(&candidate_kmer, &10);

        let fkc = kc.filtered_view(1);

        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let mut node_lookup: HashMap<u64, NodeIndex> = HashMap::new();

        let node = graph.add_node(DBNode {
            sub_kmer: start_sub_kmer,
            is_start: true,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        node_lookup.insert(start_sub_kmer, node);

        // With max_length=2500 and k=7, expected_linear_nodes = 2493
        // MIN_EXTENSION_FRACTION = 0.1, so min_nodes = 249
        // We only have ~2 nodes, which is well below 249 and terminated
        let params = make_params();

        evaluate_seeds(&mut graph, &node_lookup, &fkc, &params, 1, &[]);

        assert!(
            graph[node].is_terminal,
            "Small dead-end seed should be abandoned"
        );
    }

    #[test]
    fn test_dual_role_seed_skipped() {
        let k: usize = 7;
        let kc = KmerCounts::new(&k);
        let fkc = kc.filtered_view(1);

        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let mut node_lookup: HashMap<u64, NodeIndex> = HashMap::new();

        let dual_sub_kmer: u64 = 42;
        let node = graph.add_node(DBNode {
            sub_kmer: dual_sub_kmer,
            is_start: true,
            is_end: true,
            is_terminal: false,
            visited: false,
        });
        node_lookup.insert(dual_sub_kmer, node);

        let params = make_params();

        evaluate_seeds(&mut graph, &node_lookup, &fkc, &params, 1, &[]);

        // Dual-role node should never be marked terminal
        assert!(
            !graph[node].is_terminal,
            "Dual-role seed should not be evaluated or abandoned"
        );
        assert!(
            !graph[node].visited,
            "Dual-role seed should not be marked visited"
        );
    }
}
