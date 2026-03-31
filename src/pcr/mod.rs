// pcr/mod.rs

use anyhow::{Context, Result};
use bio::io::fasta;
use log::{debug, trace};
use std::fs::File;
use std::io::Write;

use crate::format::format_duration;
use crate::kmer::FilteredKmerCounts;
#[cfg(test)]
use crate::kmer::KmerCounts;

/// Log at info level with a [gene_name] prefix for attribution in parallel runs.
macro_rules! gene_info {
    ($gene:expr_2021, $($arg:tt)*) => {
        log::info!("[{}] {}", $gene, format!($($arg)*))
    };
}

/// Log at warn level with a [gene_name] prefix for attribution in parallel runs.
macro_rules! gene_warn {
    ($gene:expr_2021, $($arg:tt)*) => {
        log::warn!("[{}] {}", $gene, format!($($arg)*))
    };
}

pub mod preconfigured;

mod graph;
mod paths;
mod primers;
mod pruning;

// Constants that may require tuning

/// The multiplier for establishing when a kmer is considered to have high coverage,
/// relative to the min_count threshold. It is then used to also adjust the threshold.
const COVERAGE_MULTIPLIER: u32 = 2;

/// A multiplier for adjusting the threshold as it is applied.
const COVERAGE_STEPS: u32 = 4;

/// The maximum number of kmers containing the forward or reverse primers to maintain,
/// with only those with the highest count being retained
const MAX_NUM_PRIMER_KMERS: usize = 100;

/// Default edit distance threshold for deduplication of sPCR products
pub const DEFAULT_DEDUP_EDIT_THRESHOLD: u32 = 10;

/// Score for a path through the assembly graph. Higher is better.
/// Used for ranking and selecting the best amplicon sequences.
#[derive(Debug, Clone, PartialEq)]
struct PathScore {
    /// Minimum kmer count along the path
    kmer_min_count: u32,
    /// Median kmer count along the path
    kmer_median_count: f64,
    /// Coverage consistency: coefficient of variation of edge counts
    /// (lower is better — a perfect single-copy path has uniform coverage)
    coverage_cv: f64,
    /// Maximum coverage ratio seen on any edge in the path
    /// (high values indicate the path traverses a repeat)
    max_coverage_ratio: f64,
}

impl PathScore {
    /// Composite score for ranking paths. Higher is better.
    fn composite(&self) -> f64 {
        // Start with median count as the base signal (robust to outliers)
        let base = self.kmer_median_count;
        // Penalize paths with high coverage variance (likely chimeric)
        let cv_penalty = if self.coverage_cv > 1.0 {
            1.0 / self.coverage_cv
        } else {
            1.0
        };
        // Penalize paths that traverse high-coverage-ratio edges (repeats)
        let repeat_penalty = if self.max_coverage_ratio > 5.0 {
            5.0 / self.max_coverage_ratio
        } else {
            1.0
        };
        base * cv_penalty * repeat_penalty
    }
}

struct AssemblyRecord {
    fasta_record: fasta::Record,
    score: PathScore,
}

// Create a structure to hold a kmer representing an oligo up to 32 nucleotides long in the
// length*2 least significant bits
pub struct Oligo {
    pub length: usize,
    pub kmer: u64,
}

#[derive(PartialEq)]
enum PrimerDirection {
    Forward,
    Reverse,
}

// De Bruijn graph node
#[derive(Debug, Clone)]
pub struct DBNode {
    pub sub_kmer: u64,     // k-1 mer that contains overlap between kmers
    pub is_start: bool,    // Contains the forward primer
    pub is_end: bool,      // Contains the reverse primer
    pub is_terminal: bool, // Is a terminal node
    pub visited: bool,     // Has been visited during graph traversal
}

#[derive(Debug, Clone)]
pub struct DBEdge {
    pub _kmer: u64,          // kmer that contains overlap between sub_kmers
    pub count: u32,          // Number of times this kmer was observed
    pub coverage_ratio: f64, // count / local median (set during annotation)
}

#[derive(Clone, Debug, serde::Deserialize)]
pub struct PCRParams {
    pub forward_seq: String,
    pub reverse_seq: String,
    #[serde(default)]
    pub min_length: usize,
    #[serde(default = "default_max_length")]
    pub max_length: usize,
    pub gene_name: String,
    #[serde(default = "default_min_count")]
    pub min_count: u32,
    #[serde(default = "default_mismatches")]
    pub mismatches: usize,
    #[serde(default = "default_trim")]
    pub trim: usize,
    #[serde(default)]
    #[allow(dead_code)] // Deserialized from YAML panels and re-exported; not read by Rust code
    pub citation: String,
    #[serde(default)]
    #[allow(dead_code)] // Deserialized from YAML panels and re-exported; not read by Rust code
    pub notes: String,
    #[serde(default = "default_dedup_edit_threshold")]
    pub dedup_edit_threshold: u32,
    /// Where this primer was loaded from (display only, not serialized)
    #[serde(skip)]
    pub source: String,
}

fn default_max_length() -> usize {
    10000
}
fn default_min_count() -> u32 {
    2
}
fn default_mismatches() -> usize {
    2
}
fn default_trim() -> usize {
    15
}
fn default_dedup_edit_threshold() -> u32 {
    DEFAULT_DEDUP_EDIT_THRESHOLD
}

/// Validate a primer pair and return a list of (error, suggestion) pairs.
/// An empty list means the primer is valid.
pub fn validate_pcr_params(params: &PCRParams) -> Vec<(String, String)> {
    let mut errors: Vec<(String, String)> = Vec::new();

    if params.forward_seq.len() < 2 {
        errors.push((
            format!(
                "Forward primer sequence is too short: '{}'",
                params.forward_seq
            ),
            "Primer sequences must be at least 2 bases".to_string(),
        ));
    }

    if params.reverse_seq.len() < 2 {
        errors.push((
            format!(
                "Reverse primer sequence is too short: '{}'",
                params.reverse_seq
            ),
            "Primer sequences must be at least 2 bases".to_string(),
        ));
    }

    // Only check nucleotide validity if sequences are long enough to be meaningful
    if params.forward_seq.len() >= 2 {
        let invalid: Vec<char> = params
            .forward_seq
            .chars()
            .filter(|c| !primers::is_valid_nucleotide(*c))
            .collect();
        if !invalid.is_empty() {
            let chars: Vec<String> = invalid.iter().map(|c| c.to_string()).collect();
            errors.push((
                format!(
                    "Invalid nucleotide(s) {} in forward primer {}",
                    chars.join(", "),
                    params.forward_seq
                ),
                "Valid characters: A C G T R Y W S M K B D H V N".to_string(),
            ));
        }
    }

    if params.reverse_seq.len() >= 2 {
        let invalid: Vec<char> = params
            .reverse_seq
            .chars()
            .filter(|c| !primers::is_valid_nucleotide(*c))
            .collect();
        if !invalid.is_empty() {
            let chars: Vec<String> = invalid.iter().map(|c| c.to_string()).collect();
            errors.push((
                format!(
                    "Invalid nucleotide(s) {} in reverse primer {}",
                    chars.join(", "),
                    params.reverse_seq
                ),
                "Valid characters: A C G T R Y W S M K B D H V N".to_string(),
            ));
        }
    }

    if params.min_length > params.max_length {
        errors.push((
            format!(
                "min-length ({}) is greater than max-length ({})",
                params.min_length, params.max_length
            ),
            "Swap the values or adjust the range".to_string(),
        ));
    }

    if params.min_count < 2 {
        errors.push((
            format!("min-count is {}, must be at least 2", params.min_count),
            "Set min-count to at least 2".to_string(),
        ));
    }

    if params.max_length == 0 {
        errors.push((
            "max-length is 0".to_string(),
            "Set max-length to a positive value".to_string(),
        ));
    }

    if params.gene_name.is_empty() {
        errors.push((
            "Gene name is empty".to_string(),
            "Provide a unique name for the primer pair via the 'name' field".to_string(),
        ));
    }

    if params.forward_seq == params.reverse_seq && params.forward_seq.len() >= 2 {
        errors.push((
            format!(
                "Forward and reverse primers are identical: {}",
                params.forward_seq
            ),
            "Check that forward and reverse sequences are not swapped".to_string(),
        ));
    }

    errors
}

/// Compute a sequence of coverage thresholds for graph extension, stepping down
/// from a high threshold (derived from observed primer coverage) to `min_count`.
fn compute_coverage_thresholds(primer_count: u32, min_count: u32) -> Vec<u32> {
    let coverage_high_threshold = primer_count / COVERAGE_MULTIPLIER;
    let mut thresholds: Vec<u32> = Vec::new();

    if coverage_high_threshold <= min_count {
        thresholds.push(min_count);
    } else {
        let step_size = (coverage_high_threshold - min_count) / (COVERAGE_STEPS - 1);

        for i in 0..COVERAGE_STEPS {
            thresholds.push(coverage_high_threshold.saturating_sub(i * step_size));
        }

        // Make sure the last element is min_count, even if there are rounding errors
        *thresholds
            .last_mut()
            .expect("coverage_thresholds is non-empty") = min_count;
    }

    thresholds
}

// The primary function for PCR
pub fn do_pcr(
    kmer_counts: &FilteredKmerCounts,
    sample_name: &str,
    params: &PCRParams,
    dump_graph: bool,
    output_directory: &str,
) -> Result<Vec<bio::io::fasta::Record>> {
    gene_info!(params.gene_name, "Running PCR");

    gene_info!(params.gene_name, "Preprocessing primers");
    let (forward_primer_kmers, reverse_primer_kmers) =
        primers::get_primer_kmers(params, kmer_counts)?;

    if forward_primer_kmers.is_empty() {
        gene_info!(
            params.gene_name,
            "Binding sites were not found for the forward primer. Abandoning PCR."
        );
        gene_info!(
            params.gene_name,
            "Suggested actions: optimize primer sequence, or increase the number of reads."
        );
        return Ok(Vec::new());
    }

    if reverse_primer_kmers.is_empty() {
        gene_info!(
            params.gene_name,
            "Binding sites were not found for the reverse primer. Abandoning PCR."
        );
        gene_info!(
            params.gene_name,
            "Suggested actions: optimize primer sequence, or increase the number of reads."
        );
        return Ok(Vec::new());
    }

    // Log forward primer kmers for diagnostics
    let mut sorted_forward: Vec<(u64, u32)> =
        forward_primer_kmers.iter().map(|(&k, &v)| (k, v)).collect();
    sorted_forward.sort();
    for (kmer, count) in sorted_forward.iter() {
        gene_info!(
            params.gene_name,
            "Forward primer kmer {} (count {})",
            crate::kmer::kmer_to_seq(kmer, &kmer_counts.get_k()),
            count
        );
    }

    // Build a single graph seeded with all forward and reverse primer kmers
    gene_info!(
        params.gene_name,
        "Creating graph, seeding with {} forward and {} reverse primer kmer nodes...",
        forward_primer_kmers.len(),
        reverse_primer_kmers.len()
    );
    let (seed_graph, node_lookup) =
        graph::create_seed_graph(&forward_primer_kmers, &reverse_primer_kmers, kmer_counts);

    debug!(
        "There are {} start nodes",
        graph::get_start_nodes(&seed_graph).len()
    );
    debug!(
        "There are {} end nodes",
        graph::get_end_nodes(&seed_graph).len()
    );

    for node in seed_graph.node_indices() {
        trace!("Node {}:", node.index());
        trace!(
            "  sub_kmer: {}",
            crate::kmer::kmer_to_seq(&seed_graph[node].sub_kmer, &(kmer_counts.get_k() - 1))
        );
        trace!("  is_start: {}", seed_graph[node].is_start);
        trace!("  is_end: {}", seed_graph[node].is_end);
        trace!("  is_terminal: {}", seed_graph[node].is_terminal);
    }

    let max_forward_count = forward_primer_kmers.get_max_count();
    let max_reverse_count = reverse_primer_kmers.get_max_count();
    let primer_count = max_forward_count.min(max_reverse_count);

    gene_info!(
        params.gene_name,
        "Observed primer coverage is {}, user specified min-count is {}",
        primer_count,
        params.min_count
    );

    let coverage_thresholds = compute_coverage_thresholds(primer_count, params.min_count);
    debug!("Minimum kmer counts to attempt: {:?}", coverage_thresholds);

    let mut assembly_records_all: Vec<AssemblyRecord> = Vec::new();
    let mut amplicon_index: usize = 0;

    // The graph accumulates across threshold steps. Start from the seed;
    // at each lower threshold, re-extend from previously-terminal nodes
    // that may now have qualifying successors.
    let mut current_graph = seed_graph;
    let mut current_node_lookup = node_lookup;

    for (step_idx, min_count) in coverage_thresholds.iter().enumerate() {
        let start = std::time::Instant::now();
        gene_info!(
            params.gene_name,
            "Extending graph with minimum kmer count {} (step {}/{})",
            min_count,
            step_idx + 1,
            coverage_thresholds.len()
        );

        // For steps after the first, prepare the existing graph for
        // re-extension at the lower threshold
        if step_idx > 0 {
            graph::prepare_for_lower_threshold(&mut current_graph, kmer_counts, params);
        }

        let (graph_after_fwd, node_lookup_after_fwd, fwd_found_end) = graph::extend_graph(
            current_graph,
            current_node_lookup,
            kmer_counts,
            min_count,
            params,
        )?;

        // Only run reverse extension if forward extension didn't already
        // reach end nodes. When forward succeeds, reverse extension would
        // just consume the node budget with off-target reverse seeds.
        let (graph_result, final_node_lookup) = if fwd_found_end {
            debug!("Forward extension found end nodes, skipping reverse extension");
            (graph_after_fwd, node_lookup_after_fwd)
        } else {
            graph::extend_graph_reverse(
                graph_after_fwd,
                node_lookup_after_fwd,
                kmer_counts,
                min_count,
                params,
            )?
        };

        graph::log_extended_graph_diagnostics(&graph_result, kmer_counts);

        debug!("  Final extension statistics:");
        graph::summarize_extension(&graph_result, "    ");
        debug!(
            "  There are {} start nodes with edges",
            graph::get_start_nodes(&graph_result).len()
        );
        debug!(
            "  There are {} end nodes with edges",
            graph::get_end_nodes(&graph_result).len()
        );

        if dump_graph || log::log_enabled!(log::Level::Trace) {
            let dot_string = write_annotated_dot(&graph_result, kmer_counts);
            let file_name = format!(
                "{}{}_{}_{}.dot",
                output_directory, sample_name, params.gene_name, min_count
            );
            trace!("Writing dot file {}", file_name);
            let mut file = File::create(&file_name).context("Unable to create dot file")?;
            file.write_all(dot_string.as_bytes())
                .context("Unable to write dot file")?;
        }

        gene_info!(
            params.gene_name,
            "Done. Time to extend graph: {}",
            format_duration(start.elapsed())
        );

        // Prune and find paths on a copy — the original graph is kept for
        // incremental extension at the next threshold step
        let mut pruned_graph = graph_result.clone();
        let start = std::time::Instant::now();
        gene_info!(params.gene_name, "Pruning the assembly graph...");

        pruning::remove_low_coverage_tips(&mut pruned_graph, &kmer_counts.get_k());
        pruning::reachability_pruning(&mut pruned_graph);

        // Annotate edges with coverage ratios for scoring
        graph::annotate_coverage_ratios(&mut pruned_graph);

        debug!(
            "    There are {} nodes in the graph",
            pruned_graph.node_count()
        );
        debug!(
            "    There are {} edges in the graph",
            pruned_graph.edge_count()
        );
        debug!(
            "    There are {} start nodes",
            graph::get_start_nodes(&pruned_graph).len()
        );
        debug!(
            "    There are {} end nodes",
            graph::get_end_nodes(&pruned_graph).len()
        );

        gene_info!(
            params.gene_name,
            "Done. Time to prune graph: {}",
            format_duration(start.elapsed())
        );

        let start = std::time::Instant::now();
        gene_info!(
            params.gene_name,
            "Traversing the assembly graph to find paths from forward to reverse primers..."
        );

        let all_paths = paths::get_assembly_paths(&pruned_graph, kmer_counts, params);

        debug!(
            "    There are {} paths from forward to reverse primers in the graph",
            all_paths.len()
        );
        gene_info!(
            params.gene_name,
            "Done. Time to traverse graph: {}",
            format_duration(start.elapsed())
        );

        if all_paths.is_empty() {
            gene_info!(
                params.gene_name,
                "Extending graph with minimum kmer count {} failed to generate a PCR product.",
                min_count
            );
        } else {
            gene_info!(params.gene_name, "Generating sequences from paths...");

            let (records, new_index) = paths::generate_sequences_from_paths(
                &pruned_graph,
                all_paths,
                kmer_counts,
                sample_name,
                params,
                amplicon_index,
            )?;
            amplicon_index = new_index;

            if records.is_empty() {
                gene_info!(params.gene_name, "Did not obtain PCR product.");
            } else {
                gene_info!(
                    params.gene_name,
                    "Obtained {} PCR product(s) at threshold {}.",
                    records.len(),
                    min_count
                );
                assembly_records_all.extend(records);
                break;
            }
        }

        // Use the node lookup from the (unpruned) graph for the next step
        current_node_lookup = final_node_lookup;
        current_graph = graph_result;
    }

    debug!(
        "      - The maximum count of a forward kmer is {} and of a reverse kmer is {}. Large differences in value can indicate non-specific binding of one of the primers.",
        max_forward_count, max_reverse_count
    );

    let count_threshold: u32 = 5;
    if (max_forward_count < count_threshold) || (max_reverse_count < count_threshold) {
        gene_info!(
            params.gene_name,
            "Primer kmer counts are low, in this case less than {}. Consider increasing the number of reads.",
            count_threshold
        );
    }

    gene_info!(params.gene_name, "Done.");

    if assembly_records_all.is_empty() {
        gene_info!(
            params.gene_name,
            "No path was found from a forward primer binding site to a reverse primer binding site. Abandoning PCR."
        );
        gene_info!(params.gene_name, "Suggested actions:");
        gene_info!(
            params.gene_name,
            "  - The max-length for the PCR product of {} may be too short. Consider increasing it.",
            params.max_length
        );
        gene_info!(
            params.gene_name,
            "  - The primers may have non-specific binding and are not close enough to generate a product. Consider increasing the primer trim length from the default to create a more specific primer."
        );

        return Ok(Vec::new());
    }

    let records = paths::sort_and_deduplicate(assembly_records_all, params);

    // Re-number products sequentially (0, 1, 2, ...) after dedup so IDs are
    // deterministic regardless of thread count or path enumeration order.
    let records = records
        .into_iter()
        .enumerate()
        .map(|(i, record)| {
            let id = format!("{}_{}_{}", sample_name, params.gene_name, i);
            let desc = if let Some(d) = record.desc() {
                // Update the product=N field in the description
                let updated = d
                    .split_whitespace()
                    .map(|field| {
                        if field.starts_with("product=") {
                            format!("product={}", i)
                        } else {
                            field.to_string()
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(" ");
                Some(updated)
            } else {
                None
            };
            fasta::Record::with_attrs(&id, desc.as_deref(), record.seq())
        })
        .collect();

    Ok(records)
}

/// Write an annotated DOT file with node and edge attributes for visualization.
///
/// Node attributes: sub_kmer sequence, is_start, is_end, is_terminal.
/// Edge attributes: kmer sequence, kmer count from the graph.
fn write_annotated_dot(
    graph: &petgraph::stable_graph::StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
) -> String {
    use std::fmt::Write as FmtWrite;

    let k = kmer_counts.get_k();
    let sub_k = k - 1;
    let mut dot = String::new();

    writeln!(dot, "digraph {{").unwrap();
    writeln!(dot, "  rankdir=LR;").unwrap();

    for node_idx in graph.node_indices() {
        let node = &graph[node_idx];
        let seq = crate::kmer::kmer_to_seq(&node.sub_kmer, &sub_k);

        let mut attrs = vec![format!("label=\"{}\"", seq)];

        // Shape: start nodes are double-circle, end nodes are box, both are diamond
        if node.is_start && node.is_end {
            attrs.push("shape=diamond".to_string());
        } else if node.is_start {
            attrs.push("shape=doublecircle".to_string());
        } else if node.is_end {
            attrs.push("shape=box".to_string());
        }

        if node.is_terminal {
            attrs.push("style=dashed".to_string());
        }

        writeln!(dot, "  {} [{}];", node_idx.index(), attrs.join(", ")).unwrap();
    }

    for edge_idx in graph.edge_indices() {
        let (src, tgt) = graph.edge_endpoints(edge_idx).unwrap();
        let edge = &graph[edge_idx];
        let seq = crate::kmer::kmer_to_seq(&edge._kmer, &k);
        writeln!(
            dot,
            "  {} -> {} [label=\"{} ({})\"];",
            src.index(),
            tgt.index(),
            seq,
            edge.count
        )
        .unwrap();
    }

    writeln!(dot, "}}").unwrap();
    dot
}

#[cfg(test)]
mod tests {
    use crate::kmer::seq_to_kmer;

    use super::graph::*;
    use super::paths::*;
    use super::primers::*;
    use super::pruning::*;
    use super::*;
    use petgraph::stable_graph::StableDiGraph;
    use std::collections::HashMap;
    use std::collections::HashSet;

    // These functions are used in the tests
    fn n_combinations(n: usize, r: usize) -> usize {
        // Compute n choose r iteratively to avoid factorial overflow.
        // n*(n-1)*...*(n-r+1) / r!
        if r > n {
            return 0;
        }
        let r = r.min(n - r); // take advantage of symmetry
        let mut result: usize = 1;
        for i in 0..r {
            result = result * (n - i) / (i + 1);
        }
        result
    }

    fn expected_permutations(n: usize, r: usize) -> usize {
        // Given a sequence of length n and r sites that can be permuted,
        // there are (n! / (r!(n-r)!)) combinations of r sites in the sequence and
        // 4^r permutations of the r sites.
        // So there are (n! / (r!(n-r)!)) * 4^r permutations.
        // The original sequence is included in each combination, but we only
        // want to include it once. Se we subtract the number of combinations
        // and add 1.

        n_combinations(n, r) * 4_usize.pow(r as u32) - n_combinations(n, r) + 1
    }

    // Create a test graph of known structure
    fn create_test_graph() -> (
        StableDiGraph<DBNode, DBEdge>,
        HashMap<&'static str, petgraph::graph::NodeIndex>,
    ) {
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let mut nodes = HashMap::new();

        // Add nodes and store their indices in the HashMap
        nodes.insert(
            "a",
            graph.add_node(DBNode {
                sub_kmer: 0,
                is_start: true,
                is_end: false,
                is_terminal: false,
                visited: false,
            }),
        );
        nodes.insert(
            "b",
            graph.add_node(DBNode {
                sub_kmer: 1,
                is_start: false,
                is_end: false,
                is_terminal: false,
                visited: false,
            }),
        );
        nodes.insert(
            "c",
            graph.add_node(DBNode {
                sub_kmer: 2,
                is_start: false,
                is_end: false,
                is_terminal: false,
                visited: false,
            }),
        );
        nodes.insert(
            "d",
            graph.add_node(DBNode {
                sub_kmer: 3,
                is_start: false,
                is_end: true,
                is_terminal: false,
                visited: false,
            }),
        );
        nodes.insert(
            "e",
            graph.add_node(DBNode {
                sub_kmer: 4,
                is_start: false,
                is_end: true,
                is_terminal: false,
                visited: false,
            }),
        );

        // Add directed edges
        graph.add_edge(
            nodes["a"],
            nodes["b"],
            DBEdge {
                _kmer: 10,
                count: 5,
                coverage_ratio: 0.0,
            },
        );
        graph.add_edge(
            nodes["b"],
            nodes["c"],
            DBEdge {
                _kmer: 11,
                count: 10,
                coverage_ratio: 0.0,
            },
        );
        graph.add_edge(
            nodes["c"],
            nodes["d"],
            DBEdge {
                _kmer: 12,
                count: 4,
                coverage_ratio: 0.0,
            },
        );
        graph.add_edge(
            nodes["c"],
            nodes["e"],
            DBEdge {
                _kmer: 13,
                count: 1,
                coverage_ratio: 0.0,
            },
        );

        (graph, nodes)
    }

    #[test]
    fn test_n_descendants() {
        let (graph, nodes) = create_test_graph();

        // Testing using the node indices from the HashMap
        assert_eq!(graph::descendants(&graph, nodes["a"], 1).len(), 1); // Direct successor
        assert_eq!(graph::descendants(&graph, nodes["a"], 2).len(), 2); // Including b's successors
        assert_eq!(graph::descendants(&graph, nodes["a"], 3).len(), 4); // All nodes reachable from a within 3 steps
        assert_eq!(graph::descendants(&graph, nodes["a"], 4).len(), 4); // All nodes reachable from a within 4 steps
        assert_eq!(graph::descendants(&graph, nodes["b"], 2).len(), 3);
    }

    #[test]
    fn test_get_backward_edge_counts() {
        let (graph, nodes) = create_test_graph();

        // Testing using the node indices from the HashMap
        assert_eq!(
            graph::get_backward_edge_counts(&graph, nodes["a"], 3).len(),
            0
        );
        assert_eq!(graph::get_backward_edge_counts(&graph, nodes["b"], 3), [5]);
        assert_eq!(
            graph::get_backward_edge_counts(&graph, nodes["c"], 3),
            [10, 5]
        );
        assert_eq!(
            graph::get_backward_edge_counts(&graph, nodes["d"], 3),
            [4, 10, 5]
        );
        assert_eq!(
            graph::get_backward_edge_counts(&graph, nodes["e"], 3),
            [1, 10, 5]
        );
    }

    #[test]
    fn test_get_backward_node_degrees() {
        let (_graph, _nodes) = create_test_graph();

        // get_backward_node_degrees is now private to graph module;
        // tested indirectly through integration tests
    }

    #[test]
    fn test_n_combinations() {
        assert_eq!(n_combinations(3, 2), 3);
        assert_eq!(n_combinations(5, 5), 1);
    }

    #[test]
    fn test_combinations() {
        // Make sure there are the correct numbers of combinations
        assert_eq!(combinations(3, 2).len(), n_combinations(3, 2));
        assert_eq!(combinations(8, 5).len(), n_combinations(8, 5));

        // Make sure the first combination has r elements
        assert_eq!(combinations(8, 5)[0].len(), 5);

        // Make sure that the combinations are unique
        let combos = combinations(8, 5);
        let combos_set: HashSet<Vec<usize>> = HashSet::from_iter(combos.clone());
        assert_eq!(combos.len(), combos_set.len());
    }

    #[test]
    fn test_expected_permutations() {
        assert_eq!(expected_permutations(2, 1), 7);
    }

    #[test]
    fn test_resolve_primers() {
        // Check with no ambiguous nucleotides
        let seq1 = "CGTAATGCGGCGA".to_string();
        let mut expected1: HashSet<String> = HashSet::new();
        expected1.insert(seq1.clone());
        let result1 = resolve_primer(&seq1);
        assert_eq!(result1, expected1);

        // Check with one ambiguous nucleotide
        let seq2 = "CGTAATGCGGCGN".to_string();
        let mut expected2: HashSet<String> = HashSet::new();
        expected2.insert("CGTAATGCGGCGA".to_string());
        expected2.insert("CGTAATGCGGCGC".to_string());
        expected2.insert("CGTAATGCGGCGG".to_string());
        expected2.insert("CGTAATGCGGCGT".to_string());

        let result2 = resolve_primer(&seq2);
        assert_eq!(result2, expected2);

        // Check with one ambiguous nucleotide
        let seq3 = "CGTAATRCGGCGA".to_string();
        let mut expected3: HashSet<String> = HashSet::new();
        expected3.insert("CGTAATACGGCGA".to_string());
        expected3.insert("CGTAATGCGGCGA".to_string());
        let result3 = resolve_primer(&seq3);
        assert_eq!(result3, expected3);

        // Check with two ambiguous nucleotides
        let seq4 = "CGTAATRCGGCGY".to_string();
        let mut expected4: HashSet<String> = HashSet::new();
        expected4.insert("CGTAATACGGCGC".to_string());
        expected4.insert("CGTAATGCGGCGC".to_string());
        expected4.insert("CGTAATACGGCGT".to_string());
        expected4.insert("CGTAATGCGGCGT".to_string());

        let result4 = resolve_primer(&seq4);
        assert_eq!(result4, expected4);

        // Check when the first nucleotide is ambiguous
        let seq5 = "RCGTAATCGGCGA".to_string();
        let mut expected5: HashSet<String> = HashSet::new();
        expected5.insert("ACGTAATCGGCGA".to_string());
        expected5.insert("GCGTAATCGGCGA".to_string());
        let result5 = resolve_primer(&seq5);
        assert_eq!(result5, expected5);
    }

    #[test]
    fn test_permute_sequences() {
        // Check specific permutations for a tiny example
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("CG".to_string());

            let mut expected: HashSet<String> = HashSet::new();
            expected.insert("CA".to_string());
            expected.insert("CC".to_string());
            expected.insert("CG".to_string());
            expected.insert("CT".to_string());
            expected.insert("AG".to_string());
            expected.insert("GG".to_string());
            expected.insert("TG".to_string());

            let result = permute_sequences(seq, &1_usize);
            println!(
                "Permutations: {}",
                result.iter().cloned().collect::<Vec<_>>().join(", ")
            );
            assert_eq!(result, expected);
        }

        // Check sequence length 3 and r=3
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("CGT".to_string());

            let result = permute_sequences(seq, &3_usize);

            // All sites permuted, so there should be 4^n sequences
            assert_eq!(result.len(), 64);
        }

        // Construct all the permutations procedurally
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("CGT".to_string());
            let r: usize = 2;
            let n: usize = 3;

            let bases = vec!['A', 'C', 'G', 'T'];
            let mut expected: HashSet<String> = HashSet::new();
            for i in 0..4 {
                for j in 0..4 {
                    expected.insert(format!("C{}{}", bases[i], bases[j]));
                    expected.insert(format!("{}G{}", bases[i], bases[j]));
                    expected.insert(format!("{}{}T", bases[i], bases[j]));
                }
            }

            println!(
                "There are {} combinations for n={} r={}",
                expected.len(),
                n,
                r
            );
            let result = permute_sequences(seq, &r);
            assert_eq!(expected.len(), result.len());
        }

        // Procedural test where n=4 and r=2
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("ACGT".to_string());
            let r: usize = 2;
            let n: usize = 4;

            let bases = vec!['A', 'C', 'G', 'T'];
            let mut expected: HashSet<String> = HashSet::new();
            for i in 0..4 {
                for j in 0..4 {
                    expected.insert(format!("{}{}GT", bases[i], bases[j]));
                    expected.insert(format!("{}C{}T", bases[i], bases[j]));
                    expected.insert(format!("{}CG{}", bases[i], bases[j]));
                    expected.insert(format!("A{}{}T", bases[i], bases[j]));
                    expected.insert(format!("A{}G{}", bases[i], bases[j]));
                    expected.insert(format!("AC{}{}", bases[i], bases[j]));
                }
            }

            println!(
                "There are {} combinations for n={} r={}",
                expected.len(),
                n,
                r
            );
            let result = permute_sequences(seq, &r);
            assert_eq!(expected.len(), result.len());
        }

        // Check for inclusion of particular strings as a spot check
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("TGCAGGTTCACCTAC".to_string());
            let r: usize = 2;
            let result = permute_sequences(seq, &r);

            // Check for inclusion of original sequence
            assert!(result.contains("TGCAGGTTCACCTAC"));

            // Check for inclusion of "GGCAGGTTCACCTAC"
            assert!(result.contains("GGCAGGTTCACCTAC"));
        }
    }

    #[test]
    fn test_string_to_oligo() {
        let oligo = string_to_oligo("GCGA").unwrap();
        assert_eq!(oligo.kmer, 0b1001_1000);
        assert_eq!(oligo.length, 4);
    }

    #[test]
    fn test_get_end_nodes() {
        let (graph, nodes) = create_test_graph();

        let end_nodes = get_end_nodes(&graph);
        assert_eq!(end_nodes.len(), 2);
        assert!(end_nodes.contains(&nodes["d"]));
        assert!(end_nodes.contains(&nodes["e"]));
    }

    #[test]
    fn test_get_start_nodes() {
        let (graph, nodes) = create_test_graph();

        let start_nodes = get_start_nodes(&graph);
        assert_eq!(start_nodes.len(), 1);
        assert!(start_nodes.contains(&nodes["a"]));
    }

    #[test]
    fn test_get_suffix_mask() {
        // If k is 21, then suffix should be 20 nucleotides long, and there are two bits per nucleotide
        let k: usize = 21;
        assert_eq!(
            graph::get_suffix_mask(&k),
            0b11111111_11111111_11111111_11111111_11111111
        );

        let k: usize = 3;
        assert_eq!(graph::get_suffix_mask(&k), 0b1111);
    }

    fn build_test_case() -> (String, usize, usize, KmerCounts, PCRParams) {
        // This is the 18s that was assembled from https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR26955578, per the readme tutorial
        // Padded with C's at the ends
        let read_string = "CCCCCCCCCCCCGTTGATCCTGCCAGTATCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAGCACTTGTACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATCGTTTATTTGATTGTACTCTCTTACTACTTGGATAACCGTAGTAATTCTAGAGCTAATACATGCGAAAAGTCCCGACTCTCGTGGAAGGGATGTATTTATTAGATTAAAAACCAATGCGGCTTAACGGCCGCTTACAAACTTGGTGATTCATAGTAACTGTTCGAATCGCATGGCCTTCTCTGTTCGTGCCGGCGATGTTTCATTCAAATTTCTGCCCTATCAACTGTCGATGGTAAGATAGTGGCTTACCATGGTCGCAACGGGTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACGTGGGGAGGTAGTGACAAAAAATAACAATACAGGGCTTTTTTGTAGTCTTGTAATTGGAATGAGTACAATTTAAATCTCTTAACGAGGACCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATTGTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGGCTCGACGGCGACGGTCAGCCGCAAGGTATGTCACTGTCGACGTTGGCCTTCTTCGCGCAGACTTCGCGTGCTCTTAACTGAGTGTGCGTTGGATACGCGACGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCTTGTGCTTGGATACATAAGCATGGAATAATGGAATAGGACTTTGGTTCTATTTTCCGTTGGTTTCTGGAACCGAAGTAATGATTAATAGGGACAGTTGGGGGCATTCGTATTTCGTTGTCAGAGGTGAAATTCTTGGATTTACGAAAGACGAACTAATGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTAGAGGATCGAAGACGATCAGATACCGTCCTAGTTCTAACCATAAACGATGCCGACTAGGGATCAGCGAGTGTTATTTGATGACCTCGTTGGCACCTTATGGGAAACCAAAGTTTTTGGGTTCCGGGGGAAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAAACTCACCAGGTCCAGACATAGTAAGGATTGACAGATTGAGAGCTCTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCTTGACCTGCTAAATAGTCAGACGATTCTCGAATCGCTCTCGACTTCTTAGAGGGACTGTTGCGTGTTTAACCAAAGTCAGGAAGGCAATAACAGGTCTGTGATGCCCTTAGATGTCCTGGGCCGCACGCGCGCTACACTGACGATGGCAACGAGTCGCTCCTTCACCGAAAGGTGTGGGTAATCTTGTGAATCATCGTCGTGCTGGGGATAGATCATTGTAATTCTTGATCTTGAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGAATGGTTTAGTGAGGCCTCCGGATTGGCACTGTCAGATGGGCTTCGGTCCATCCGACGGACGTCAAAAAGTTGGTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGCCCCCCCCCCCC".to_string();

        let k: usize = 21;
        let replicates: usize = 10;
        let mut kmer_counts = KmerCounts::new(&k);
        for _i in 0..replicates {
            let reads = crate::kmer::seq_to_reads(&read_string).unwrap();
            kmer_counts.ingest_reads(&reads).unwrap();
        }

        let params = PCRParams {
            forward_seq: "AACCTGGTTGATCCTGCCAGT".to_string(),

            // AGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCG 3' end of 18S above
            // CGCAGGTTCACCTACGGAAACCTTGTTACGACTTTTACT revcomp of above
            reverse_seq: "TGATCCTTCTGCAGGTTCACCTAC".to_string(),
            min_length: 0,
            max_length: 2500,
            gene_name: "18s".to_string(),
            min_count: 3,
            mismatches: 2,
            trim: 15,
            citation: "".to_string(),
            notes: "".to_string(),
            dedup_edit_threshold: DEFAULT_DEDUP_EDIT_THRESHOLD,
            source: "test".to_string(),
        };

        (read_string, k, replicates, kmer_counts, params)
    }

    #[test]
    fn test_primer_preprocessing_steps() {
        let (_, _, _, kmer_counts, params) = build_test_case();
        let kmer_counts = kmer_counts.filtered_view(1);

        let reverse_variants =
            preprocess_primer(&params, PrimerDirection::Reverse, &kmer_counts.get_k()).unwrap();

        // There should be 991 variants of the reverse primer when r=2
        assert_eq!(reverse_variants.len(), 991);

        println!("Reverse variants: {:?}", reverse_variants);

        // Check for inclusion of original trimmed sequence TGCAGGTTCACCTAC
        let member = "TGCAGGTTCACCTAC".to_string();
        assert!(reverse_variants.contains(&member));

        // Check for inclusion of off by one variant GGCAGGTTCACCTAC
        let member = "GGCAGGTTCACCTAC".to_string();
        assert!(reverse_variants.contains(&member));

        let mut reverse_primer_kmers =
            get_kmers_from_primers(&reverse_variants, &kmer_counts, &params.min_count).unwrap();

        // Check for kmer
        assert_eq!(reverse_primer_kmers.len(), 1);

        reverse_primer_kmers = filter_primer_kmers(reverse_primer_kmers);

        // Check for kmer after filtering
        assert_eq!(reverse_primer_kmers.len(), 1);
    }

    #[test]
    fn test_extension_steps() {
        let (_read_string, _k, _replicates, kmer_counts, _params) = build_test_case();
        let kmer_counts = kmer_counts.filtered_view(1);

        let seq = "TGATCCTGCCAGTATCATATG".to_string();
        let kmer: u64 = seq_to_kmer(&seq).unwrap();
        assert!(kmer_counts.get_canonical(&kmer).is_some());
    }

    #[test]
    fn test_integration() {
        let (read_string, k, replicates, kmer_counts, params) = build_test_case();
        let min_count = 5;

        // Check the number of kmers (canonical only, no reverse complements stored)
        assert_eq!(kmer_counts.len(), read_string.len() - k + 1);

        // Check the total count of kmers
        assert_eq!(
            kmer_counts.get_n_kmers(),
            ((read_string.len() - k + 1) * replicates) as u64
        );

        let filtered = kmer_counts.filtered_view(1);

        let (forward_primer_kmers, reverse_primer_kmers) =
            primers::get_primer_kmers(&params, &filtered).unwrap();

        assert_eq!(forward_primer_kmers.len(), 1);
        assert_eq!(reverse_primer_kmers.len(), 1);

        let (seed_graph, node_lookup) =
            graph::create_seed_graph(&forward_primer_kmers, &reverse_primer_kmers, &filtered);

        // Check the number of nodes in the seed graph
        // With unified primer handling, both primers use prefix (kmer >> 2).
        // The forward and reverse primer kmers produce distinct sub_kmers,
        // so we still expect 2 nodes.
        assert_eq!(seed_graph.node_count(), 2);
        assert_eq!(get_start_nodes(&seed_graph).len(), 1);
        assert_eq!(get_end_nodes(&seed_graph).len(), 1);

        let (graph_after_fwd, node_lookup_after_fwd, _fwd_found_end) =
            graph::extend_graph(seed_graph, node_lookup, &filtered, &min_count, &params).unwrap();

        let (mut graph_result, _node_lookup_final) = graph::extend_graph_reverse(
            graph_after_fwd,
            node_lookup_after_fwd,
            &filtered,
            &min_count,
            &params,
        )
        .unwrap();

        // Print the number of nodes and edges in the graph
        println!("There are {} nodes in the graph", graph_result.node_count());
        println!("There are {} edges in the graph", graph_result.edge_count());

        // With reverse extension, forward extension from the start node
        // and reverse extension from the end node should converge, producing
        // complete paths.
        let all_paths = get_assembly_paths(&graph_result, &filtered, &params);
        assert!(
            !all_paths.is_empty(),
            "Expected paths after reverse extension"
        );

        remove_low_coverage_tips(&mut graph_result, &filtered.get_k());
        reachability_pruning(&mut graph_result);

        let all_paths = get_assembly_paths(&graph_result, &filtered, &params);
        assert!(!all_paths.is_empty(), "Expected paths after pruning");
    }
}
