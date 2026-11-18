use anyhow::{Context, Result};
use petgraph::Direction;
use petgraph::graph::EdgeIndex;
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use std::sync::Arc;

use super::memory::buffers::BudgetSet;
use super::{DBEdge, DBNode};
use crate::io::{Mate, ReadRecord};
use crate::kmer::budget::ResourceBudget;
use crate::kmer::budget::collections::{BudgetMap, BudgetVec};
use crate::kmer::encoding::revcomp_kmer;

#[cfg(test)]
pub(super) mod reference;
#[cfg(test)]
use ahash::AHashMap;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct EdgeReadSupport {
    pub read_support_total: u32,
    pub read_support_unambiguous: u32,
}

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct BranchLink {
    pub incoming_edge: EdgeIndex,
    pub outgoing_edge: EdgeIndex,
}

#[cfg_attr(test, derive(Debug))]
#[allow(dead_code)]
pub struct PairedEndLink {
    pub r1_edges: BudgetVec<EdgeIndex>,
    pub r2_edges: BudgetVec<EdgeIndex>,
}

#[cfg_attr(test, derive(Debug))]
pub struct ThreadingAnnotations {
    pub edge_support: BudgetMap<EdgeIndex, EdgeReadSupport>,
    pub branch_links: BudgetMap<BranchLink, u32>,
    pub paired_links: BudgetVec<PairedEndLink>,
}

impl ThreadingAnnotations {
    pub(super) fn new(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self {
            edge_support: BudgetMap::new(budget),
            branch_links: BudgetMap::new(budget),
            paired_links: BudgetVec::new(budget),
        }
    }
}

#[derive(Clone, Copy)]
struct PositionedKmer {
    kmer: u64,
    start: usize,
}

struct ReadRun {
    edges: BudgetVec<EdgeIndex>,
}

struct ReadEvidence {
    total_edges: BudgetSet<EdgeIndex>,
    unambiguous_edges: BudgetSet<EdgeIndex>,
    branch_links: BudgetSet<BranchLink>,
    ordered_edges: BudgetVec<EdgeIndex>,
}

impl ReadEvidence {
    fn new(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self {
            total_edges: BudgetSet::new(budget),
            unambiguous_edges: BudgetSet::new(budget),
            branch_links: BudgetSet::new(budget),
            ordered_edges: BudgetVec::new(budget),
        }
    }

    fn merge(&mut self, other: Self) -> Result<()> {
        for &edge in other.total_edges.iter() {
            self.total_edges.insert(edge)?;
        }
        for &edge in other.unambiguous_edges.iter() {
            self.unambiguous_edges.insert(edge)?;
        }
        for &link in other.branch_links.iter() {
            self.branch_links.insert(link)?;
        }
        let mut seen = BudgetSet::new(self.ordered_edges.budget());
        for &edge in self.ordered_edges.iter() {
            seen.insert(edge)?;
        }
        for edge in other.ordered_edges {
            if seen.insert(edge)? {
                self.ordered_edges.push(edge)?;
            }
        }
        Ok(())
    }

    fn has_same_local_support(&self, other: &Self) -> bool {
        fn same<T: Eq + std::hash::Hash>(a: &BudgetSet<T>, b: &BudgetSet<T>) -> bool {
            a.len() == b.len() && a.iter().all(|key| b.contains(key))
        }
        same(&self.total_edges, &other.total_edges)
            && same(&self.unambiguous_edges, &other.unambiguous_edges)
            && same(&self.branch_links, &other.branch_links)
            && *self.ordered_edges == *other.ordered_edges
    }
}

#[derive(Default)]
struct PairEvidence {
    r1: Option<ReadEvidence>,
    r2: Option<ReadEvidence>,
}

#[cfg(test)]
pub fn thread_reads(
    graph: &StableDiGraph<DBNode, DBEdge>,
    reads: &[ReadRecord],
    k: usize,
) -> ThreadingAnnotations {
    thread_reads_budgeted(graph, reads, k, None).unwrap()
}

/// Map independent reads under the graph's target allowance.
pub(super) fn thread_reads_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    reads: &[ReadRecord],
    k: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<ThreadingAnnotations> {
    let mut annotations = ThreadingAnnotations::new(budget);
    let edge_lookup = build_edge_lookup(graph, budget)?;
    for read in reads {
        let evidence = map_read(graph, &edge_lookup, &read.sequence, k)?;
        apply_evidence(&mut annotations, &evidence)?;
    }
    Ok(annotations)
}

#[cfg(test)]
pub fn thread_reads_paired(
    graph: &StableDiGraph<DBNode, DBEdge>,
    reads: &[ReadRecord],
    k: usize,
) -> ThreadingAnnotations {
    thread_reads_paired_budgeted(graph, reads, k, None).unwrap()
}

/// Preserve arbitrary legacy fragment order and duplicate-mate merging.
pub(super) fn thread_reads_paired_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    reads: &[ReadRecord],
    k: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<ThreadingAnnotations> {
    let mut annotations = ThreadingAnnotations::new(budget);
    let edge_lookup = build_edge_lookup(graph, budget)?;
    let mut pairs = BudgetMap::new(budget);
    for read in reads {
        let evidence = map_read(graph, &edge_lookup, &read.sequence, k)?;
        if evidence.total_edges.is_empty() {
            continue;
        }
        if read.mate == Mate::Unpaired {
            apply_evidence(&mut annotations, &evidence)?;
            continue;
        }
        let fragment = read.index / 2;
        if pairs.get(&fragment).is_none() {
            pairs.insert(fragment, PairEvidence::default())?;
        }
        let pair = pairs
            .get_mut(&fragment)
            .context("Missing fragment evidence after insertion")?;
        if read.mate == Mate::R1 {
            merge_mate(&mut pair.r1, evidence)?;
        } else {
            merge_mate(&mut pair.r2, evidence)?;
        }
    }
    for (_, pair) in pairs {
        apply_pair(&mut annotations, pair)?;
    }
    Ok(annotations)
}

#[cfg(test)]
pub(crate) fn thread_replay(
    graph: &StableDiGraph<DBNode, DBEdge>,
    replay: &crate::replay::ReplayStore,
    k: usize,
) -> Result<ThreadingAnnotations> {
    thread_replay_budgeted(graph, replay, k, None)
}

/// Stream verified original records, retaining at most one pending fragment.
pub(super) fn thread_replay_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    replay: &crate::replay::ReplayStore,
    k: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<ThreadingAnnotations> {
    let mut annotations = ThreadingAnnotations::new(budget);
    let edge_lookup = build_edge_lookup(graph, budget)?;
    let mut pending = PairEvidence::default();
    let mut fragment = None;
    replay.visit(&mut |record| {
        if replay.paired() && fragment != Some(record.fragment) {
            apply_pair(&mut annotations, std::mem::take(&mut pending))?;
            fragment = Some(record.fragment);
        }
        let sequence = std::str::from_utf8(record.fields[1])?;
        let evidence = map_read(graph, &edge_lookup, sequence, k)?;
        if !evidence.total_edges.is_empty() {
            match record.mate {
                Mate::Unpaired => apply_evidence(&mut annotations, &evidence)?,
                Mate::R1 => merge_mate(&mut pending.r1, evidence)?,
                Mate::R2 => merge_mate(&mut pending.r2, evidence)?,
            }
        }
        Ok(())
    })?;
    apply_pair(&mut annotations, pending)?;
    Ok(annotations)
}

fn apply_pair(annotations: &mut ThreadingAnnotations, pair: PairEvidence) -> Result<()> {
    let budget = annotations.edge_support.budget().cloned();
    let copy = |evidence: &Option<ReadEvidence>| -> Result<BudgetVec<EdgeIndex>> {
        match evidence {
            Some(evidence) => evidence.ordered_edges.try_clone(),
            None => Ok(BudgetVec::new(budget.as_ref())),
        }
    };
    let r1_edges = copy(&pair.r1)?;
    let r2_edges = copy(&pair.r2)?;
    let mut evidence = pair
        .r1
        .unwrap_or_else(|| ReadEvidence::new(budget.as_ref()));
    if let Some(r2) = pair.r2 {
        evidence.merge(r2)?;
    }
    apply_evidence(annotations, &evidence)?;
    if !r1_edges.is_empty() && !r2_edges.is_empty() {
        annotations
            .paired_links
            .push(PairedEndLink { r1_edges, r2_edges })?;
    }
    Ok(())
}

fn merge_mate(destination: &mut Option<ReadEvidence>, evidence: ReadEvidence) -> Result<()> {
    if let Some(existing) = destination {
        existing.merge(evidence)?;
    } else {
        *destination = Some(evidence);
    }
    Ok(())
}

fn build_edge_lookup(
    graph: &StableDiGraph<DBNode, DBEdge>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetMap<u64, EdgeIndex>> {
    let mut lookup = BudgetMap::new(budget);
    for edge in graph.edge_references() {
        let kmer = super::graph::reconstruct_edge_kmer(graph, edge.id());
        anyhow::ensure!(
            lookup.get(&kmer).is_none(),
            "Directional graph kmers must be unique"
        );
        lookup.insert(kmer, edge.id())?;
    }
    Ok(lookup)
}

#[cfg(test)]
fn positioned_directional_kmers(sequence: &str, k: usize) -> BudgetVec<PositionedKmer> {
    positioned_directional_kmers_budgeted(sequence, k, None).unwrap()
}

fn positioned_directional_kmers_budgeted(
    sequence: &str,
    k: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<PositionedKmer>> {
    if k == 0 || k >= 32 {
        return Ok(BudgetVec::new(budget));
    }
    let mask = (1u64 << (2 * k)) - 1;
    let mut frame = 0u64;
    let mut valid_bases = 0usize;
    let mut windows = BudgetVec::with_capacity(sequence.len().saturating_sub(k - 1), budget)?;
    for (position, nucleotide) in sequence.bytes().enumerate() {
        let base = match nucleotide {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                frame = 0;
                valid_bases = 0;
                continue;
            }
        };
        frame = ((frame << 2) | base) & mask;
        valid_bases += 1;
        if valid_bases >= k {
            windows.push(PositionedKmer {
                kmer: frame,
                start: position + 1 - k,
            })?;
        }
    }
    Ok(windows)
}

fn reverse_normalized_kmers(
    windows: &[PositionedKmer],
    sequence_length: usize,
    k: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<PositionedKmer>> {
    BudgetVec::collect(
        windows.iter().rev().map(|window| PositionedKmer {
            kmer: revcomp_kmer(&window.kmer, &k),
            start: sequence_length - k - window.start,
        }),
        budget,
    )
}

fn map_read(
    graph: &StableDiGraph<DBNode, DBEdge>,
    edge_lookup: &BudgetMap<u64, EdgeIndex>,
    sequence: &str,
    k: usize,
) -> Result<ReadEvidence> {
    let budget = edge_lookup.budget();
    let forward_windows = positioned_directional_kmers_budgeted(sequence, k, budget)?;
    let reverse_windows = reverse_normalized_kmers(&forward_windows, sequence.len(), k, budget)?;
    let forward_runs = find_contiguous_runs(&forward_windows, edge_lookup, graph)?;
    let reverse_runs = find_contiguous_runs(&reverse_windows, edge_lookup, graph)?;
    let forward_span = longest_run(&forward_runs);
    let reverse_span = longest_run(&reverse_runs);
    match forward_span.cmp(&reverse_span) {
        std::cmp::Ordering::Greater => evidence_from_runs(graph, &forward_runs, budget),
        std::cmp::Ordering::Less => evidence_from_runs(graph, &reverse_runs, budget),
        std::cmp::Ordering::Equal => {
            let forward_evidence = evidence_from_runs(graph, &forward_runs, budget)?;
            let reverse_evidence = evidence_from_runs(graph, &reverse_runs, budget)?;
            if forward_evidence.has_same_local_support(&reverse_evidence) {
                Ok(forward_evidence)
            } else {
                common_evidence(forward_evidence, reverse_evidence, budget)
            }
        }
    }
}

fn longest_run(runs: &[ReadRun]) -> usize {
    runs.iter().map(|run| run.edges.len()).max().unwrap_or(0)
}

fn find_contiguous_runs(
    windows: &[PositionedKmer],
    edge_lookup: &BudgetMap<u64, EdgeIndex>,
    graph: &StableDiGraph<DBNode, DBEdge>,
) -> Result<BudgetVec<ReadRun>> {
    let budget = edge_lookup.budget();
    let mut runs = BudgetVec::new(budget);
    let mut current_edges = BudgetVec::new(budget);
    let mut previous_start = None;
    for window in windows {
        let Some(&edge) = edge_lookup.get(&window.kmer) else {
            flush_run(&mut runs, &mut current_edges)?;
            previous_start = None;
            continue;
        };
        let position_is_contiguous = previous_start.is_some_and(|start| window.start == start + 1);
        let graph_is_contiguous = if let Some(previous_edge) = current_edges.last() {
            let (_, previous_target) = graph
                .edge_endpoints(*previous_edge)
                .context("Missing previous threading edge")?;
            let (current_source, _) = graph
                .edge_endpoints(edge)
                .context("Missing current threading edge")?;
            previous_target == current_source
        } else {
            false
        };
        if !current_edges.is_empty() && !(position_is_contiguous && graph_is_contiguous) {
            flush_run(&mut runs, &mut current_edges)?;
        }
        current_edges.push(edge)?;
        previous_start = Some(window.start);
    }
    flush_run(&mut runs, &mut current_edges)?;
    Ok(runs)
}

fn flush_run(
    runs: &mut BudgetVec<ReadRun>,
    current_edges: &mut BudgetVec<EdgeIndex>,
) -> Result<()> {
    if !current_edges.is_empty() {
        let empty = BudgetVec::new(current_edges.budget());
        runs.push(ReadRun {
            edges: std::mem::replace(current_edges, empty),
        })?;
    }
    Ok(())
}

fn evidence_from_runs(
    graph: &StableDiGraph<DBNode, DBEdge>,
    runs: &[ReadRun],
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<ReadEvidence> {
    let mut evidence = ReadEvidence::new(budget);
    let mut ordered_seen = BudgetSet::new(budget);
    for run in runs {
        let is_unambiguous = is_run_unambiguous(graph, &run.edges)?;
        for &edge in run.edges.iter() {
            evidence.total_edges.insert(edge)?;
            if is_unambiguous {
                evidence.unambiguous_edges.insert(edge)?;
            }
            if ordered_seen.insert(edge)? {
                evidence.ordered_edges.push(edge)?;
            }
        }
        branch_links_for_run(graph, &run.edges, &mut evidence.branch_links)?;
    }
    Ok(evidence)
}

fn common_evidence(
    forward: ReadEvidence,
    reverse: ReadEvidence,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<ReadEvidence> {
    let mut common = ReadEvidence::new(budget);
    for &edge in forward.total_edges.iter() {
        if reverse.total_edges.contains(&edge) {
            common.total_edges.insert(edge)?;
        }
    }
    for &link in forward.branch_links.iter() {
        if reverse.branch_links.contains(&link) {
            common.branch_links.insert(link)?;
        }
    }
    for edge in forward.ordered_edges {
        if common.total_edges.contains(&edge) {
            common.ordered_edges.push(edge)?;
        }
    }
    Ok(common)
}

fn apply_evidence(annotations: &mut ThreadingAnnotations, evidence: &ReadEvidence) -> Result<()> {
    for &edge in evidence.total_edges.iter() {
        let mut support = annotations
            .edge_support
            .get(&edge)
            .cloned()
            .unwrap_or_default();
        support.read_support_total = support
            .read_support_total
            .checked_add(1)
            .context("Read support total overflow")?;
        annotations.edge_support.insert(edge, support)?;
    }
    for &edge in evidence.unambiguous_edges.iter() {
        let mut support = annotations
            .edge_support
            .get(&edge)
            .cloned()
            .unwrap_or_default();
        support.read_support_unambiguous = support
            .read_support_unambiguous
            .checked_add(1)
            .context("Unambiguous read support overflow")?;
        annotations.edge_support.insert(edge, support)?;
    }
    for &link in evidence.branch_links.iter() {
        let count = annotations
            .branch_links
            .get(&link)
            .copied()
            .unwrap_or(0)
            .checked_add(1)
            .context("Branch-link support overflow")?;
        annotations.branch_links.insert(link, count)?;
    }
    Ok(())
}

fn is_run_unambiguous(graph: &StableDiGraph<DBNode, DBEdge>, edges: &[EdgeIndex]) -> Result<bool> {
    if edges.len() < 2 {
        return Ok(true);
    }
    for window in edges.windows(2) {
        let (_, node) = graph
            .edge_endpoints(window[0])
            .context("Missing unambiguous-run edge")?;
        if graph.neighbors_directed(node, Direction::Incoming).count() > 1
            || graph.neighbors_directed(node, Direction::Outgoing).count() > 1
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn branch_links_for_run(
    graph: &StableDiGraph<DBNode, DBEdge>,
    edges: &[EdgeIndex],
    links: &mut BudgetSet<BranchLink>,
) -> Result<()> {
    for window in edges.windows(2) {
        let incoming_edge = window[0];
        let outgoing_edge = window[1];
        let (_, node) = graph
            .edge_endpoints(incoming_edge)
            .context("Missing branch-link edge")?;
        let in_degree = graph.neighbors_directed(node, Direction::Incoming).count();
        let out_degree = graph.neighbors_directed(node, Direction::Outgoing).count();
        if in_degree > 1 || out_degree > 1 {
            links.insert(BranchLink {
                incoming_edge,
                outgoing_edge,
            })?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graph::NodeIndex;

    fn encode(sequence: &str) -> u64 {
        crate::kmer::encoding::seq_to_kmer(sequence).unwrap()
    }

    fn add_sequence(
        graph: &mut StableDiGraph<DBNode, DBEdge>,
        nodes: &mut AHashMap<u64, NodeIndex>,
        sequence: &str,
        k: usize,
    ) -> Vec<EdgeIndex> {
        let mut edges = Vec::new();
        for start in 0..=sequence.len() - k {
            let kmer = encode(&sequence[start..start + k]);
            let source_kmer = kmer >> 2;
            let target_kmer = kmer & ((1u64 << (2 * (k - 1))) - 1);
            let source = *nodes.entry(source_kmer).or_insert_with(|| {
                graph.add_node(DBNode {
                    sub_kmer: source_kmer,
                    is_start: false,
                    is_end: false,
                })
            });
            let target = *nodes.entry(target_kmer).or_insert_with(|| {
                graph.add_node(DBNode {
                    sub_kmer: target_kmer,
                    is_start: false,
                    is_end: false,
                })
            });
            let edge = graph.find_edge(source, target).unwrap_or_else(|| {
                graph.add_edge(
                    source,
                    target,
                    DBEdge {
                        count: 10,
                        coverage_ratio: 1.0,
                    },
                )
            });
            edges.push(edge);
        }
        edges
    }

    fn graph_from_sequences(sequences: &[&str], k: usize) -> StableDiGraph<DBNode, DBEdge> {
        let mut graph = StableDiGraph::new();
        let mut nodes = AHashMap::new();
        for sequence in sequences {
            add_sequence(&mut graph, &mut nodes, sequence, k);
        }
        graph
    }

    fn read(sequence: &str, index: u64, mate: Mate) -> ReadRecord {
        ReadRecord {
            sequence: sequence.to_string(),
            index,
            mate,
        }
    }

    fn supported_kmers(
        graph: &StableDiGraph<DBNode, DBEdge>,
        annotations: &ThreadingAnnotations,
        k: usize,
    ) -> AHashMap<String, EdgeReadSupport> {
        annotations
            .edge_support
            .iter()
            .map(|(edge, support)| {
                let kmer = super::super::graph::reconstruct_edge_kmer(graph, *edge);
                (
                    crate::kmer::encoding::kmer_to_seq(&kmer, &k),
                    support.clone(),
                )
            })
            .collect()
    }

    #[test]
    fn forward_and_reverse_reads_produce_equal_graph_evidence() {
        let graph = graph_from_sequences(&["AACGATTCCG", "AACGACTCCG"], 5);
        let forward = thread_reads(&graph, &[read("AACGATTCCG", 0, Mate::Unpaired)], 5);
        let reverse = thread_reads(&graph, &[read("CGGAATCGTT", 0, Mate::Unpaired)], 5);

        assert_eq!(
            supported_kmers(&graph, &forward, 5),
            supported_kmers(&graph, &reverse, 5)
        );
        assert_eq!(forward.branch_links, reverse.branch_links);
        assert!(!forward.branch_links.is_empty());
    }

    #[test]
    fn directional_primer_edge_maps_from_both_read_strands() {
        let graph = graph_from_sequences(&["TTTGA"], 5);
        let forward = thread_reads(&graph, &[read("TTTGA", 0, Mate::Unpaired)], 5);
        let reverse = thread_reads(&graph, &[read("TCAAA", 0, Mate::Unpaired)], 5);

        assert_eq!(
            supported_kmers(&graph, &forward, 5),
            supported_kmers(&graph, &reverse, 5)
        );
        assert!(supported_kmers(&graph, &forward, 5).contains_key("TTTGA"));
    }

    #[test]
    fn unique_flanks_resolve_a_read_spanning_an_inverted_motif() {
        let sequence = "GGAACGATCGTTACC";
        let reverse_complement = "GGTAACGATCGTTCC";
        let graph = graph_from_sequences(&[sequence], 5);
        let forward = thread_reads(&graph, &[read(sequence, 0, Mate::Unpaired)], 5);
        let reverse = thread_reads(&graph, &[read(reverse_complement, 0, Mate::Unpaired)], 5);

        assert_eq!(
            supported_kmers(&graph, &forward, 5),
            supported_kmers(&graph, &reverse, 5)
        );
        assert_eq!(forward.edge_support.len(), sequence.len() - 4);
    }

    #[test]
    fn gaps_split_directional_runs() {
        let graph = graph_from_sequences(&["AACGATT", "AACGACT"], 5);
        for sequence in ["AACGANACGAT", "AACGAXACGAT"] {
            let annotations = thread_reads(&graph, &[read(sequence, 0, Mate::Unpaired)], 5);
            assert!(annotations.branch_links.is_empty());
            assert_eq!(annotations.edge_support.len(), 2);
        }
    }

    #[test]
    fn interior_read_without_primer_kmer_contributes_support() {
        let graph = graph_from_sequences(&["TTTGAAACGATTCCGAAAAA"], 5);
        let annotations = thread_reads(&graph, &[read("AACGATTCCG", 0, Mate::Unpaired)], 5);

        assert!(annotations.edge_support.len() > 1);
        assert!(!supported_kmers(&graph, &annotations, 5).contains_key("TTTGA"));
    }

    #[test]
    fn longest_directional_mapping_beats_incidental_opposite_hit() {
        let graph = graph_from_sequences(&["AACGATTCCG", "CGGAA"], 5);
        let annotations = thread_reads(&graph, &[read("AACGATTCCG", 0, Mate::Unpaired)], 5);
        let support = supported_kmers(&graph, &annotations, 5);

        assert_eq!(support.len(), 6);
        assert!(!support.contains_key("CGGAA"));
    }

    #[test]
    fn tied_opposite_paths_do_not_choose_an_arbitrary_arm() {
        for sequences in [["AACGATTCCG", "CGGAATCGTT"], ["CGGAATCGTT", "AACGATTCCG"]] {
            let graph = graph_from_sequences(&sequences, 5);
            let annotations = thread_reads(&graph, &[read("AACGATTCCG", 0, Mate::Unpaired)], 5);
            assert!(annotations.edge_support.is_empty());
            assert!(annotations.branch_links.is_empty());
        }
    }

    #[test]
    fn tied_opposite_paths_retain_only_common_local_evidence() {
        let graph = graph_from_sequences(&["GAACGATCGTTA", "TAACGATCGTTC"], 5);
        let annotations = thread_reads(&graph, &[read("GAACGATCGTTA", 0, Mate::Unpaired)], 5);
        let support = supported_kmers(&graph, &annotations, 5);

        assert_eq!(support.len(), 6);
        assert!(!support.contains_key("GAACG"));
        assert!(!support.contains_key("CGTTA"));
        assert!(!support.contains_key("TAACG"));
        assert!(!support.contains_key("CGTTC"));
        assert!(
            support
                .values()
                .all(|edge_support| edge_support.read_support_unambiguous == 0)
        );
    }

    #[test]
    fn repeated_edge_and_branch_evidence_count_once_per_read() {
        let graph = graph_from_sequences(&["AACGATT", "AACGACT"], 5);
        let annotations = thread_reads(&graph, &[read("AACGATTNAACGATT", 0, Mate::Unpaired)], 5);

        assert!(
            annotations
                .edge_support
                .values()
                .all(|support| support.read_support_total == 1)
        );
        assert_eq!(
            annotations
                .branch_links
                .values()
                .copied()
                .collect::<Vec<_>>(),
            vec![1]
        );
    }

    #[test]
    fn paired_mates_count_support_once_per_fragment() {
        let graph = graph_from_sequences(&["AACGATT", "AACGACT"], 5);
        let paired = thread_reads_paired(
            &graph,
            &[read("AACGATT", 0, Mate::R1), read("AACGATT", 1, Mate::R2)],
            5,
        );
        let unpaired = thread_reads(
            &graph,
            &[
                read("AACGATT", 0, Mate::Unpaired),
                read("AACGATT", 1, Mate::Unpaired),
            ],
            5,
        );

        assert!(paired.edge_support.values().all(
            |support| support.read_support_total == 1 && support.read_support_unambiguous <= 1
        ));
        assert!(
            unpaired
                .edge_support
                .values()
                .all(|support| support.read_support_total == 2)
        );
        assert_eq!(
            paired.branch_links.values().copied().collect::<Vec<_>>(),
            vec![1]
        );
        assert_eq!(
            unpaired.branch_links.values().copied().collect::<Vec<_>>(),
            vec![2]
        );
        assert_eq!(paired.paired_links.len(), 1);
    }

    #[test]
    fn replay_matches_slice_evidence_for_mates_reverse_strands_gaps_and_unmatched_tail() {
        use crate::kmer::budget::{MemoryClass, ResourceBudget, ResourceLimits};
        use crate::replay::{ReplayRecord, ReplayWriter};
        use sha2::{Digest, Sha256};
        let graph = graph_from_sequences(&["AACGATTCCG", "AACGACTCCG"], 5);
        let directory = tempfile::tempdir().unwrap();
        let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
        classes[MemoryClass::ReplayBuffer as usize] = 4096;
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 4096,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 1 << 20,
            temporary_open_files: 128,
        })
        .unwrap();
        for paired in [false, true] {
            let reads: Vec<_> = ["AACGATT", "AANCGATT", "CGGAATCGTT", "TTTTTTT", "AACGATT"]
                .into_iter()
                .enumerate()
                .map(|(index, sequence)| {
                    read(
                        sequence,
                        index as u64,
                        if !paired {
                            Mate::Unpaired
                        } else if index % 2 == 0 {
                            Mate::R1
                        } else {
                            Mate::R2
                        },
                    )
                })
                .collect();
            let expected = if paired {
                thread_reads_paired(&graph, &reads, 5)
            } else {
                thread_reads(&graph, &reads, 5)
            };
            let mut writer =
                ReplayWriter::new(directory.path(), &budget, 1024, paired, reads.len() as u64)
                    .unwrap();
            let mut digest = Sha256::new();
            digest.update(b"sharkmer-selected-fastq-v1\0");
            for read in &reads {
                let quality = vec![b'J'; read.sequence.len()];
                let fields = [
                    b"@original identifier".as_slice(),
                    read.sequence.as_bytes(),
                    b"+original",
                    quality.as_slice(),
                ];
                let source = "one-shot stdin";
                digest.update(read.index.to_le_bytes());
                digest.update([match read.mate {
                    Mate::Unpaired => 0,
                    Mate::R1 => 1,
                    Mate::R2 => 2,
                }]);
                for field in std::iter::once(source.as_bytes()).chain(fields) {
                    digest.update((field.len() as u64).to_le_bytes());
                    digest.update(field);
                }
                writer
                    .add(ReplayRecord {
                        index: read.index,
                        fragment: if paired { read.index / 2 } else { read.index },
                        mate: read.mate,
                        source,
                        fields,
                    })
                    .unwrap();
            }
            let replay = writer.finish(digest.finalize().into()).unwrap();
            let tiny = evidence_budget(1);
            let error = thread_replay_budgeted(&graph, &replay, 5, Some(&tiny))
                .err()
                .unwrap();
            assert!(
                error
                    .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                    .is_some()
            );
            assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
            let target = evidence_budget(1 << 20);
            let bounded = thread_replay_budgeted(&graph, &replay, 5, Some(&target)).unwrap();
            let observed = thread_replay(&graph, &replay, 5).unwrap();
            assert_eq!(evidence_signature(&bounded), evidence_signature(&observed));
            assert!(target.snapshot().reserved_memory_bytes > 0);
            drop(bounded);
            assert_eq!(target.snapshot().reserved_memory_bytes, 0);
            assert_eq!(observed.edge_support, expected.edge_support);
            assert_eq!(observed.branch_links, expected.branch_links);
            let links = |annotations: ThreadingAnnotations| {
                let mut result: Vec<_> = annotations
                    .paired_links
                    .into_iter()
                    .map(|link| {
                        (
                            link.r1_edges
                                .iter()
                                .map(|edge| edge.index())
                                .collect::<Vec<_>>(),
                            link.r2_edges
                                .iter()
                                .map(|edge| edge.index())
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect();
                result.sort();
                result
            };
            assert_eq!(links(observed), links(expected));
            replay.close().unwrap();
        }
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    }

    #[test]
    fn paired_background_reads_produce_no_fragment_evidence() {
        let graph = graph_from_sequences(&["AACGATT"], 5);
        let annotations = thread_reads_paired(
            &graph,
            &[
                read("TTTTTTT", 0, Mate::R1),
                read("CCCCCCC", 1, Mate::R2),
                read("GGGGGGG", 2, Mate::R1),
                read("TATATAT", 3, Mate::R2),
            ],
            5,
        );

        assert!(annotations.edge_support.is_empty());
        assert!(annotations.branch_links.is_empty());
        assert!(annotations.paired_links.is_empty());
    }

    #[test]
    fn positioned_windows_preserve_primer_orientation_and_gap_offsets() {
        let windows = positioned_directional_kmers("TTTGANACGAT", 5);
        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].kmer, encode("TTTGA"));
        assert_eq!(windows[0].start, 0);
        assert_eq!(windows[1].kmer, encode("ACGAT"));
        assert_eq!(windows[1].start, 6);
    }
    type EvidenceSignature = (
        std::collections::BTreeMap<usize, (u32, u32)>,
        std::collections::BTreeMap<(usize, usize), u32>,
        Vec<(Vec<usize>, Vec<usize>)>,
    );

    fn evidence_signature(annotations: &ThreadingAnnotations) -> EvidenceSignature {
        let mut pairs: Vec<_> = annotations
            .paired_links
            .iter()
            .map(|pair| {
                (
                    pair.r1_edges.iter().map(|e| e.index()).collect(),
                    pair.r2_edges.iter().map(|e| e.index()).collect(),
                )
            })
            .collect();
        pairs.sort();
        (
            annotations
                .edge_support
                .iter()
                .map(|(e, s)| {
                    (
                        e.index(),
                        (s.read_support_total, s.read_support_unambiguous),
                    )
                })
                .collect(),
            annotations
                .branch_links
                .iter()
                .map(|(l, c)| ((l.incoming_edge.index(), l.outgoing_edge.index()), *c))
                .collect(),
            pairs,
        )
    }

    fn reference_signature(annotations: &reference::ThreadingAnnotations) -> EvidenceSignature {
        let mut pairs: Vec<_> = annotations
            .paired_links
            .iter()
            .map(|pair| {
                (
                    pair.r1_edges.iter().map(|e| e.index()).collect(),
                    pair.r2_edges.iter().map(|e| e.index()).collect(),
                )
            })
            .collect();
        pairs.sort();
        (
            annotations
                .edge_support
                .iter()
                .map(|(e, s)| {
                    (
                        e.index(),
                        (s.read_support_total, s.read_support_unambiguous),
                    )
                })
                .collect(),
            annotations
                .branch_links
                .iter()
                .map(|(l, c)| ((l.incoming_edge.index(), l.outgoing_edge.index()), *c))
                .collect(),
            pairs,
        )
    }

    fn evidence_budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn owned_threading_matches_frozen_oracle_for_strands_gaps_duplicates_and_fragment_order() {
        let templates = ["AACGATTCCG", "AACGACTCCG", "GAACGATCGTTA", "TAACGATCGTTC"];
        let sequences = [
            "AACGATTCCG",
            "CGGAATCGTT",
            "AANCGATTCCG",
            "AACGATTNAACGATT",
            "GAACGATCGTTA",
            "TAACGATCGTTC",
            "NNNNNN",
            "aacgATTCCG",
            "A",
        ];
        for case in 0..36 {
            let k = [3, 5, 7][case % 3];
            let mut graph = graph_from_sequences(&templates, k);
            let hole = graph.add_node(DBNode {
                sub_kmer: u64::MAX,
                is_start: false,
                is_end: false,
            });
            graph.remove_node(hole);
            if case % 2 == 0 {
                let edge = graph.edge_indices().nth(2).unwrap();
                graph.remove_edge(edge);
            }
            let mut reads: Vec<_> = (0..27)
                .map(|index| {
                    read(
                        sequences[(index * 7 + case) % sequences.len()],
                        (index % 13) as u64,
                        match index % 3 {
                            0 => Mate::Unpaired,
                            1 => Mate::R1,
                            _ => Mate::R2,
                        },
                    )
                })
                .collect();
            reads.rotate_left(case % 27);
            for paired in [false, true] {
                let budget = evidence_budget(1 << 20);
                let actual = if paired {
                    thread_reads_paired_budgeted(&graph, &reads, k, Some(&budget))
                } else {
                    thread_reads_budgeted(&graph, &reads, k, Some(&budget))
                }
                .unwrap();
                let expected = if paired {
                    reference::thread_reads_paired(&graph, &reads, k)
                } else {
                    reference::thread_reads(&graph, &reads, k)
                };
                assert_eq!(
                    evidence_signature(&actual),
                    reference_signature(&expected),
                    "case {case}, paired {paired}"
                );
                drop(actual);
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            }
        }
    }

    #[test]
    fn threading_refusal_at_successive_allocations_releases_all_nested_owners() {
        let graph = graph_from_sequences(&["AACGATTCCG", "AACGACTCCG"], 5);
        let reads = vec![
            read("AACGATTNAACGATTCCG", 0, Mate::R1),
            read("CGGAATCGTT", 1, Mate::R2),
            read("AACGACTCCG", 2, Mate::R1),
            read("AACGACTCCG", 3, Mate::R2),
        ];
        for paired in [false, true] {
            let expected = if paired {
                reference::thread_reads_paired(&graph, &reads, 5)
            } else {
                reference::thread_reads(&graph, &reads, 5)
            };
            let (mut refused, mut succeeded) = (0, 0);
            for bytes in (1..=128).map(|n| n * 128).chain([1 << 20]) {
                let budget = evidence_budget(bytes);
                let result = if paired {
                    thread_reads_paired_budgeted(&graph, &reads, 5, Some(&budget))
                } else {
                    thread_reads_budgeted(&graph, &reads, 5, Some(&budget))
                };
                match result {
                    Ok(annotations) => {
                        assert_eq!(
                            evidence_signature(&annotations),
                            reference_signature(&expected)
                        );
                        assert!(budget.snapshot().reserved_memory_bytes > 0);
                        succeeded += 1;
                    }
                    Err(error) => {
                        assert!(
                            error
                                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                                .is_some(),
                            "{error:#}"
                        );
                        refused += 1;
                    }
                }
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0, "budget {bytes}");
            }
            assert!(refused > 4 && succeeded > 0);
        }
    }

    #[test]
    fn support_counter_overflow_is_an_error_without_wrapping_the_existing_value() {
        let budget = evidence_budget(1 << 20);
        let edge = EdgeIndex::new(0);
        let link = BranchLink {
            incoming_edge: edge,
            outgoing_edge: EdgeIndex::new(1),
        };
        for kind in 0..3 {
            let mut annotations = ThreadingAnnotations::new(Some(&budget));
            let mut evidence = ReadEvidence::new(Some(&budget));
            match kind {
                0 => {
                    annotations
                        .edge_support
                        .insert(
                            edge,
                            EdgeReadSupport {
                                read_support_total: u32::MAX,
                                read_support_unambiguous: 0,
                            },
                        )
                        .unwrap();
                    evidence.total_edges.insert(edge).unwrap();
                }
                1 => {
                    annotations
                        .edge_support
                        .insert(
                            edge,
                            EdgeReadSupport {
                                read_support_total: 0,
                                read_support_unambiguous: u32::MAX,
                            },
                        )
                        .unwrap();
                    evidence.unambiguous_edges.insert(edge).unwrap();
                }
                _ => {
                    annotations.branch_links.insert(link, u32::MAX).unwrap();
                    evidence.branch_links.insert(link).unwrap();
                }
            }
            let before = evidence_signature(&annotations);
            let error = apply_evidence(&mut annotations, &evidence).unwrap_err();
            assert!(error.to_string().contains("overflow"));
            assert!(
                error
                    .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                    .is_none()
            );
            assert_eq!(evidence_signature(&annotations), before);
            drop(evidence);
            drop(annotations);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        }
    }
}
