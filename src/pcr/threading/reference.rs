//! Frozen pre-ownership threading algorithm for differential tests only.
#![allow(dead_code)]

use ahash::{AHashMap, AHashSet};
use petgraph::Direction;
use petgraph::graph::EdgeIndex;
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};

use super::{DBEdge, DBNode};
use crate::io::{Mate, ReadRecord};
use crate::kmer::encoding::revcomp_kmer;

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

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PairedEndLink {
    pub r1_edges: Vec<EdgeIndex>,
    pub r2_edges: Vec<EdgeIndex>,
}

#[derive(Debug, Clone)]
pub struct ThreadingAnnotations {
    pub edge_support: AHashMap<EdgeIndex, EdgeReadSupport>,
    pub branch_links: AHashMap<BranchLink, u32>,
    pub paired_links: Vec<PairedEndLink>,
}

impl ThreadingAnnotations {
    fn new() -> Self {
        Self {
            edge_support: AHashMap::new(),
            branch_links: AHashMap::new(),
            paired_links: Vec::new(),
        }
    }
}

#[derive(Clone, Copy)]
struct PositionedKmer {
    kmer: u64,
    start: usize,
}

#[derive(Clone)]
struct ReadRun {
    edges: Vec<EdgeIndex>,
}

#[derive(Clone, Default)]
struct ReadEvidence {
    total_edges: AHashSet<EdgeIndex>,
    unambiguous_edges: AHashSet<EdgeIndex>,
    branch_links: AHashSet<BranchLink>,
    ordered_edges: Vec<EdgeIndex>,
}

impl ReadEvidence {
    fn merge(&mut self, other: Self) {
        self.total_edges.extend(other.total_edges);
        self.unambiguous_edges.extend(other.unambiguous_edges);
        self.branch_links.extend(other.branch_links);
        let mut seen: AHashSet<EdgeIndex> = self.ordered_edges.iter().copied().collect();
        for edge in other.ordered_edges {
            if seen.insert(edge) {
                self.ordered_edges.push(edge);
            }
        }
    }

    fn has_same_local_support(&self, other: &Self) -> bool {
        self.total_edges == other.total_edges
            && self.unambiguous_edges == other.unambiguous_edges
            && self.branch_links == other.branch_links
            && self.ordered_edges == other.ordered_edges
    }
}

#[derive(Default)]
struct PairEvidence {
    r1: Option<ReadEvidence>,
    r2: Option<ReadEvidence>,
}

pub fn thread_reads(
    graph: &StableDiGraph<DBNode, DBEdge>,
    reads: &[ReadRecord],
    k: usize,
) -> ThreadingAnnotations {
    let mut annotations = ThreadingAnnotations::new();
    let edge_lookup = build_edge_lookup(graph);

    for read in reads {
        let evidence = map_read(graph, &edge_lookup, &read.sequence, k);
        apply_evidence(&mut annotations, &evidence);
    }

    annotations
}

pub fn thread_reads_paired(
    graph: &StableDiGraph<DBNode, DBEdge>,
    reads: &[ReadRecord],
    k: usize,
) -> ThreadingAnnotations {
    let mut annotations = ThreadingAnnotations::new();
    let edge_lookup = build_edge_lookup(graph);
    let mut pairs: AHashMap<u64, PairEvidence> = AHashMap::new();

    for read in reads {
        let evidence = map_read(graph, &edge_lookup, &read.sequence, k);
        if evidence.total_edges.is_empty() {
            continue;
        }
        match read.mate {
            Mate::Unpaired => apply_evidence(&mut annotations, &evidence),
            Mate::R1 => merge_mate(&mut pairs.entry(read.index / 2).or_default().r1, evidence),
            Mate::R2 => merge_mate(&mut pairs.entry(read.index / 2).or_default().r2, evidence),
        }
    }

    for pair in pairs.into_values() {
        let r1_edges = pair
            .r1
            .as_ref()
            .map(|evidence| evidence.ordered_edges.clone())
            .unwrap_or_default();
        let r2_edges = pair
            .r2
            .as_ref()
            .map(|evidence| evidence.ordered_edges.clone())
            .unwrap_or_default();
        let mut fragment_evidence = pair.r1.unwrap_or_default();
        if let Some(r2_evidence) = pair.r2 {
            fragment_evidence.merge(r2_evidence);
        }
        apply_evidence(&mut annotations, &fragment_evidence);
        if !r1_edges.is_empty() && !r2_edges.is_empty() {
            annotations
                .paired_links
                .push(PairedEndLink { r1_edges, r2_edges });
        }
    }

    annotations
}

/// Stream verified original records, retaining at most one pending fragment.
pub(crate) fn thread_replay(
    graph: &StableDiGraph<DBNode, DBEdge>,
    replay: &crate::replay::ReplayStore,
    k: usize,
) -> anyhow::Result<ThreadingAnnotations> {
    let mut annotations = ThreadingAnnotations::new();
    let edge_lookup = build_edge_lookup(graph);
    let mut pending = PairEvidence::default();
    let mut fragment = None;
    replay.visit(&mut |record| {
        if replay.paired() && fragment != Some(record.fragment) {
            apply_pair(&mut annotations, std::mem::take(&mut pending));
            fragment = Some(record.fragment);
        }
        let sequence = std::str::from_utf8(record.fields[1])?;
        let evidence = map_read(graph, &edge_lookup, sequence, k);
        if !evidence.total_edges.is_empty() {
            match record.mate {
                Mate::Unpaired => apply_evidence(&mut annotations, &evidence),
                Mate::R1 => merge_mate(&mut pending.r1, evidence),
                Mate::R2 => merge_mate(&mut pending.r2, evidence),
            }
        }
        Ok(())
    })?;
    apply_pair(&mut annotations, pending);
    Ok(annotations)
}

fn apply_pair(annotations: &mut ThreadingAnnotations, pair: PairEvidence) {
    let r1_edges = pair
        .r1
        .as_ref()
        .map(|evidence| evidence.ordered_edges.clone())
        .unwrap_or_default();
    let r2_edges = pair
        .r2
        .as_ref()
        .map(|evidence| evidence.ordered_edges.clone())
        .unwrap_or_default();
    let mut evidence = pair.r1.unwrap_or_default();
    if let Some(r2) = pair.r2 {
        evidence.merge(r2);
    }
    apply_evidence(annotations, &evidence);
    if !r1_edges.is_empty() && !r2_edges.is_empty() {
        annotations
            .paired_links
            .push(PairedEndLink { r1_edges, r2_edges });
    }
}

fn merge_mate(destination: &mut Option<ReadEvidence>, evidence: ReadEvidence) {
    if let Some(existing) = destination {
        existing.merge(evidence);
    } else {
        *destination = Some(evidence);
    }
}

fn build_edge_lookup(graph: &StableDiGraph<DBNode, DBEdge>) -> AHashMap<u64, EdgeIndex> {
    let mut lookup = AHashMap::new();
    for edge in graph.edge_references() {
        let kmer = crate::pcr::graph::reconstruct_edge_kmer(graph, edge.id());
        let previous = lookup.insert(kmer, edge.id());
        debug_assert!(previous.is_none(), "directional graph kmers must be unique");
    }
    lookup
}

fn positioned_directional_kmers(sequence: &str, k: usize) -> Vec<PositionedKmer> {
    if k == 0 || k >= 32 {
        return Vec::new();
    }
    let mask = (1u64 << (2 * k)) - 1;
    let mut frame = 0u64;
    let mut valid_bases = 0usize;
    let mut windows = Vec::with_capacity(sequence.len().saturating_sub(k - 1));

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
            });
        }
    }

    windows
}

fn reverse_normalized_kmers(
    windows: &[PositionedKmer],
    sequence_length: usize,
    k: usize,
) -> Vec<PositionedKmer> {
    windows
        .iter()
        .rev()
        .map(|window| PositionedKmer {
            kmer: revcomp_kmer(&window.kmer, &k),
            start: sequence_length - k - window.start,
        })
        .collect()
}

fn map_read(
    graph: &StableDiGraph<DBNode, DBEdge>,
    edge_lookup: &AHashMap<u64, EdgeIndex>,
    sequence: &str,
    k: usize,
) -> ReadEvidence {
    let forward_windows = positioned_directional_kmers(sequence, k);
    let reverse_windows = reverse_normalized_kmers(&forward_windows, sequence.len(), k);
    let forward_runs = find_contiguous_runs(&forward_windows, edge_lookup, graph);
    let reverse_runs = find_contiguous_runs(&reverse_windows, edge_lookup, graph);
    let forward_span = longest_run(&forward_runs);
    let reverse_span = longest_run(&reverse_runs);

    match forward_span.cmp(&reverse_span) {
        std::cmp::Ordering::Greater => evidence_from_runs(graph, &forward_runs),
        std::cmp::Ordering::Less => evidence_from_runs(graph, &reverse_runs),
        std::cmp::Ordering::Equal => {
            let forward_evidence = evidence_from_runs(graph, &forward_runs);
            let reverse_evidence = evidence_from_runs(graph, &reverse_runs);
            if forward_evidence.has_same_local_support(&reverse_evidence) {
                forward_evidence
            } else {
                common_evidence(forward_evidence, reverse_evidence)
            }
        }
    }
}

fn longest_run(runs: &[ReadRun]) -> usize {
    runs.iter().map(|run| run.edges.len()).max().unwrap_or(0)
}

fn find_contiguous_runs(
    windows: &[PositionedKmer],
    edge_lookup: &AHashMap<u64, EdgeIndex>,
    graph: &StableDiGraph<DBNode, DBEdge>,
) -> Vec<ReadRun> {
    let mut runs = Vec::new();
    let mut current_edges = Vec::new();
    let mut previous_start = None;

    for window in windows {
        let Some(&edge) = edge_lookup.get(&window.kmer) else {
            flush_run(&mut runs, &mut current_edges);
            previous_start = None;
            continue;
        };

        let position_is_contiguous = previous_start.is_some_and(|start| window.start == start + 1);
        let graph_is_contiguous = current_edges.last().is_some_and(|previous_edge| {
            let (_, previous_target) = graph
                .edge_endpoints(*previous_edge)
                .expect("edge must exist in graph");
            let (current_source, _) = graph
                .edge_endpoints(edge)
                .expect("edge must exist in graph");
            previous_target == current_source
        });

        if !current_edges.is_empty() && !(position_is_contiguous && graph_is_contiguous) {
            flush_run(&mut runs, &mut current_edges);
        }
        current_edges.push(edge);
        previous_start = Some(window.start);
    }

    flush_run(&mut runs, &mut current_edges);
    runs
}

fn flush_run(runs: &mut Vec<ReadRun>, current_edges: &mut Vec<EdgeIndex>) {
    if !current_edges.is_empty() {
        runs.push(ReadRun {
            edges: std::mem::take(current_edges),
        });
    }
}

fn evidence_from_runs(graph: &StableDiGraph<DBNode, DBEdge>, runs: &[ReadRun]) -> ReadEvidence {
    let mut evidence = ReadEvidence::default();
    let mut ordered_seen = AHashSet::new();

    for run in runs {
        let is_unambiguous = is_run_unambiguous(graph, &run.edges);
        for &edge in &run.edges {
            evidence.total_edges.insert(edge);
            if is_unambiguous {
                evidence.unambiguous_edges.insert(edge);
            }
            if ordered_seen.insert(edge) {
                evidence.ordered_edges.push(edge);
            }
        }
        evidence
            .branch_links
            .extend(branch_links_for_run(graph, &run.edges));
    }

    evidence
}

fn common_evidence(forward: ReadEvidence, reverse: ReadEvidence) -> ReadEvidence {
    let total_edges: AHashSet<EdgeIndex> = forward
        .total_edges
        .intersection(&reverse.total_edges)
        .copied()
        .collect();
    let branch_links = forward
        .branch_links
        .intersection(&reverse.branch_links)
        .copied()
        .collect();
    let ordered_edges = forward
        .ordered_edges
        .into_iter()
        .filter(|edge| total_edges.contains(edge))
        .collect();
    ReadEvidence {
        total_edges,
        unambiguous_edges: AHashSet::new(),
        branch_links,
        ordered_edges,
    }
}

fn apply_evidence(annotations: &mut ThreadingAnnotations, evidence: &ReadEvidence) {
    for &edge in &evidence.total_edges {
        annotations
            .edge_support
            .entry(edge)
            .or_default()
            .read_support_total += 1;
    }
    for &edge in &evidence.unambiguous_edges {
        annotations
            .edge_support
            .entry(edge)
            .or_default()
            .read_support_unambiguous += 1;
    }
    for &link in &evidence.branch_links {
        *annotations.branch_links.entry(link).or_insert(0) += 1;
    }
}

fn is_run_unambiguous(graph: &StableDiGraph<DBNode, DBEdge>, edges: &[EdgeIndex]) -> bool {
    if edges.len() < 2 {
        return true;
    }
    edges.windows(2).all(|window| {
        let (_, node) = graph.edge_endpoints(window[0]).expect("edge must exist");
        graph.neighbors_directed(node, Direction::Incoming).count() <= 1
            && graph.neighbors_directed(node, Direction::Outgoing).count() <= 1
    })
}

fn branch_links_for_run(
    graph: &StableDiGraph<DBNode, DBEdge>,
    edges: &[EdgeIndex],
) -> AHashSet<BranchLink> {
    let mut links = AHashSet::new();
    for window in edges.windows(2) {
        let incoming_edge = window[0];
        let outgoing_edge = window[1];
        let (_, node) = graph
            .edge_endpoints(incoming_edge)
            .expect("edge must exist");
        let in_degree = graph.neighbors_directed(node, Direction::Incoming).count();
        let out_degree = graph.neighbors_directed(node, Direction::Outgoing).count();
        if in_degree > 1 || out_degree > 1 {
            links.insert(BranchLink {
                incoming_edge,
                outgoing_edge,
            });
        }
    }
    links
}
