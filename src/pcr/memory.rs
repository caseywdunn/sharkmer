//! Owned graph backing buffers charged before allocation and cloning.
//!
//! This owner covers primary graph arrays; `lookup` and `buffers` charge
//! node maps and extension/pruning/repeat workspaces to the same ceiling.
//! Primer candidates, DFS/path buffers, read threading and bubble resolution
//! also share this ceiling. Returned sequence/diagnostic/results use separate
//! parent ownership; other library/process headroom still needs accounting.

use super::{DBEdge, DBNode};
use crate::kmer::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use anyhow::{Context, Result, ensure};
use petgraph::graph::{Edge, EdgeIndex, Node, NodeIndex};
use petgraph::stable_graph::StableDiGraph;
use std::ops::{Deref, Index, IndexMut};
use std::sync::Arc;

pub(crate) mod admission;
pub(crate) mod buffers;
pub(crate) mod lookup;

type Graph = StableDiGraph<DBNode, DBEdge>;

/// A stable graph whose backing-buffer reservation outlives its allocation.
/// Read access is borrowed; growing mutable graph access never escapes.
pub(super) struct AssemblyGraph {
    // Field order closes the allocations before their reservation on all paths.
    graph: Graph,
    reservation: Option<MemoryReservation>,
    charged_bytes: u64,
    budget: Option<Arc<ResourceBudget>>,
}

impl AssemblyGraph {
    /// Construct an empty graph, optionally sharing a real assembly budget.
    pub(super) fn new(budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        let reservation = budget
            .map(|budget| budget.reserve(MemoryClass::TargetGraphs, 0))
            .transpose()?;
        Ok(Self {
            graph: Graph::new(),
            reservation,
            charged_bytes: 0,
            budget: budget.cloned(),
        })
    }

    /// Borrow this graph's target ceiling for auxiliary owners.
    pub(super) fn budget(&self) -> Option<&Arc<ResourceBudget>> {
        self.budget.as_ref()
    }

    /// Immutable access for generic petgraph algorithms.
    pub(super) fn as_graph(&self) -> &Graph {
        &self.graph
    }

    fn bytes(nodes: usize, edges: usize) -> Result<u64> {
        // petgraph 0.6.5 StableGraph stores Graph<Option<N>, Option<E>>.
        let nodes = nodes
            .checked_mul(size_of::<Node<Option<DBNode>>>())
            .context("Graph node allocation size overflow")?;
        let edges = edges
            .checked_mul(size_of::<Edge<Option<DBEdge>>>())
            .context("Graph edge allocation size overflow")?;
        u64::try_from(
            nodes
                .checked_add(edges)
                .context("Graph allocation size overflow")?,
        )
        .context("Graph allocation size exceeds the supported range")
    }

    /// Live array capacity, excluding allocator/process and traversal overhead.
    pub(super) fn buffer_bytes(&self) -> Result<u64> {
        let (nodes, edges) = self.graph.capacity();
        Self::bytes(nodes, edges)
    }

    fn reserve_growth(&mut self, nodes: bool) -> Result<Option<(usize, usize)>> {
        if self.reservation.is_none() {
            return Ok(None);
        }
        let (node_capacity, edge_capacity) = self.graph.capacity();
        ensure!(
            self.buffer_bytes()? == self.charged_bytes,
            "Graph allocation differs from its reservation"
        );
        let (count, capacity) = if nodes {
            (self.graph.node_count(), node_capacity)
        } else {
            (self.graph.edge_count(), edge_capacity)
        };
        if count < capacity {
            // A vacant slot or spare capacity can be used without allocation.
            return Ok(None);
        }
        ensure!(
            count < NodeIndex::<u32>::end().index(),
            "Graph index range exhausted"
        );
        // For these nonzero element sizes, Vec grows from four slots and then
        // doubles. Charge both old and replacement buffers before petgraph's
        // push; check the observed capacity before accepting the operation.
        let new_capacity = capacity
            .checked_mul(2)
            .context("Graph capacity overflow")?
            .max(4);
        let replacement = if nodes {
            Self::bytes(new_capacity, 0)?
        } else {
            Self::bytes(0, new_capacity)?
        };
        let peak = self
            .charged_bytes
            .checked_add(replacement)
            .context("Graph resize peak overflow")?;
        self.reservation
            .as_mut()
            .context("Graph reservation is missing")?
            .resize(peak)?;
        Ok(Some(if nodes {
            (new_capacity, edge_capacity)
        } else {
            (node_capacity, new_capacity)
        }))
    }

    fn finish_growth(&mut self, expected: Option<(usize, usize)>) -> Result<()> {
        let Some((nodes, edges)) = expected else {
            return Ok(());
        };
        let actual = self.graph.capacity();
        ensure!(
            actual.0 <= nodes && actual.1 <= edges,
            "Graph growth exceeded its audited capacity bound"
        );
        let bytes = self.buffer_bytes()?;
        self.reservation
            .as_mut()
            .context("Graph reservation is missing")?
            .resize(bytes)?;
        self.charged_bytes = bytes;
        Ok(())
    }

    /// Add one node, checking old/new simultaneous allocation before growth.
    pub(super) fn add_node(&mut self, node: DBNode) -> Result<NodeIndex> {
        let growth = self.reserve_growth(true)?;
        let index = self.graph.add_node(node);
        self.finish_growth(growth)?;
        Ok(index)
    }

    /// Add one edge; vacant slots need no new backing allocation.
    pub(super) fn add_edge(
        &mut self,
        source: NodeIndex,
        target: NodeIndex,
        edge: DBEdge,
    ) -> Result<EdgeIndex> {
        let growth = self.reserve_growth(false)?;
        let index = self.graph.add_edge(source, target, edge);
        self.finish_growth(growth)?;
        Ok(index)
    }

    /// Clone into independently charged buffers, preserving holes and IDs.
    pub(super) fn try_clone(&self) -> Result<Self> {
        let budget = self.budget.as_ref();
        let mut reservation = budget
            .map(|budget| budget.reserve(MemoryClass::TargetGraphs, self.buffer_bytes()?))
            .transpose()?;
        let graph = self.graph.clone();
        let (nodes, edges) = graph.capacity();
        let bytes = Self::bytes(nodes, edges)?;
        ensure!(
            bytes <= self.buffer_bytes()?,
            "Graph clone exceeded its source capacity bound"
        );
        if let Some(charge) = &mut reservation {
            charge.resize(bytes)?;
        }
        Ok(Self {
            graph,
            reservation,
            charged_bytes: bytes,
            budget: budget.cloned(),
        })
    }

    /// Removal retains the backing capacity and its ownership charge.
    pub(super) fn remove_node(&mut self, node: NodeIndex) -> Option<DBNode> {
        self.graph.remove_node(node)
    }

    /// Removal retains the backing capacity and its ownership charge.
    pub(super) fn remove_edge(&mut self, edge: EdgeIndex) -> Option<DBEdge> {
        self.graph.remove_edge(edge)
    }

    /// Existing pruning only removes nodes/edges and changes scalar weights.
    pub(super) fn prune_tips(&mut self, k: &usize, tip_coverage_fraction: f64) -> Result<()> {
        super::pruning::remove_low_coverage_tips_budgeted(
            &mut self.graph,
            k,
            tip_coverage_fraction,
            self.budget.as_ref(),
        )
    }

    /// Existing reachability pruning only removes nodes and incident edges.
    pub(super) fn prune_reachability(&mut self) -> Result<()> {
        super::pruning::reachability_pruning_budgeted(&mut self.graph, self.budget.as_ref())
    }

    /// Coverage annotation changes scalar edge weights without graph growth.
    pub(super) fn annotate_coverage(&mut self) -> Result<()> {
        super::graph::annotate_coverage_ratios(&mut self.graph, self.budget.as_ref())
    }
}

impl From<Graph> for AssemblyGraph {
    fn from(graph: Graph) -> Self {
        // Explicit unbudgeted adapter for legacy tests/callers. A graph that
        // has already allocated cannot be adopted as precharged ownership.
        Self {
            graph,
            reservation: None,
            charged_bytes: 0,
            budget: None,
        }
    }
}

impl Deref for AssemblyGraph {
    type Target = Graph;
    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl Index<NodeIndex> for AssemblyGraph {
    type Output = DBNode;
    fn index(&self, index: NodeIndex) -> &Self::Output {
        &self.graph[index]
    }
}

impl IndexMut<NodeIndex> for AssemblyGraph {
    fn index_mut(&mut self, index: NodeIndex) -> &mut Self::Output {
        &mut self.graph[index]
    }
}

impl Index<EdgeIndex> for AssemblyGraph {
    type Output = DBEdge;
    fn index(&self, index: EdgeIndex) -> &Self::Output {
        &self.graph[index]
    }
}

impl IndexMut<EdgeIndex> for AssemblyGraph {
    fn index_mut(&mut self, index: EdgeIndex) -> &mut Self::Output {
        &mut self.graph[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{ResourceExhausted, ResourceKind, ResourceLimits};

    fn resources(total: u64, graph: u64) -> Arc<ResourceBudget> {
        let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
        classes[MemoryClass::TargetGraphs as usize] = graph;
        classes[MemoryClass::CountSlots as usize] = total;
        ResourceBudget::new(ResourceLimits {
            memory_bytes: total,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    fn node(value: u64) -> DBNode {
        DBNode {
            sub_kmer: value,
            is_start: false,
            is_end: false,
        }
    }
    fn edge() -> DBEdge {
        DBEdge {
            count: 1,
            coverage_ratio: 1.0,
        }
    }

    #[test]
    fn resize_refuses_before_mutation_and_retry_preserves_nodes_and_charge() {
        let bytes = AssemblyGraph::bytes(4, 0).unwrap();
        let budget = resources(4 * bytes, 4 * bytes);
        let mut graph = AssemblyGraph::new(Some(&budget)).unwrap();
        for value in 0..4 {
            graph.add_node(node(value)).unwrap();
        }
        let other = budget
            .reserve(MemoryClass::TargetGraphs, 2 * bytes)
            .unwrap();
        let before = budget.snapshot();
        let error = graph.add_node(node(4)).unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::StageMemory(MemoryClass::TargetGraphs)
        );
        assert_eq!(graph.node_count(), 4);
        assert_eq!(graph.capacity().0, 4);
        assert_eq!(graph[NodeIndex::new(3)].sub_kmer, 3);
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            before.reserved_memory_bytes
        );
        assert_eq!(
            budget.snapshot().peak_reserved_memory_bytes,
            before.peak_reserved_memory_bytes
        );
        drop(other);
        graph.add_node(node(4)).unwrap();
        assert_eq!(graph.buffer_bytes().unwrap(), 2 * bytes);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 2 * bytes);
        assert_eq!(budget.snapshot().peak_reserved_memory_bytes, 3 * bytes);
        drop(graph);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn clones_and_reused_holes_preserve_stable_ids_and_independent_ownership() {
        let budget = resources(1 << 20, 1 << 20);
        let mut graph = AssemblyGraph::new(Some(&budget)).unwrap();
        let nodes: Vec<_> = (0..17)
            .map(|value| graph.add_node(node(value)).unwrap())
            .collect();
        let edges: Vec<_> = nodes
            .windows(2)
            .map(|pair| graph.add_edge(pair[0], pair[1], edge()).unwrap())
            .collect();
        graph.remove_node(nodes[3]);
        graph.remove_edge(edges[8]);
        let bytes = graph.buffer_bytes().unwrap();
        let mut cloned = graph.try_clone().unwrap();
        assert!(!cloned.contains_node(nodes[3]));
        assert_eq!(cloned[nodes[16]].sub_kmer, 16);
        assert!(cloned.edge_weight(edges[8]).is_none());
        assert_eq!(cloned.add_node(node(42)).unwrap(), nodes[3]);
        cloned[nodes[16]].is_end = true;
        assert!(!graph[nodes[16]].is_end);
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            bytes + cloned.buffer_bytes().unwrap()
        );
        drop(graph);
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            cloned.buffer_bytes().unwrap()
        );
        drop(cloned);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn graph_and_counter_reservations_share_the_same_run_limit() {
        let bytes = AssemblyGraph::bytes(4, 0).unwrap();
        let budget = resources(4 * bytes, 4 * bytes);
        let counts = budget.reserve(MemoryClass::CountSlots, 2 * bytes).unwrap();
        let mut graph = AssemblyGraph::new(Some(&budget)).unwrap();
        for value in 0..4 {
            graph.add_node(node(value)).unwrap();
        }
        let error = graph.add_node(node(4)).unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::RunMemory
        );
        assert_eq!(graph.node_count(), 4);
        drop(counts);
        graph.add_node(node(4)).unwrap();
        drop(graph);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn rejected_clone_preserves_the_original_and_retry_uses_the_same_ledger() {
        let bytes = AssemblyGraph::bytes(4, 0).unwrap();
        let budget = resources(4 * bytes, 4 * bytes);
        let mut graph = AssemblyGraph::new(Some(&budget)).unwrap();
        for value in 0..4 {
            graph.add_node(node(value)).unwrap();
        }
        let other = budget
            .reserve(MemoryClass::TargetGraphs, 3 * bytes)
            .unwrap();
        let before = budget.snapshot();
        let error = graph.try_clone().err().unwrap();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::StageMemory(MemoryClass::TargetGraphs)
        );
        assert_eq!(graph.node_count(), 4);
        assert_eq!(graph[NodeIndex::new(3)].sub_kmer, 3);
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            before.reserved_memory_bytes
        );
        assert_eq!(
            budget.snapshot().peak_reserved_memory_bytes,
            before.peak_reserved_memory_bytes
        );
        drop(other);
        let cloned = graph.try_clone().unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 2 * bytes);
        drop(graph);
        drop(cloned);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }
}
