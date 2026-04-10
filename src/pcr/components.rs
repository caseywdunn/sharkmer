// pcr/components.rs — Component tracking for graph extension
//
// LiveComponentTracker: dynamic union-find that tracks components during
// graph extension. Components merge when extension from one seed reaches
// a node belonging to another seed's component.

use std::collections::HashMap;
use std::collections::VecDeque;

use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;

use super::seed_eval::SeedEvalResult;
use super::{DBEdge, DBNode};

// ---------------------------------------------------------------------------
// LiveComponentTracker — dynamic component tracking during graph extension
// ---------------------------------------------------------------------------

/// Extension direction for a frontier node.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ExtensionDir {
    Forward,
    Reverse,
}

/// Per-component state maintained during prioritized extension.
#[allow(dead_code)]
#[derive(Debug)]
pub(super) struct ComponentInfo {
    pub start_seeds: Vec<NodeIndex>,
    pub end_seeds: Vec<NodeIndex>,
    /// Nodes added to the graph by this component's extension.
    pub nodes_used: usize,
    /// Priority for the extension heap (higher = extend first).
    pub priority: f64,
    /// Forward frontier: nodes to extend rightward.
    pub frontier_fwd: VecDeque<NodeIndex>,
    /// Reverse frontier: nodes to extend leftward.
    pub frontier_rev: VecDeque<NodeIndex>,
    /// True once this component contains both start and end seeds.
    pub is_bridged: bool,
    /// True if the component has been abandoned (excessive branching).
    pub is_abandoned: bool,
    /// Branching events observed during the current slice.
    pub slice_branching_events: usize,
    /// Nodes processed during the current slice.
    pub slice_nodes_processed: usize,
}

#[allow(dead_code)]
impl ComponentInfo {
    /// Whether the component has any frontier nodes remaining.
    pub fn has_frontier(&self) -> bool {
        !self.frontier_fwd.is_empty() || !self.frontier_rev.is_empty()
    }

    /// Reset per-slice counters at the start of each time slice.
    pub fn reset_slice_counters(&mut self) {
        self.slice_branching_events = 0;
        self.slice_nodes_processed = 0;
    }

    /// Branching ratio for the current slice.
    pub fn slice_branching_ratio(&self) -> f64 {
        if self.slice_nodes_processed == 0 {
            return 0.0;
        }
        self.slice_branching_events as f64 / self.slice_nodes_processed as f64
    }
}

/// Live component tracker using union-find. Tracks component identity for
/// every node in the graph and merges components when extensions overlap.
#[allow(dead_code)]
pub(super) struct LiveComponentTracker {
    uf: UnionFind,
    /// Component info keyed by the current union-find root.
    /// After a merge, the old root's entry is removed and absorbed into the new root.
    info: HashMap<NodeIndex, ComponentInfo>,
    /// Maps every node to the component it was registered to. Used to find
    /// the component root for any node in O(α(n)) via `uf.find()`.
    /// All nodes (seeds + extension-created) are registered here.
    node_registered: HashMap<NodeIndex, ()>,
}

#[allow(dead_code)]
impl LiveComponentTracker {
    /// Create a tracker from the seed graph. Each seed node becomes its own
    /// singleton component. Forward seeds get `frontier_fwd` populated;
    /// reverse seeds get `frontier_rev` populated.
    pub fn from_seed_graph(
        graph: &StableDiGraph<DBNode, DBEdge>,
        forward_primer_kmers: &crate::kmer::KmerCounts,
        reverse_primer_kmers: &crate::kmer::KmerCounts,
    ) -> Self {
        let mut uf = UnionFind::new();
        let mut info: HashMap<NodeIndex, ComponentInfo> = HashMap::new();
        let mut node_registered: HashMap<NodeIndex, ()> = HashMap::new();

        for node in graph.node_indices() {
            let data = &graph[node];
            if !data.is_start && !data.is_end {
                continue; // skip non-seed nodes (shouldn't exist at this point)
            }

            uf.make_set(node);
            node_registered.insert(node, ());

            // Determine priority from primer kmer count.
            // The seed's sub_kmer was derived from a primer kmer; look up
            // the max count among matching primer kmers.
            let mut priority = 0.0f64;
            if data.is_start {
                // Forward seed: sub_kmer = kmer >> 2, so any kmer with this prefix matches
                for (kmer, count) in forward_primer_kmers.iter() {
                    if (kmer >> 2) == data.sub_kmer {
                        priority = priority.max(*count as f64);
                    }
                }
            }
            if data.is_end {
                let k = if !forward_primer_kmers.is_empty() {
                    forward_primer_kmers.get_k()
                } else {
                    reverse_primer_kmers.get_k()
                };
                let suffix_mask = (1u64 << (2 * (k - 1))) - 1;
                for (kmer, count) in reverse_primer_kmers.iter() {
                    let rc = crate::kmer::revcomp_kmer(kmer, &k);
                    if (rc & suffix_mask) == data.sub_kmer {
                        priority = priority.max(*count as f64);
                    }
                }
            }

            let mut frontier_fwd = VecDeque::new();
            let mut frontier_rev = VecDeque::new();
            if data.is_start {
                frontier_fwd.push_back(node);
            }
            if data.is_end {
                frontier_rev.push_back(node);
            }

            let is_bridged = data.is_start && data.is_end;

            info.insert(
                node,
                ComponentInfo {
                    start_seeds: if data.is_start {
                        vec![node]
                    } else {
                        Vec::new()
                    },
                    end_seeds: if data.is_end { vec![node] } else { Vec::new() },
                    nodes_used: 0,
                    priority,
                    frontier_fwd,
                    frontier_rev,
                    is_bridged,
                    is_abandoned: false,
                    slice_branching_events: 0,
                    slice_nodes_processed: 0,
                },
            );
        }

        LiveComponentTracker {
            uf,
            info,
            node_registered,
        }
    }

    /// Register a newly created graph node as belonging to the given
    /// parent component.
    pub fn register_node(&mut self, node: NodeIndex, parent: NodeIndex) {
        let root = self.uf.find(parent);
        self.uf.make_set(node);
        self.uf.union(node, root);
        self.node_registered.insert(node, ());

        // Update node count on the new root (may have changed after union)
        let new_root = self.uf.find(node);
        if new_root != root {
            // Root changed — move info from old root to new root
            if let Some(mut ci) = self.info.remove(&root) {
                ci.nodes_used += 1;
                self.info.insert(new_root, ci);
            }
        } else if let Some(ci) = self.info.get_mut(&root) {
            ci.nodes_used += 1;
        }
    }

    /// Get the component root for any registered node.
    pub fn component_of(&mut self, node: NodeIndex) -> NodeIndex {
        self.uf.find(node)
    }

    /// Merge two components when an edge is added between nodes in different
    /// components. Returns `true` if the merge created a newly bridged
    /// component (both start and end seeds).
    pub fn merge(&mut self, node_a: NodeIndex, node_b: NodeIndex) -> bool {
        let root_a = self.uf.find(node_a);
        let root_b = self.uf.find(node_b);

        if root_a == root_b {
            return false; // already same component
        }

        // Perform union
        self.uf.union(root_a, root_b);
        let new_root = self.uf.find(root_a); // either root_a or root_b

        // Determine which was absorbed
        let (winner, loser) = if new_root == root_a {
            (root_a, root_b)
        } else {
            (root_b, root_a)
        };

        // Absorb loser's info into winner's
        if let Some(loser_info) = self.info.remove(&loser) {
            if let Some(winner_info) = self.info.get_mut(&winner) {
                winner_info.start_seeds.extend(loser_info.start_seeds);
                winner_info.end_seeds.extend(loser_info.end_seeds);
                winner_info.nodes_used += loser_info.nodes_used;
                winner_info.priority = winner_info.priority.max(loser_info.priority);
                winner_info.frontier_fwd.extend(loser_info.frontier_fwd);
                winner_info.frontier_rev.extend(loser_info.frontier_rev);
                let was_bridged = winner_info.is_bridged;
                winner_info.is_bridged =
                    !winner_info.start_seeds.is_empty() && !winner_info.end_seeds.is_empty();
                // Preserve non-abandoned status if either component was active
                if !loser_info.is_abandoned {
                    winner_info.is_abandoned = false;
                }
                return !was_bridged && winner_info.is_bridged;
            }
        }
        false
    }

    /// Check whether a component is bridged (has both start and end seeds).
    pub fn is_bridged(&mut self, node: NodeIndex) -> bool {
        let root = self.uf.find(node);
        self.info.get(&root).is_some_and(|ci| ci.is_bridged)
    }

    /// Get mutable reference to a component's info by any member node.
    pub fn info_mut(&mut self, node: NodeIndex) -> Option<&mut ComponentInfo> {
        let root = self.uf.find(node);
        self.info.get_mut(&root)
    }

    /// Get immutable reference to a component's info by any member node.
    pub fn info(&mut self, node: NodeIndex) -> Option<&ComponentInfo> {
        let root = self.uf.find(node);
        self.info.get(&root)
    }

    /// Whether a node has been registered in the tracker.
    pub fn is_registered(&self, node: NodeIndex) -> bool {
        self.node_registered.contains_key(&node)
    }

    /// Iterate all component roots and their info.
    pub fn components(&self) -> impl Iterator<Item = (&NodeIndex, &ComponentInfo)> {
        self.info.iter()
    }

    /// Number of active (non-abandoned) components.
    pub fn active_component_count(&self) -> usize {
        self.info.values().filter(|ci| !ci.is_abandoned).count()
    }

    /// Abandon a component (mark it as abandoned, clear frontiers).
    pub fn abandon(&mut self, node: NodeIndex) {
        let root = self.uf.find(node);
        if let Some(ci) = self.info.get_mut(&root) {
            ci.is_abandoned = true;
            ci.frontier_fwd.clear();
            ci.frontier_rev.clear();
        }
    }

    /// Get all component roots sorted by priority (highest first).
    pub fn roots_by_priority(&self) -> Vec<NodeIndex> {
        let mut roots: Vec<(NodeIndex, f64)> = self
            .info
            .iter()
            .filter(|(_, ci)| !ci.is_abandoned && ci.has_frontier())
            .map(|(&root, ci)| (root, ci.priority))
            .collect();
        roots.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        roots.into_iter().map(|(root, _)| root).collect()
    }
}

// ---------------------------------------------------------------------------
// Legacy component API (to be removed in Phase 4)
// ---------------------------------------------------------------------------

/// A connected group of seed nodes identified during seed evaluation.
pub(super) struct SeedComponent {
    /// Forward (start) seed node indices in this component
    pub start_seeds: Vec<NodeIndex>,
    /// Reverse (end) seed node indices in this component
    pub end_seeds: Vec<NodeIndex>,
    /// Whether seed eval found a forward seed that reached a reverse seed
    pub has_connected_seeds: bool,
    /// Average branching ratio across seeds in this component (lower is better)
    pub avg_branching_ratio: f64,
    /// Sum of seed eval node counts across seeds in this component
    pub total_extension_nodes: usize,
    /// Max mean edge coverage across seeds in this component (higher is better).
    /// On-target paths explore high-coverage regions; off-target hits have low counts.
    pub max_mean_edge_count: f64,
    /// Min coverage CV across seeds in this component (lower is better).
    /// Uniform coverage = real amplicon; high CV = wandering across repeat boundaries.
    pub min_coverage_cv: f64,
    /// Allocated node budget for extension
    pub node_budget: usize,
    /// Priority score for ordering (higher = extend first)
    pub priority_score: f64,
}

impl SeedComponent {
    /// All seed node indices in this component (start + end).
    pub fn all_seeds(&self) -> Vec<NodeIndex> {
        let mut seeds = self.start_seeds.clone();
        seeds.extend_from_slice(&self.end_seeds);
        seeds
    }

    /// Number of seeds (start + end) in this component.
    pub fn seed_count(&self) -> usize {
        self.start_seeds.len() + self.end_seeds.len()
    }
}

// --- Union-Find for grouping seeds into components ---

struct UnionFind {
    parent: HashMap<NodeIndex, NodeIndex>,
    rank: HashMap<NodeIndex, usize>,
}

impl UnionFind {
    fn new() -> Self {
        UnionFind {
            parent: HashMap::new(),
            rank: HashMap::new(),
        }
    }

    fn make_set(&mut self, x: NodeIndex) {
        self.parent.entry(x).or_insert(x);
        self.rank.entry(x).or_insert(0);
    }

    fn find(&mut self, x: NodeIndex) -> NodeIndex {
        // Iterative path find: walk to root, then compress
        let mut current = x;
        while self.parent[&current] != current {
            current = self.parent[&current];
        }
        let root = current;
        // Path compression: point all nodes on the path directly to root
        current = x;
        while self.parent[&current] != root {
            let next = self.parent[&current];
            self.parent.insert(current, root);
            current = next;
        }
        root
    }

    fn union(&mut self, x: NodeIndex, y: NodeIndex) {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx == ry {
            return;
        }
        let rank_x = self.rank[&rx];
        let rank_y = self.rank[&ry];
        if rank_x < rank_y {
            self.parent.insert(rx, ry);
        } else if rank_x > rank_y {
            self.parent.insert(ry, rx);
        } else {
            self.parent.insert(ry, rx);
            self.rank.insert(rx, rank_x + 1);
        }
    }
}

/// Identify connected components among surviving seed nodes using the
/// connectivity information from seed evaluation.
///
/// Two seeds are in the same component if a forward seed reached a reverse
/// seed during bounded exploration (directly or transitively). Seeds with
/// no connections form singleton components.
pub(super) fn identify_components(
    eval: &SeedEvalResult,
    graph: &StableDiGraph<DBNode, DBEdge>,
) -> Vec<SeedComponent> {
    let mut uf = UnionFind::new();

    // Register all surviving seeds
    for (node_idx, _is_start, _metrics) in &eval.seed_metrics {
        uf.make_set(*node_idx);
    }

    // Union connected seeds
    for (fwd, rev) in &eval.connections {
        uf.make_set(*fwd);
        uf.make_set(*rev);
        uf.union(*fwd, *rev);
    }

    // Group seeds by component root
    let mut component_map: HashMap<NodeIndex, Vec<NodeIndex>> = HashMap::new();
    let surviving_nodes: Vec<NodeIndex> = eval.seed_metrics.iter().map(|(n, _, _)| *n).collect();
    for &node in &surviving_nodes {
        let root = uf.find(node);
        component_map.entry(root).or_default().push(node);
    }

    // Build SeedComponent structs
    let mut components: Vec<SeedComponent> = Vec::new();

    // Track which seeds have connections
    let mut connected_seeds: std::collections::HashSet<NodeIndex> =
        std::collections::HashSet::new();
    for (fwd, rev) in &eval.connections {
        connected_seeds.insert(*fwd);
        connected_seeds.insert(*rev);
    }

    // Build a lookup from NodeIndex to metrics
    let metrics_map: HashMap<NodeIndex, (bool, &super::seed_eval::SeedMetrics)> = eval
        .seed_metrics
        .iter()
        .map(|(node, is_start, metrics)| (*node, (*is_start, metrics)))
        .collect();

    for members in component_map.values() {
        let mut start_seeds = Vec::new();
        let mut end_seeds = Vec::new();
        let mut has_connected = false;
        let mut total_branching_ratio = 0.0;
        let mut total_extension_nodes = 0usize;
        let mut max_mean_edge_count = 0.0f64;
        let mut min_coverage_cv = f64::INFINITY;
        let mut seed_count = 0usize;

        for &node in members {
            if connected_seeds.contains(&node) {
                has_connected = true;
            }

            if let Some(&(is_start, metrics)) = metrics_map.get(&node) {
                if is_start {
                    start_seeds.push(node);
                } else {
                    end_seeds.push(node);
                }
                let ratio = if metrics.node_count > 0 {
                    metrics.branching_events as f64 / metrics.node_count as f64
                } else {
                    0.0
                };
                total_branching_ratio += ratio;
                total_extension_nodes += metrics.node_count;
                if metrics.mean_edge_count > max_mean_edge_count {
                    max_mean_edge_count = metrics.mean_edge_count;
                }
                if metrics.coverage_cv < min_coverage_cv && metrics.node_count > 0 {
                    min_coverage_cv = metrics.coverage_cv;
                }
                seed_count += 1;
            } else {
                // Node in connections but not in seed_metrics (shouldn't happen,
                // but handle gracefully by checking graph)
                if graph[node].is_start {
                    start_seeds.push(node);
                }
                if graph[node].is_end {
                    end_seeds.push(node);
                }
            }
        }

        let avg_branching_ratio = if seed_count > 0 {
            total_branching_ratio / seed_count as f64
        } else {
            0.0
        };
        if !min_coverage_cv.is_finite() {
            min_coverage_cv = 0.0;
        }

        components.push(SeedComponent {
            start_seeds,
            end_seeds,
            has_connected_seeds: has_connected,
            avg_branching_ratio,
            total_extension_nodes,
            max_mean_edge_count,
            min_coverage_cv,
            node_budget: 0,      // set by allocate_budgets
            priority_score: 0.0, // set by prioritize_components
        });
    }

    components
}

/// Compute a priority score for a component. Higher = extend first.
///
/// The score is a weighted sum of factors. The connected_bonus dominates,
/// ensuring connected components always come first. Other factors break
/// ties and rank non-connected components meaningfully.
fn compute_priority_score(component: &SeedComponent) -> f64 {
    let connected_bonus = if component.has_connected_seeds {
        1000.0
    } else {
        0.0
    };

    // More seeds = more primer variants found = likely real target
    let seed_count_score = component.seed_count() as f64;

    // Lower branching = cleaner extension = more likely on-target
    let branching_penalty = -10.0 * component.avg_branching_ratio;

    // Larger total extension during seed eval = more promising
    let extension_score = (component.total_extension_nodes as f64).sqrt();

    // Higher mean edge coverage = path is in real amplicon territory, not
    // a low-coverage off-target region. Log-scale to compress the range
    // (coverage can span 2–2000+). ln(max_mean) gives ~0 at count 1,
    // ~7 at count 1000.
    let coverage_score = if component.max_mean_edge_count > 1.0 {
        component.max_mean_edge_count.ln() * 5.0
    } else {
        0.0
    };

    // Lower coefficient of variation = uniform coverage = real amplicon.
    // High CV means the path wanders between coverage regimes (e.g., into
    // a repeat). Penalty scales linearly; typical good paths have CV < 0.5,
    // wandering paths often have CV > 1.0.
    let uniformity_penalty = -10.0 * component.min_coverage_cv;

    connected_bonus
        + seed_count_score
        + branching_penalty
        + extension_score
        + coverage_score
        + uniformity_penalty
}

/// Sort components by priority score (highest first) and compute scores.
pub(super) fn prioritize_components(components: &mut [SeedComponent]) {
    for comp in components.iter_mut() {
        comp.priority_score = compute_priority_score(comp);
    }
    components.sort_by(|a, b| {
        b.priority_score
            .partial_cmp(&a.priority_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
}

/// Allocate node budgets to components proportional to seed count, with a
/// minimum floor per component.
///
/// If applying the requested floor to every component would cause their
/// summed budgets to exceed `total_budget`, the floor is scaled down to the
/// largest value that still fits (i.e. `total_budget / num_components`).
/// This preserves fairness across components while guaranteeing that the
/// sum of per-component budgets is `<= total_budget`.
pub(super) fn allocate_budgets(
    components: &mut [SeedComponent],
    total_budget: usize,
    min_per_component: usize,
) {
    if components.is_empty() {
        return;
    }

    let num_components = components.len();
    let effective_min = if min_per_component.saturating_mul(num_components) > total_budget {
        total_budget / num_components
    } else {
        min_per_component
    };

    let total_seeds: usize = components.iter().map(|c| c.seed_count().max(1)).sum();

    for comp in components.iter_mut() {
        let proportional =
            (comp.seed_count().max(1) as f64 / total_seeds as f64 * total_budget as f64) as usize;
        comp.node_budget = proportional.max(effective_min).min(total_budget);
    }

    // Final safety: if floating-point rounding pushed the sum above budget
    // (e.g. from several proportional allocations each rounding up), trim
    // from the lowest-priority components until the sum fits. Components
    // are ordered highest-priority first by `rank_components`.
    let sum: usize = components.iter().map(|c| c.node_budget).sum();
    if sum > total_budget {
        let mut excess = sum - total_budget;
        for comp in components.iter_mut().rev() {
            if excess == 0 {
                break;
            }
            let trim = excess.min(comp.node_budget.saturating_sub(1));
            comp.node_budget -= trim;
            excess -= trim;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pcr::seed_eval::{SeedEvalResult, SeedMetrics};

    fn make_metrics(node_count: usize, branching: usize, reached: bool) -> SeedMetrics {
        SeedMetrics {
            node_count,
            branching_events: branching,
            budget_exhausted: false,
            reached_opposite: reached,
            terminated: false,
            mean_edge_count: 10.0,
            coverage_cv: 0.3,
        }
    }

    #[test]
    fn test_single_connected_component() {
        let fwd = NodeIndex::new(0);
        let rev = NodeIndex::new(1);

        let eval = SeedEvalResult {
            seed_metrics: vec![
                (fwd, true, make_metrics(100, 2, true)),
                (rev, false, make_metrics(80, 1, true)),
            ],
            connections: vec![(fwd, rev)],
        };

        // Minimal graph just for the node data
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        graph.add_node(DBNode {
            sub_kmer: 0,
            is_start: true,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        graph.add_node(DBNode {
            sub_kmer: 1,
            is_start: false,
            is_end: true,
            is_terminal: false,
            visited: false,
        });

        let components = identify_components(&eval, &graph);
        assert_eq!(components.len(), 1);
        assert!(components[0].has_connected_seeds);
        assert_eq!(components[0].start_seeds.len(), 1);
        assert_eq!(components[0].end_seeds.len(), 1);
    }

    #[test]
    fn test_disconnected_singletons() {
        let fwd = NodeIndex::new(0);
        let rev = NodeIndex::new(1);

        let eval = SeedEvalResult {
            seed_metrics: vec![
                (fwd, true, make_metrics(50, 0, false)),
                (rev, false, make_metrics(30, 0, false)),
            ],
            connections: vec![], // no connections
        };

        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        graph.add_node(DBNode {
            sub_kmer: 0,
            is_start: true,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        graph.add_node(DBNode {
            sub_kmer: 1,
            is_start: false,
            is_end: true,
            is_terminal: false,
            visited: false,
        });

        let components = identify_components(&eval, &graph);
        assert_eq!(components.len(), 2);
        assert!(!components[0].has_connected_seeds);
        assert!(!components[1].has_connected_seeds);
    }

    #[test]
    fn test_connected_components_prioritized_first() {
        let fwd1 = NodeIndex::new(0);
        let rev1 = NodeIndex::new(1);
        let fwd2 = NodeIndex::new(2);

        let eval = SeedEvalResult {
            seed_metrics: vec![
                (fwd1, true, make_metrics(100, 2, true)),
                (rev1, false, make_metrics(80, 1, true)),
                (fwd2, true, make_metrics(200, 0, false)), // bigger but not connected
            ],
            connections: vec![(fwd1, rev1)],
        };

        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        for i in 0..3 {
            graph.add_node(DBNode {
                sub_kmer: i,
                is_start: i != 1,
                is_end: i == 1,
                is_terminal: false,
                visited: false,
            });
        }

        let mut components = identify_components(&eval, &graph);
        prioritize_components(&mut components);

        // Connected component should be first despite fwd2 having more nodes
        assert!(components[0].has_connected_seeds);
        assert!(!components[1].has_connected_seeds);
        assert!(components[0].priority_score > components[1].priority_score);
    }

    #[test]
    fn test_budget_allocation() {
        let fwd = NodeIndex::new(0);
        let rev = NodeIndex::new(1);
        let fwd2 = NodeIndex::new(2);

        let eval = SeedEvalResult {
            seed_metrics: vec![
                (fwd, true, make_metrics(100, 0, false)),
                (rev, false, make_metrics(50, 0, false)),
                (fwd2, true, make_metrics(50, 0, false)),
            ],
            connections: vec![(fwd, rev)],
        };

        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        for i in 0..3 {
            graph.add_node(DBNode {
                sub_kmer: i,
                is_start: i != 1,
                is_end: i == 1,
                is_terminal: false,
                visited: false,
            });
        }

        let mut components = identify_components(&eval, &graph);
        allocate_budgets(&mut components, 50000, 1000);

        // All components should have at least min_per_component
        for comp in &components {
            assert!(comp.node_budget >= 1000);
        }
    }

    /// When the requested min_per_component × num_components exceeds
    /// total_budget, the summed budget must not overflow total_budget.
    #[test]
    fn test_budget_allocation_respects_total() {
        let mut eval = SeedEvalResult {
            seed_metrics: Vec::new(),
            connections: Vec::new(),
        };
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        // 20 disconnected forward-only seed components.
        for i in 0..20u64 {
            let n = graph.add_node(DBNode {
                sub_kmer: i,
                is_start: true,
                is_end: false,
                is_terminal: false,
                visited: false,
            });
            eval.seed_metrics
                .push((n, true, make_metrics(10, 0, false)));
        }

        let mut components = identify_components(&eval, &graph);
        // total_budget=10000, min_per_component=1000, 20 components:
        // naive floor would give 20*1000 = 20000 > 10000.
        allocate_budgets(&mut components, 10000, 1000);

        let sum: usize = components.iter().map(|c| c.node_budget).sum();
        assert!(
            sum <= 10000,
            "sum of per-component budgets ({}) exceeds total_budget (10000)",
            sum
        );
    }

    #[test]
    fn test_empty_eval() {
        let eval = SeedEvalResult {
            seed_metrics: vec![],
            connections: vec![],
        };
        let graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let components = identify_components(&eval, &graph);
        assert!(components.is_empty());
    }

    // --- LiveComponentTracker tests ---

    fn make_seed_graph_for_tracker() -> (
        StableDiGraph<DBNode, DBEdge>,
        crate::kmer::KmerCounts,
        crate::kmer::KmerCounts,
    ) {
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let k = 5usize;
        let mut fwd_kmers = crate::kmer::KmerCounts::new(&k);
        let mut rev_kmers = crate::kmer::KmerCounts::new(&k);

        // Forward seed: sub_kmer = 0b00011011 = 27 (kmer ACGTX >> 2 = ACGT)
        let fwd_kmer: u64 = 0b0001101100; // ACGTA, k=5
        fwd_kmers.insert(&fwd_kmer, &10);
        graph.add_node(DBNode {
            sub_kmer: fwd_kmer >> 2, // ACGT = 27
            is_start: true,
            is_end: false,
            is_terminal: false,
            visited: false,
        });

        // Reverse seed: sub_kmer derived from revcomp
        // Use a different kmer so seeds are separate nodes
        let rev_kmer: u64 = 0b1100110011; // TGCTG, k=5
        rev_kmers.insert(&rev_kmer, &5);
        let rc = crate::kmer::revcomp_kmer(&rev_kmer, &k);
        let suffix_mask = (1u64 << (2 * (k - 1))) - 1;
        graph.add_node(DBNode {
            sub_kmer: rc & suffix_mask,
            is_start: false,
            is_end: true,
            is_terminal: false,
            visited: false,
        });

        (graph, fwd_kmers, rev_kmers)
    }

    #[test]
    fn test_tracker_init_two_seeds() {
        let (graph, fwd, rev) = make_seed_graph_for_tracker();
        let tracker = LiveComponentTracker::from_seed_graph(&graph, &fwd, &rev);

        // Two singleton components
        assert_eq!(tracker.info.len(), 2);
        assert_eq!(tracker.active_component_count(), 2);

        // Neither should be bridged
        for (_, ci) in tracker.components() {
            assert!(!ci.is_bridged);
            assert!(!ci.is_abandoned);
        }
    }

    #[test]
    fn test_tracker_register_and_merge() {
        let (graph, fwd, rev) = make_seed_graph_for_tracker();
        let mut tracker = LiveComponentTracker::from_seed_graph(&graph, &fwd, &rev);

        let fwd_node = NodeIndex::new(0);
        let rev_node = NodeIndex::new(1);

        // Register a new node belonging to the forward component
        let new_node = NodeIndex::new(2);
        tracker.uf.make_set(new_node); // need to make set before register
        tracker.register_node(new_node, fwd_node);

        // Forward component should have 1 node used
        let fwd_root = tracker.component_of(fwd_node);
        assert_eq!(tracker.info[&fwd_root].nodes_used, 1);

        // Now merge: new_node connects to rev_node's component
        let newly_bridged = tracker.merge(new_node, rev_node);

        // Should be bridged now (has both start and end seeds)
        assert!(newly_bridged);

        // All three nodes should be in the same component
        let root_a = tracker.component_of(fwd_node);
        let root_b = tracker.component_of(rev_node);
        let root_c = tracker.component_of(new_node);
        assert_eq!(root_a, root_b);
        assert_eq!(root_b, root_c);

        // The merged component should have both seeds
        let ci = &tracker.info[&root_a];
        assert!(!ci.start_seeds.is_empty());
        assert!(!ci.end_seeds.is_empty());
        assert!(ci.is_bridged);
    }

    #[test]
    fn test_tracker_merge_same_component_noop() {
        let (graph, fwd, rev) = make_seed_graph_for_tracker();
        let mut tracker = LiveComponentTracker::from_seed_graph(&graph, &fwd, &rev);

        let fwd_node = NodeIndex::new(0);

        // Merging a node with itself should be a no-op
        let bridged = tracker.merge(fwd_node, fwd_node);
        assert!(!bridged);
        assert_eq!(tracker.info.len(), 2); // still two components
    }

    #[test]
    fn test_tracker_abandon() {
        let (graph, fwd, rev) = make_seed_graph_for_tracker();
        let mut tracker = LiveComponentTracker::from_seed_graph(&graph, &fwd, &rev);

        let fwd_node = NodeIndex::new(0);
        tracker.abandon(fwd_node);

        let root = tracker.component_of(fwd_node);
        assert!(tracker.info[&root].is_abandoned);
        assert!(!tracker.info[&root].has_frontier());
        assert_eq!(tracker.active_component_count(), 1);
    }

    #[test]
    fn test_tracker_priority_ordering() {
        let (graph, fwd, rev) = make_seed_graph_for_tracker();
        let tracker = LiveComponentTracker::from_seed_graph(&graph, &fwd, &rev);

        let roots = tracker.roots_by_priority();
        assert_eq!(roots.len(), 2);

        // Forward seed had count=10, reverse had count=5
        // So forward component should be first (higher priority)
        let first_priority = tracker.info[&roots[0]].priority;
        let second_priority = tracker.info[&roots[1]].priority;
        assert!(
            first_priority >= second_priority,
            "Higher priority should come first: {} >= {}",
            first_priority,
            second_priority
        );
    }

    #[test]
    fn test_tracker_dual_role_seed_bridged() {
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let k = 5usize;
        let mut fwd_kmers = crate::kmer::KmerCounts::new(&k);
        let rev_kmers = crate::kmer::KmerCounts::new(&k);

        // A node that is both start and end (overlapping primers)
        let kmer: u64 = 0b0001101100; // ACGTA
        fwd_kmers.insert(&kmer, &10);
        let sub_kmer = kmer >> 2;
        // Find a reverse kmer whose revcomp suffix equals sub_kmer
        // For simplicity, just insert directly and add the node as both
        graph.add_node(DBNode {
            sub_kmer,
            is_start: true,
            is_end: true,
            is_terminal: false,
            visited: false,
        });
        // We need rev_kmers to have something for the is_end lookup
        // but the sub_kmer matching is complex; just test that the
        // is_start && is_end node is detected as bridged
        let tracker = LiveComponentTracker::from_seed_graph(&graph, &fwd_kmers, &rev_kmers);

        assert_eq!(tracker.info.len(), 1);
        for (_, ci) in tracker.components() {
            assert!(ci.is_bridged, "Dual-role seed should be bridged");
        }
    }
}
