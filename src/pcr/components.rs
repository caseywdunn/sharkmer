// pcr/components.rs — Seed component identification, prioritization, and budget allocation

use std::collections::HashMap;

use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;

use super::seed_eval::SeedEvalResult;
use super::{DBEdge, DBNode};

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
        let p = self.parent[&x];
        if p != x {
            let root = self.find(p);
            self.parent.insert(x, root);
            root
        } else {
            x
        }
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
        comp.node_budget = proportional.max(effective_min);
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
}
