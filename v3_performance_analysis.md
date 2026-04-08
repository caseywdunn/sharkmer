# v3.0 Performance Analysis: Drosophila melanogaster (insecta panel)

Direct head-to-head comparison of v2.0.0 vs v3.0.0-rc on the same machine,
same input file (SRR31887760, 1M reads cached locally), same hardware
(Linux aarch64, 47GB RAM).

## Runtime comparison

### Varying k at default budgets (v2=50K fixed, v3=100K dynamic)

| k  | v2 (50K) | v3 (100K) | v3/v2  |
|----|----------|-----------|--------|
| 31 | 38s      | 39s       | 1.0x   |
| 21 | 47s      | 13m 21s   | 17x    |
| 19 | 47s      | 19m 39s   | 25x    |

### Varying k at equal budget (both 50K)

| k  | v2 (50K) | v3 (50K) | v3/v2  |
|----|----------|----------|--------|
| 31 | 38s      | 34s      | 0.9x   |
| 21 | 47s      | 2m 32s   | 3.2x   |
| 19 | 47s      | 1m 29s   | 1.9x   |

**Key finding:** At k=31 with equal budget, v3 is the same speed as v2. The
slowdown is entirely driven by the interaction between lower k and v3's
per-component architecture.

## Gene recovery comparison

| Gene    | v2 k=31 | v2 k=21 | v2 k=19 | v3 k=31 | v3 k=21 | v3 k=19 (100K) |
|---------|---------|---------|---------|---------|---------|-----------------|
| 12S     | YES     | YES     | -       | -       | -       | -               |
| 18S     | YES     | YES     | YES     | YES     | YES     | YES             |
| 18S-v2  | YES     | YES     | YES     | YES     | YES     | YES             |
| 28S     | -       | -       | -       | -       | -       | **YES**         |
| CO1-v2  | YES     | YES     | YES     | YES     | YES     | YES             |
| CO2     | YES     | YES     | YES     | YES     | YES     | YES             |
| CO2-v2  | YES     | YES     | YES     | YES     | YES     | YES             |
| ND1     | YES     | YES     | YES     | -       | -       | -               |
| ND5     | YES     | -       | -       | -       | -       | -               |
| Yp2     | -       | -       | -       | -       | -       | **YES**         |
| **Total** | **8** | **7**   | **6**   | **5**   | **5**   | **7**           |

v3 at k=19/100K recovers 7 genes — same count as v2 at k=21, but different
genes. v3 gains 28S and Yp2 (via reverse extension and larger budget) but
loses 12S, ND1, ND5 (node budget exceeded at all k values — a graph
traversal behavioral change, not a budget or k issue).

## Root cause analysis

### Why v3 is slow at lower k

Every successful gene was found at **component 1** (the highest-priority
component). One gene (Yp2, k=19 only) was found at component 4. No
successful gene ever needed more than 4 component attempts.

Failed genes exhaust **every component** before giving up:

| Gene    | Components | Connected | Tried | Outcome              |
|---------|------------|-----------|-------|----------------------|
| ND4     | 377        | 0         | 377   | node budget exceeded |
| ND1     | 339        | 1         | 339   | node budget exceeded |
| CO2*    | 267        | 2         | 267   | node budget exceeded |
| CytB    | 247        | 1         | 247   | node budget exceeded |
| ND5     | 246        | 1         | 246   | node budget exceeded |
| CO1     | 244        | 1         | 244   | node budget exceeded |
| 12S     | 190        | 0         | 190   | node budget exceeded |
| ...     | ...        | ...       | ...   | ...                  |

*CO2 has 267 components total but succeeds at component 1 in first-product
mode, so it only tries 1. The 267 figure is from an all-components run.

Across all 16 failed genes at k=21/50K: **1,873 total component attempts**
wasted. Each attempt involves forward extension, reachability pruning, DFS
path finding, reverse extension, pruning, and DFS again — 6 graph operations
per component. Individual operations are fast (~0ms each) but 1,873 x 6 =
~11,000 operations adds up.

### Why so many components?

The component count is driven by the number of surviving seed nodes after
seed evaluation. The number of seed nodes is driven by the number of primer
kmers entering `create_seed_graph`. This is supposed to be capped at 100 per
direction by `--max-primer-kmers`, but the cap has a tie-breaking flaw:

```rust
if counts.len() > max_primer_kmers {
    top_count_cutoff = &counts[max_primer_kmers - 1];
}
// ... keeps all kmers with count >= top_count_cutoff
```

When many primer kmers share the same count (common with degenerate primers
at low coverage where most matches sit at `min_count=2`), the cutoff equals
that shared count and **all** kmers pass. ND4 has 332 reverse primer kmers
instead of the intended 100 cap.

### Why this didn't matter in v2

v2 has no seed evaluation and no components. All primer kmers seed a single
graph, extension runs once until the 50K node cap, and the gene is done.
Whether there are 100 or 332 seed nodes is irrelevant — they're all starting
points in the same walk. The leaky cap was a latent bug that v3's
per-component architecture turned into a performance problem.

## Architecture comparison

### v2 per-gene flow (one pass)

```
Seed all primer kmers → Extend single graph → Hit budget → Done
```

Cost: O(budget) regardless of primer kmer count.

### v3 per-gene flow (per-seed)

```
Seed all primer kmers into graph
→ Evaluate each seed individually (500 nodes, bounded exploration)
→ Surviving seeds become "components" via union-find
  (seeds are merged only if a forward seed reached a reverse seed
   during bounded eval, proving they're part of the same amplicon —
   but with off-target matches this almost never happens, so nearly
   every seed becomes a singleton component)
→ For each component (almost always = one surviving seed):
    → Forward extend (budget/N nodes)
    → Reachability prune (2x BFS)
    → DFS path search
    → Reverse extend (budget/N nodes)
    → Reachability prune (2x BFS)
    → DFS path search
→ Stop at first success or exhaust all components
```

The "component" system is effectively a prioritized list of individual
seeds, each getting its own full extension cycle. For ND4: 377 components,
0 connected — every one is a singleton seed getting its own extension.

Cost: O(N × budget/N) = O(budget) in total nodes explored, but with O(N)
overhead from pruning and path-finding per seed. When N=377 (ND4), this
overhead dominates even though each individual operation is fast.

## Implications for v3.0 release

### Fix 1: Enforce hard cap on primer kmers (bug fix)

The `filter_primer_kmers` function should enforce a strict count limit rather
than a count-value threshold. When ties exceed the cap, break ties
deterministically (e.g., by kmer value) and keep exactly `max_primer_kmers`.

Expected impact: reduces component counts by roughly 2-3x for the worst
offenders (ND4 from 377 to ~150-200), proportional runtime improvement.

### Fix 2: Cap disconnected component attempts (new parameter)

For `first-product` mode, stop trying disconnected singleton components
after a small number of failures (e.g., 5). Always try all connected
components (where forward seeds reached reverse seeds during evaluation),
since these have strong signal. A cap of 5 would retain every successful
gene recovery observed in testing while eliminating ~95% of wasted component
attempts.

Expected impact: at k=21/50K budget, reduces wasted attempts from ~1,873 to
~80 (16 failed genes x 5 max attempts). This alone would bring v3 runtime
close to v2 at the same k and budget.

### Combined effect

With both fixes, v3 at k=19 should run in under 2 minutes on this dataset
(vs current 20 minutes), while retaining all gene recovery improvements
(28S, Yp2). The residual overhead vs v2 would be the legitimate cost of
seed evaluation and reverse extension — the features that enable the new
gene recoveries.

### What the fixes don't address

12S, ND1, and ND5 are genuine regressions that occur at all k values and
all budget levels in v3. These are behavioral changes in v3's graph
traversal, not consequences of the performance problem. They likely need
investigation into why v3's component-based extension fails to find paths
that v2's single-graph extension finds — possibly the per-component budget
fragmentation (50000/377 = 132 nodes per component) starves the on-target
component of budget that v2 would give it in full.
