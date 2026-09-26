#!/usr/bin/env python3
"""Post-outcome independent DP and ledger interpretation; never edits the assay."""
from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import time

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[2]


def sha(value):
    return hashlib.sha256(value.encode() if isinstance(value, str) else value).hexdigest()


def rc(value):
    return value.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@lru_cache(maxsize=None)
def independent_dp(reference, observed):
    """Unbanded rectangular DP with independently counted optimal paths."""
    scores = list(range(len(observed) + 1))
    ways = [1] * (len(observed) + 1)
    for i, reference_base in enumerate(reference, 1):
        next_scores, next_ways = [i], [1]
        for j, observed_base in enumerate(observed, 1):
            alternatives = ((scores[j-1] + (reference_base != observed_base), ways[j-1]),
                            (scores[j] + 1, ways[j]), (next_scores[j-1] + 1, next_ways[j-1]))
            best = min(item[0] for item in alternatives)
            next_scores.append(best)
            next_ways.append(min(2, sum(item[1] for item in alternatives if item[0] == best)))
        scores, ways = next_scores, next_ways
    return scores[-1], ways[-1] > 1


@lru_cache(maxsize=None)
def minimum_substring_distance(query, target):
    """Global query, free target ends; <=4 implies substring length within4."""
    previous = [0] * (len(target) + 1)
    for i, base in enumerate(query, 1):
        row = [i]
        for j, other in enumerate(target, 1):
            row.append(min(previous[j-1] + (base != other), previous[j] + 1, row[j-1] + 1))
        previous = row
    return min(previous)


def validate_script(reference, observed, match):
    i = j = 0
    counts = Counter()
    pieces = re.findall(r"(\d+)([=XID])", match["cigar"])
    assert "".join(count + op for count, op in pieces) == match["cigar"]
    for count, operation in pieces:
        for _ in range(int(count)):
            counts[operation] += 1
            if operation in ("=", "X"):
                assert (reference[i] == observed[j]) == (operation == "=")
                i += 1
                j += 1
            elif operation == "I":
                j += 1
            else:
                i += 1
    assert (i, j) == (len(reference), len(observed))
    assert [counts[op] for op in "XID"] == [match[key] for key in ("substitutions", "insertions", "deletions")]
    assert sum(counts[op] for op in "XID") == match["distance"]


def identity(match):
    return match["class_id"], match["strand"], tuple(match["read_interval"])


def main():
    started = time.monotonic()
    freeze = json.loads((BASE / "freeze.json").read_text())
    output = Path(freeze["run_dir"]) / "recruitment.json"
    raw = output.read_bytes()
    result = json.loads(raw)
    manifest = json.loads((BASE / "inputs/manifest.json").read_text())
    regions = {region["id"]: region for region in manifest["regions"]}
    decoys = {item["id"]: item["sequence"] for item in manifest["decoys"]}
    counts = Counter()
    rows = {key: [] for key in regions}
    best_details = []
    for read in result["read_evidence"]:
        sequence, quality = read["sequence"], read["quality"]
        oriented = {"+": (sequence, quality), "-": (rc(sequence), quality[::-1])}
        for observation in read["observations"]:
            region = regions[observation["region_id"]]
            classes = {item["id"]: item for item in region["classes"]}
            rows[region["id"]].append((read, observation))
            counts["observations"] += 1
            for anchor in observation["anchors"]:
                seq, qual = oriented[anchor["strand"]]
                a, b = anchor["start"], anchor["end"]
                assert anchor["q20"] == (min(map(ord, qual[a:b])) >= 53)
                assert anchor["read_interval"] == ([a, b] if anchor["strand"] == "+" else [len(seq)-b, len(seq)-a])
                legal_patterns = {value["sequence"][:21] if anchor["side"] == "left" else value["sequence"][-21:] for value in classes.values()}
                for match in anchor["flank_matches"]:
                    assert match["sequence"] in legal_patterns
                    distance, _ = independent_dp(match["sequence"], seq[a:b])
                    assert distance == match["distance"] and distance <= 2
                    counts["flank_matches_checked"] += 1
            expected_pairs = Counter((left["strand"], left["start"], left["end"], right["start"], right["end"])
                for left in observation["anchors"] if left["side"] == "left"
                for right in observation["anchors"] if right["side"] == "right" and right["strand"] == left["strand"]
                and left["end"] <= right["start"])
            actual_pairs = Counter((pair["strand"], *pair["left_interval"], *pair["right_interval"])
                                   for pair in observation["pairs"])
            assert expected_pairs == actual_pairs
            alternatives, qualified = {}, {}
            for pair in observation["pairs"]:
                seq, qual = oriented[pair["strand"]]
                a, b = pair["start"], pair["end"]
                span = seq[a:b]
                assert pair["observed_span"] == span and pair["span_length"] == b-a
                assert pair["q20"] == (min(map(ord, qual[a:b])) >= 53)
                assert pair["read_interval"] == ([a, b] if pair["strand"] == "+" else [len(seq)-b, len(seq)-a])
                expected_matches = {}
                if set(span) <= set("ACGT"):
                    for class_id, item in classes.items():
                        if abs(len(item["sequence"])-len(span)) > 4:
                            continue
                        distance, tied = independent_dp(item["sequence"], span)
                        counts["class_pair_distances_checked"] += 1
                        if distance <= 4:
                            expected_matches[class_id] = (distance, tied)
                assert {match["class_id"] for match in pair["class_matches"]} == set(expected_matches)
                for match in pair["class_matches"]:
                    assert (match["distance"], match["tied_optimal_alignments"]) == expected_matches[match["class_id"]]
                    validate_script(classes[match["class_id"]]["sequence"], span, match)
                    key = match["class_id"], pair["strand"], tuple(pair["read_interval"])
                    alternatives[key] = match["distance"]
                    if pair["q20"]:
                        qualified[key] = match["distance"]
                    counts["cigars_checked"] += 1
                counts["pairs_checked"] += 1
            for key, candidate_map in (("best_unfiltered", alternatives), ("best_q20", qualified)):
                expected = {item for item, distance in candidate_map.items() if distance == min(candidate_map.values())} if candidate_map else set()
                assert {identity(best) for best in observation[key]} == expected
                for best in observation[key]:
                    seq, qual = oriented[best["strand"]]
                    a, b = best["oriented_interval"]
                    span = seq[a:b]
                    item = classes[best["class_id"]]
                    assert independent_dp(item["sequence"], span) == (best["distance"], best["tied_optimal_alignments"])
                    observed_context_minima = {}
                    for decoy_id, decoy in decoys.items():
                        for strand, target in (("+", decoy), ("-", rc(decoy))):
                            minimum = minimum_substring_distance(span, target)
                            hits = [hit for hit in best["decoy_placements"] if hit["decoy_id"] == decoy_id and hit["strand"] == strand]
                            assert bool(hits) == (minimum <= 4)
                            if hits:
                                assert min(hit["distance"] for hit in hits) == minimum
                            observed_context_minima[decoy_id + ":" + strand] = minimum
                            counts["context_minima_checked"] += 1
                    for hit in best["decoy_placements"]:
                        left, right = hit["interval"]
                        target = decoys[hit["decoy_id"]][left:right]
                        if hit["strand"] == "-":
                            target = rc(target)
                        assert independent_dp(span, target)[0] == hit["distance"]
                        counts["decoy_placements_checked"] += 1
                    best_details.append({"ordinal": read["ordinal"], "read_id": read["read_id"],
                        "region_id": region["id"], "qualification": key, "class_id": best["class_id"],
                        "full_read_sequence_sha256": sha(sequence), "span_sha256": sha(span),
                        "span_length": len(span), "strand": best["strand"], "oriented_interval": [a, b],
                        "physical_interval": best["read_interval"], "distance": best["distance"],
                        "cigar": best["cigar"], "tied_optimal_alignments": best["tied_optimal_alignments"],
                        "minimum_observed_quality": min(map(ord, qual[a:b]))-33,
                        "low_quality_span_offsets": [i for i, value in enumerate(qual[a:b]) if ord(value)<53],
                        "finite_decoy_placements": len(best["decoy_placements"]),
                        "minimum_finite_decoy_distance": min(observed_context_minima.values())})
            counts["best_sets_checked"] += 2
    summary_checks = []
    for summary in result["regions"]:
        selected = rows[summary["region_id"]]
        recalc = {"recruited_records": len(selected),
                  "distinct_recruited_read_ids": len({read["read_id"] for read, _ in selected})}
        for side in ("left", "right"):
            recalc[side+"_records"] = sum(any(a["side"]==side for a in obs["anchors"]) for _, obs in selected)
            recalc[side+"_only_records"] = sum({a["side"] for a in obs["anchors"]}=={side} for _, obs in selected)
        recalc["both_flanks_records"] = sum(len({a["side"] for a in obs["anchors"]})==2 for _, obs in selected)
        recalc["ordered_pair_records"] = sum(bool(obs["pairs"]) for _, obs in selected)
        recalc["no_ordered_pair_records"] = sum(not obs["pairs"] for _, obs in selected)
        recalc["unique_records"] = sum(len(obs["best_q20"])==1 for _, obs in selected)
        recalc["ambiguous_records"] = sum(len(obs["best_q20"])>1 for _, obs in selected)
        for key, value in recalc.items():
            assert summary[key] == value, (summary["region_id"], key)
        summary_checks.append({"region_id": summary["region_id"], **recalc})
    its2 = [read for read in result["read_evidence"] if any(obs["region_id"]=="region_2" and obs["best_q20"] for obs in read["observations"])]
    assert len(its2)==2 and its2[0]["sequence"] == its2[1]["sequence"]
    payload = {"status": "pass", "scope": "Independent post-outcome verification of reported ledger alignments; not a new biological assay",
               "freeze_sha256": sha((BASE/"freeze.json").read_bytes()), "recruitment_output_sha256": sha(raw),
               "checks": dict(counts), "distinct_dp_queries": independent_dp.cache_info().currsize,
               "distinct_context_dp_queries": minimum_substring_distance.cache_info().currsize,
               "recruited_records": len(result["read_evidence"]), "region_counts_recomputed": summary_checks,
               "its_best_details": [row for row in best_details if row["region_id"] != "co1_control"],
               "its2_identical_full_read_ordinals": [read["ordinal"] for read in its2],
               "its2_full_read_sequence_sha256": sha(its2[0]["sequence"]),
               "qualified_co1_best_distance_counts": dict(Counter(str(row["distance"]) for row in best_details if row["region_id"]=="co1_control" and row["qualification"]=="best_q20")),
               "elapsed_seconds": time.monotonic()-started}
    destination = BASE/"results_review/local_alignment_review.json"
    with destination.open("x") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": payload["status"], "checks": payload["checks"], "elapsed_seconds":payload["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
