#!/usr/bin/env python3
"""Independently recount completed comparison outputs; never read FASTQ or assemble.

Reuses the reviewed transaction/metadata validator for its published schema.
Parity, historical exact/RC matches and threshold retention conservation are
recomputed here from on-disk YAML/FASTA, not copied from the producer assessment.
No graph reconstruction, biological repeat-length validation or power claim.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys

import yaml

ROOT = Path(__file__).resolve().parents[4]
BUNDLE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(BUNDLE))
from sharkmer_validate import runner
from build_snapshot import fingerprint, write_new

ORDER = [
    ("SRR27962769", "baseline_off"), ("SRR27962769", "changed_off"), ("SRR27962769", "changed_on"),
    ("SRR31887760", "changed_on"), ("SRR31887760", "baseline_off"), ("SRR31887760", "changed_off"),
    ("SRR1057608", "changed_off"), ("SRR1057608", "changed_on"), ("SRR1057608", "baseline_off"),
]
CAPS = {
    "threshold": dict(paths=32, sequence_bases=524288, marker_occurrences=4096),
    "gene": dict(paths=64, sequence_bases=1048576, marker_occurrences=8192),
    "run": dict(paths=256, sequence_bases=4194304, marker_occurrences=32768),
}
COUNTS = ["n_reads_read", "n_bases_read", "n_subreads_ingested", "n_bases_ingested", "n_kmers"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify(row):
    require(fingerprint(row["path"]) == row, f"Identity changed: {row['path']}")


def rc(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_fasta(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            require(bool(line[1:]), "Empty FASTA identifier")
            rows.append({"header": line[1:], "sequence": ""})
        else:
            require(bool(rows) and bool(re.fullmatch("[ACGT]+", line)), "Invalid FASTA sequence line")
            rows[-1]["sequence"] += line
    require(all(row["sequence"] for row in rows), "Empty FASTA product")
    return rows


def historical_matches(sequence, rows):
    return {
        "exact_oriented": [r["header"] for r in rows if r["sequence"] == sequence],
        "reverse_complement_only": [r["header"] for r in rows if r["sequence"] != sequence and r["sequence"] == rc(sequence)],
    }


def retention(entry, index):
    value = entry["uncertain_products"]
    limits = value["limits"]
    for scope, expected in CAPS.items():
        require(limits[scope] == expected, f"Changed {scope} caps")
    share = {key: n // 21 + (index < n % 21) for key, n in CAPS["run"].items()}
    require(limits["allocated_run_share_for_gene"] == share and limits["unused_share_is_not_redistributed"] is True, "Changed full-panel budget share")
    require([{k: v for k, v in row.items() if k != "retention"} for row in value["thresholds"]] == entry["threshold_diagnostics"], "Uncertain threshold search counters differ")
    previous = dict(paths=0, sequence_bases=0, marker_occurrences=0)
    summaries = {}
    observed = dropped = node_skips = 0
    for threshold in value["thresholds"]:
        require(threshold["threshold"] not in summaries, "Duplicate threshold")
        kept = threshold.get("retention")
        summaries[threshold["threshold"]] = kept
        if kept is None:
            continue
        require(kept["gene_usage_before_threshold"] == previous, "Discontinuous gene retention accounting")
        increments = dict(paths=kept["retained_paths"], sequence_bases=kept["retained_sequence_bases"], marker_occurrences=kept["retained_marker_occurrences"])
        previous = {key: previous[key] + n for key, n in increments.items()}
        require(kept["gene_usage_after_threshold"] == previous, "Retained usage increment differs")
        require(all(n <= CAPS["threshold"][key] and previous[key] <= min(share[key], CAPS["gene"][key]) for key, n in increments.items()), "Retention exceeds cap")
        losses = sum(kept[k] for k in ["dropped_oversize_sequence", "dropped_oversize_marker_set", "dropped_path_cap", "dropped_sequence_base_cap", "dropped_marker_cap"])
        require(losses == kept["observed_not_retained_paths"], "Retention drop reasons do not conserve")
        require(kept["observed_withheld_paths"] == kept["retained_paths"] + losses == threshold["withheld_candidate_paths"], "Observed/retained/withheld path counts differ")
        require(kept["retention_truncated"] == (losses > 0), "Retention truncation flag differs")
        observed += kept["observed_withheld_paths"]
        dropped += losses
        node_skips += kept["node_visit_skips"]
    require(previous["paths"] == value["retained_path_records"], "Gene retained path total differs")
    require(value["path_search_evaluated"] == any(t["path_search_evaluated"] for t in value["thresholds"]), "Path-search evaluation flag differs")
    require(value["retention_truncated"] == (dropped > 0), "Gene retention flag differs")
    require(value["search_limits_reached"] == (node_skips > 0 or any(t["node_budget_reached"] or t["dfs_limit_reached"] or t["path_quota_reached"] for t in value["thresholds"])), "Gene search-limit flag differs")
    source_usage = {key: Counter() for key in previous}
    for product in value["products"]:
        for source in product["sources"]:
            threshold = source["threshold"]
            require(threshold in summaries and summaries[threshold] is not None, "Source references an uncollected threshold")
            require(source["node_visit_limit"] == summaries[threshold]["node_visit_limit"], "Source visit bound differs from threshold")
            source_usage["paths"][threshold] += 1
            source_usage["sequence_bases"][threshold] += product["sequence_length"]
            source_usage["marker_occurrences"][threshold] += len(source["marker_occurrences"])
    fields = {"paths": "retained_paths", "sequence_bases": "retained_sequence_bases", "marker_occurrences": "retained_marker_occurrences"}
    for threshold, kept in summaries.items():
        if kept:
            for key, field in fields.items():
                require(source_usage[key][threshold] <= kept[field], "Emitted sources exceed retained payload")
                if value["excluded_standard_equivalent_path_records"] == 0:
                    require(source_usage[key][threshold] == kept[field], "Emitted sources fail complete retained-payload accounting")
    return {"gene": entry["gene_name"], "standard_products": entry["n_products"], "uncertain_products": value["unique_products"], "retained_path_records": previous["paths"], "observed_withheld_paths": observed, "dropped_path_records": dropped, "excluded_standard_equivalent_path_records": value["excluded_standard_equivalent_path_records"], "exact_duplicate_path_records": value["exact_duplicate_path_records"], "allocated_share": share, "search_limits_reached": value["search_limits_reached"], "retention_truncated": value["retention_truncated"], "node_visit_skips": node_skips}


def audit(comparison, freeze_sha):
    source_before = fingerprint(__file__)
    freeze_path = comparison / "freeze.json"
    freeze_identity = fingerprint(freeze_path)
    require(freeze_identity["sha256"] == freeze_sha, "External freeze hash differs")
    frozen = json.loads(freeze_path.read_text())
    orchestration = json.loads((comparison / "orchestration.json").read_text())
    execution = json.loads((comparison / "runs/execution.json").read_text())
    reported = json.loads((comparison / "runs/assessment.json").read_text())
    require(orchestration["status"] == execution["status"] == "complete" and reported["status"] == "pass", "Comparison did not pass")
    require(orchestration["freeze"] == execution["freeze_before"] == execution["freeze_after"] == freeze_identity, "Freeze bookends differ")
    verify(orchestration["execution"])
    verify(execution["assessment"])
    require(frozen["limits"] == dict(address_space_bytes=40 * 1024**3, timeout_seconds=1800, threads=2, k=19, chunks=0, max_reads=1000000), "Limits differ")
    require([(j["accession"], j["role"]) for j in frozen["jobs"]] == ORDER, "Frozen order differs")
    require([j["id"] for j in execution["jobs"]] == [j["id"] for j in frozen["jobs"]], "Execution jobs differ")
    for row in frozen["sources"].values():
        verify(row)
    artifacts = {freeze_path, comparison / "orchestration.json", comparison / "runs/execution.json", comparison / "runs/assessment.json"}
    observed = {}
    costs = []
    summaries = {}
    for plan, item in zip(frozen["jobs"], execution["jobs"]):
        verify(item["receipt"])
        receipt_path = Path(item["receipt"]["path"])
        receipt = json.loads(receipt_path.read_text())
        artifacts.add(receipt_path)
        require(item["status"] == receipt["status"] == "complete" and receipt["returncode"] == 0 and not receipt["timed_out"], "Failed job")
        require(all(receipt[k] == v for k, v in plan.items()), "Receipt differs from planned job")
        expected = frozen["samples"][plan["accession"]]["expected_prefix"]
        require(receipt["prefix_before"] == receipt["prefix_after"] == expected, "Consumed-prefix bookends differ")
        directory = Path(plan["output"])
        require(set(receipt["outputs"]["files"]) == {p.name for p in directory.iterdir() if p.is_file()}, "Output file inventory differs")
        for row in receipt["outputs"]["files"].values():
            verify(row)
            artifacts.add(Path(row["path"]))
        constraints = json.loads((directory / "process_constraints.json").read_text())
        require(constraints["status"] == "pass" and constraints["command"] == plan["argv"] and constraints["affinity"] == frozen["cpus"] and constraints["rlimit_as_soft_hard"] == [40 * 1024**3] * 2, "Actual process constraints differ")
        name = "uncertain_" + plan["accession"]
        stats = yaml.safe_load((directory / (name + ".stats.yaml")).read_text())
        require(stats == receipt["outputs"]["stats"], "Raw YAML differs from captured stats")
        runner._validate_output_manifest(stats, name, directory)
        require(stats["n_reads_read"] == expected["records"] and stats["n_bases_read"] == expected["bases"], "Actual consumed stats differ")
        require(len(stats["pcr_results"]) == 21 and len({e["gene_name"] for e in stats["pcr_results"]}) == 21, "Wrong gene set size")
        require(stats["stage_timings"].get("read_threading_input_s") is None, "Threading evaluated")
        rows = {"standard": [], "uncertain": [], "stats": stats}
        gene_summaries = []
        for index, entry in enumerate(stats["pcr_results"]):
            require(("uncertain_products" in entry) == (plan["role"] == "changed_on"), "Wrong flag assessment presence")
            standard = read_fasta(directory / entry["output_file"]) if entry.get("output_file") else []
            require(len(standard) == entry["n_products"] and [len(r["sequence"]) for r in standard] == entry.get("product_lengths", []), "Standard count/length differs")
            require(entry["status"] != "fail" or not standard, "Failed gene published standard output")
            rows["standard"].extend(dict(gene=entry["gene_name"], **r) for r in standard)
            if plan["role"] == "changed_on":
                value = entry["uncertain_products"]
                uncertain = read_fasta(directory / value["output_file"]) if value.get("output_file") else []
                require(uncertain == [{"header": p["id"] + " " + p["description"], "sequence": p["oriented_sequence"]} for p in value["products"]], "Uncertain FASTA differs from metadata")
                ordinary = {s for r in standard for s in [r["sequence"], rc(r["sequence"]) ]}
                require(all(r["sequence"] not in ordinary for r in uncertain), "Uncertain product equals final standard exact/RC")
                rows["uncertain"].extend(dict(gene=entry["gene_name"], **r) for r in uncertain)
                gene_summaries.append(retention(entry, index))
        require(rows["standard"] == receipt["outputs"]["standard"] and rows["uncertain"] == receipt["outputs"]["uncertain"], "Raw FASTA differs from captured rows")
        observed[(plan["accession"], plan["role"])] = rows
        if gene_summaries:
            summaries[plan["accession"]] = gene_summaries
        time_text = (directory / "time.txt").read_text()
        rss = int(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", time_text)[1]) * 1024
        require(rss == receipt["maximum_rss_bytes"], "RSS differs from GNU time")
        costs.append(dict(accession=plan["accession"], role=plan["role"], wall_seconds=receipt["wall_seconds"], maximum_rss_bytes=rss))
    comparisons = []
    for accession in ["SRR27962769", "SRR31887760", "SRR1057608"]:
        for left, right in [("baseline_off", "changed_off"), ("changed_off", "changed_on")]:
            a, b = [observed[(accession, role)] for role in [left, right]]
            for key in COUNTS + ["n_multi_kmers", "n_singleton_kmers"]:
                require((key in a["stats"]) == (key in b["stats"]) and a["stats"].get(key) == b["stats"].get(key), f"Count parity differs: {accession}/{key}")
            require(a["standard"] == b["standard"], "Ordered standard FASTA parity differs")
            cleaned = [[{k: v for k, v in e.items() if k != "uncertain_products"} for e in s["stats"]["pcr_results"]] for s in [a, b]]
            require(cleaned[0] == cleaned[1], "Preexisting gene/threshold stats parity differs")
        comparisons.append(dict(accession=accession, standard_products=len(observed[(accession,"changed_on")]["standard"]), uncertain_products=len(observed[(accession,"changed_on")]["uncertain"]), standard_parity=True))
    losses = json.loads((BUNDLE / "inputs/historical_losses.json").read_text())["rows"]
    require(len(losses) == len(reported["historical_losses"]) == 7, "Wrong loss set")
    results = []
    for loss, producer in zip(losses, reported["historical_losses"]):
        sequence = loss["sequence"]
        require(len(sequence) == loss["length"] and hashlib.sha256(sequence.encode()).hexdigest() == loss["sequence_sha256"], "Historical identity differs")
        matches = {}
        for role in ["baseline_off", "changed_off", "changed_on"]:
            matches[role] = {kind: historical_matches(sequence, [r for r in observed[(loss["accession"], role)][kind] if r["gene"] == "insecta_" + loss["gene"]]) for kind in ["standard", "uncertain"]}
        require(matches == producer["matches"], "Producer historical classification differs")
        results.append(dict(accession=loss["accession"], gene=loss["gene"], length=loss["length"], sequence_sha256=loss["sequence_sha256"], matches=matches))
    for comparison_row in comparisons:
        producer = next(r for r in reported["comparisons"] if r["accession"] == comparison_row["accession"])
        require(producer["passed"] and not producer["failures"] and producer["uncertain_count"] == comparison_row["uncertain_products"] and set(producer["standard_counts"].values()) == {comparison_row["standard_products"]}, "Producer comparison summary differs")
    bound = {str(p.relative_to(ROOT)): fingerprint(p) for p in sorted(artifacts)}
    for row in bound.values():
        verify(row)
    for row in frozen["sources"].values():
        verify(row)
    require(fingerprint(__file__) == source_before, "Auditor changed")
    return dict(status="pass", source=source_before, freeze=freeze_identity, comparisons=comparisons, historical_losses=results, gene_retention=summaries, observed_costs=costs, reviewed_artifacts=bound, limitations=["Reads were not rescanned; exact prefix identities/counts come from independently verified preexecution byte hashes and frozen per-run before/after receipts.", "Transaction and metadata schema validation reuses the reviewed production validator; retention conservation and output parity/exact-match classifications were independently recalculated.", "No graph reconstruction or biological repeat-length/read-support validation. Non-observation is restricted by unchanged search and 21-target retention shares.", "Single sequential trio per sample; runtime and RSS are descriptive only."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--freeze-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Review output already exists")
    result = audit(args.comparison.resolve(), args.freeze_sha256)
    write_new(args.output, result)
    print(json.dumps(fingerprint(args.output)))


if __name__ == "__main__":
    main()
