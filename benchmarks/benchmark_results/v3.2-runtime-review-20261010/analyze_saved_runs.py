#!/usr/bin/env python3
"""Recompute CPU, stage, resource and identity evidence for four runtime flags."""

import csv
import hashlib
import io
import json
import statistics
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SAMPLES = {"SRR17535371": "Homo sapiens", "SRR25378184": "Liriodendron tulipifera",
           "SRR8617500": "Rhopilema esculentum", "SRR1057608": "Heliconius pachinus"}
PHASES = {
    "release_to_inclusive": {"directory": "tmp/inclusive-release-20260929/comparison-v1", "state": "execution.json",
        "binary_labels": {"baseline": "released_v3.1.0", "candidate": "inclusive_ed04bde"}},
    "inclusive_to_retention": {"directory": "tmp/sequence-retention-20261001/full-v3", "state": "result.json",
        "binary_labels": {"baseline": "inclusive_ed04bde", "candidate": "sequence_retention_checked_v3"}},
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def tsv(rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: json.dumps(v, separators=(",", ":"), sort_keys=True) if isinstance(v, (dict, list)) else v for k, v in row.items()})
    return output.getvalue().encode()


def main():
    sources, jobs, machines, sets, changed_context_bindings = {}, [], {}, {}, []

    def read(path, expected=None):
        path = Path(path)
        raw = path.read_bytes()
        seen = {"sha256": digest(raw), "size_bytes": len(raw)}
        if expected and (seen["sha256"] != expected["sha256"] or seen["size_bytes"] != expected.get("size_bytes", expected.get("bytes"))):
            raise ValueError("Saved source fingerprint differs: " + str(path))
        sources[str(path)] = seen
        return raw

    prior = ROOT / "benchmarks/benchmark_results/v3.2-inclusive-release-20260929/review/resource_assessment.json"
    old_review = json.loads(read(prior))
    old_bound = {str(ROOT / key): value for key, value in old_review["reviewed_artifacts"].items()}
    for phase, spec in PHASES.items():
        directory = ROOT / spec["directory"]
        state = json.loads(read(directory / spec["state"]))
        plan = json.loads(read(directory / "plan.json"))
        if state["status"] != "pass":
            raise ValueError("Saved comparison state failed")
        machines[phase] = {"machine": state["machine"], "slurm_job_id": state["slurm_job_id"],
                           "cpus": plan["cpus"], "topology": plan["topology"],
                           "stage_timing_source": "candidate stats YAML; not available in release binary"}
        for binding in plan["bindings"].values():
            # These two live validator/test files are historical context only;
            # neither supplies runtime values or sequence records to this analysis.
            # Immutable compiled source, timing harnesses, inputs, and all raw
            # measurements below must still match their original receipts.
            relative = str(Path(binding["path"]).relative_to(ROOT)) if Path(binding["path"]).is_relative_to(ROOT) else binding["path"]
            if relative in ("scripts/sharkmer_validate/uncertain.py", "tests/test_uncertain_output.py"):
                raw = read(binding["path"])
                if digest(raw) != binding["sha256"]:
                    changed_context_bindings.append({"phase": phase, "path": binding["path"],
                        "historical_sha256": binding["sha256"], "current_sha256": digest(raw),
                        "scope": "Live validator/test changed after historical execution; not used as a source of this analysis or as compiled source."})
            else:
                read(binding["path"], binding)
        for accession, taxon in SAMPLES.items():
            for role in ("baseline", "candidate"):
                paths = sorted(directory.glob(f"*{accession}_*_{role}.json"))
                if len(paths) != 3:
                    raise ValueError("Need all three saved repeats")
                role_sets = []
                for path in paths:
                    job = json.loads(read(path, old_bound.get(str(path))))
                    if job["status"] != "pass" or job["execution"]["returncode"] != 0:
                        raise ValueError("Saved CLI invocation failed")
                    invocation, outputs = job["invocation"], job["outputs"]
                    execution = job["execution"]
                    for file in execution["files"].values():
                        read(file["path"], file)
                    timing_spec = execution["files"]["time.json"]
                    timing = json.loads(read(timing_spec["path"], timing_spec))
                    metrics = execution["metrics"]
                    if any(timing[k] != metrics[k] for k in ("wall_seconds", "user_seconds", "system_seconds", "max_rss_kib", "exit_code")):
                        raise ValueError("GNU time receipt differs from reported metrics")
                    constraint_spec = execution["files"]["process_constraints.json"]
                    constraints = json.loads(read(constraint_spec["path"], constraint_spec))
                    if constraints["status"] != "pass" or constraints["affinity"] != plan["cpus"]:
                        raise ValueError("Saved affinity constraint differs")
                    stats_file = next(f for n, f in outputs["files"].items() if n.endswith(".stats.yaml"))
                    stats_bytes = read(stats_file["path"], stats_file)
                    # Stage fields precede large per-target diagnostic lists.
                    prefix = stats_bytes.decode().split("\npcr_results:", 1)[0]
                    stats = yaml.safe_load(prefix)
                    stages = stats.get("stage_timings") or {}
                    for key in ("n_reads_read", "n_bases_read", "n_subreads_ingested", "n_bases_ingested", "n_kmers"):
                        if stats[key] != outputs["counts"][key]:
                            raise ValueError("Statistics counts differ from raw receipt")
                    binary = Path(invocation["command"][0])
                    binary_hash = digest(read(binary))
                    products = sorted((r["gene"], r["sha256"], r.get("assembly_status", "standard"), r["length"]) for r in outputs["records"])
                    for record in outputs["records"]:
                        if digest(record["sequence"].encode()) != record["sha256"]:
                            raise ValueError("Sequence hash differs")
                    role_sets.append((products, outputs["counts"]))
                    fasta = [f for n, f in outputs["files"].items() if n.endswith(".fasta")]
                    for f in fasta:
                        read(f["path"], f)
                    pipeline = stages.get("pipeline_total_s")
                    row = {"phase": phase, "accession": accession, "taxon": taxon, "role": role,
                           "binary_label": spec["binary_labels"][role], "binary_sha256": binary_hash,
                           "pair": invocation["pair_index"], "order_position": invocation["order_position"],
                           "wall_seconds": timing["wall_seconds"], "user_seconds": timing["user_seconds"],
                           "system_seconds": timing["system_seconds"], "cpu_seconds": timing["user_seconds"] + timing["system_seconds"],
                           "peak_rss_mib": timing["max_rss_kib"] / 1024,
                           "allocator_peak_bytes": outputs["allocator_peak_bytes"],
                           "read_ingest_seconds": stages.get("read_ingest_s"),
                           "count_finalize_seconds": stages.get("count_finalize_s"),
                           "pcr_including_fasta_staging_seconds": stages.get("pcr_s"),
                           "pipeline_total_seconds": pipeline,
                           "wall_minus_pipeline_seconds": timing["wall_seconds"] - pipeline if pipeline is not None else None,
                           "standard_products": sum(p[2] == "standard" for p in products),
                           "uncertain_products": sum(p[2] == "uncertain" for p in products),
                           "fasta_files": len(fasta), "fasta_bytes": sum(f["size_bytes"] for f in fasta),
                           "stats_bytes": len(stats_bytes), "input_records": outputs["counts"]["n_reads_read"],
                           "input_prefix_sha256": job["prefix_before"]["sha256"],
                           "sequence_set_sha256": digest(json.dumps(products, separators=(",", ":")).encode()),
                           "slurm_job_id": state["slurm_job_id"], "node": state["machine"]["node"],
                           "affinity": constraints["affinity"], "invocation_unix_time": constraints["unix_time"],
                           "receipt": str(path)}
                    jobs.append(row)
                if role_sets[1:] != [role_sets[0], role_sets[0]]:
                    raise ValueError("Saved repeat products/counts are unstable")
                sets[(phase, accession, role)] = role_sets[0]
    summaries, paired = [], []
    metrics = ("wall_seconds", "cpu_seconds", "peak_rss_mib", "read_ingest_seconds", "count_finalize_seconds",
               "pcr_including_fasta_staging_seconds", "pipeline_total_seconds", "wall_minus_pipeline_seconds")
    for phase in PHASES:
        for accession, taxon in SAMPLES.items():
            groups = {role: sorted((r for r in jobs if r["phase"] == phase and r["accession"] == accession and r["role"] == role), key=lambda r: r["pair"])
                      for role in ("baseline", "candidate")}
            for first, second in zip(groups["baseline"], groups["candidate"]):
                paired.append({"phase": phase, "accession": accession, "pair": first["pair"],
                               "wall_delta_s": second["wall_seconds"] - first["wall_seconds"],
                               "cpu_delta_s": second["cpu_seconds"] - first["cpu_seconds"],
                               "wall_change_pct": 100 * (second["wall_seconds"] / first["wall_seconds"] - 1)})
            for role, records in groups.items():
                item = {"phase": phase, "accession": accession, "taxon": taxon, "role": role,
                        "binary_label": records[0]["binary_label"], "binary_sha256": records[0]["binary_sha256"],
                        "standard_products": records[0]["standard_products"], "uncertain_products": records[0]["uncertain_products"],
                        "fasta_files": records[0]["fasta_files"], "fasta_bytes": records[0]["fasta_bytes"],
                        "input_records": records[0]["input_records"], "stats_bytes_median": statistics.median(r["stats_bytes"] for r in records)}
                for metric in metrics:
                    values = [r[metric] for r in records]
                    item[metric + "_median"] = statistics.median(values) if all(v is not None for v in values) else None
                    item[metric + "_values"] = values
                summaries.append(item)
    bridges = []
    for accession in SAMPLES:
        left = next(r for r in summaries if r["phase"] == "release_to_inclusive" and r["accession"] == accession and r["role"] == "candidate")
        right = next(r for r in summaries if r["phase"] == "inclusive_to_retention" and r["accession"] == accession and r["role"] == "baseline")
        if left["binary_sha256"] != right["binary_sha256"] or sets[("release_to_inclusive", accession, "candidate")] != sets[("inclusive_to_retention", accession, "baseline")]:
            raise ValueError("Cross-batch same-binary bridge identity differs")
        bridges.append({"accession": accession, "binary_sha256": left["binary_sha256"],
                        "identical_products_and_input_counts": True,
                        "earlier_wall_median_s": left["wall_seconds_median"], "later_wall_median_s": right["wall_seconds_median"],
                        "earlier_cpu_median_s": left["cpu_seconds_median"], "later_cpu_median_s": right["cpu_seconds_median"],
                        "scope": "Different nodes and times; descriptive bridge, not controlled causal comparison"})
    for name in ("benchmarks/benchmark_results/v3.2-inclusive-release-20260929/run_release_comparison.py",
                 "benchmarks/benchmark_results/v3.2-sequence-first-retention-20261001/run_comparison.py",
                 "tmp/inclusive-release-20260929/build-v1/baseline/source/src/stats.rs",
                 "tmp/inclusive-release-20260929/build-v1/candidate/source/src/stats.rs",
                 "tmp/inclusive-release-20260929/build-v1/candidate/source/src/main.rs",
                 "tmp/inclusive-release-20260929/build-v1/candidate/source/src/output.rs"):
        read(ROOT / name)
    limitations = [
        "GNU user+system CPU sums all threads and may exceed wall time; wall minus CPU is not a blocked-time measurement.",
        "pcr_s includes parallel assembly and sequential FASTA staging; it is not an isolated assembly timer.",
        "pipeline_total_s ends before statistics serialization, final transaction commit/reporting and process cleanup. Wall minus pipeline is not a measured fsync duration.",
        "Every historical invocation performs full-prefix inspection and input fingerprinting before timing, which can warm the file cache; there was no cache-drop or verified cold-cache design.",
        "All batches used local uncompressed input paths and output under the shared repository path. Historical syscall, mount-state, storage-load and scheduling profiles were not recorded.",
        "The two batches used different compute nodes, times and CPU affinities. Cross-batch same-binary observations cannot be multiplied into a controlled release-to-current speed ratio.",
        "Additional uncertain records and stats payload are real workload changes, but these timings do not isolate retention, graph search or publication costs.",
        "Three pairs per cell are descriptive evidence, not equivalence tests, release acceptance or portable performance guarantees.",
        "All four original runtime flags remain historical findings; later measurements do not erase them. No new assembly ran in this analysis."]
    payloads = {"saved_runs.tsv": tsv(jobs), "saved_summary.tsv": tsv(summaries), "paired_deltas.tsv": tsv(paired),
                "same_binary_across_batches.tsv": tsv(bridges)}
    result = {"status": "pass", "invocations": len(jobs), "samples": SAMPLES, "machines": machines,
              "runs": jobs, "summaries": summaries, "same_binary_bridges": bridges, "limitations": limitations,
              "changed_historical_context_bindings": changed_context_bindings,
              "producer_sha256": digest(Path(__file__).read_bytes()), "source_files": sources,
              "outputs": {n: {"sha256": digest(v), "size_bytes": len(v)} for n, v in payloads.items()}}
    if len(jobs) != 48 or len(summaries) != 16:
        raise ValueError("Saved comparison scope incomplete")
    payloads["saved_analysis.json"] = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    for name, contents in payloads.items():
        with (HERE / name).open("xb") as stream:
            stream.write(contents)
    print(json.dumps({"status": "pass", "invocations": len(jobs), "same_binary_bridges": bridges}))


if __name__ == "__main__":
    main()
