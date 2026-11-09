#!/usr/bin/env python3
"""Independently review final-check receipts, raw outputs, and resource results.

This reviewer does not import the execution harness or Sharkmer validators.
It reads small retained artifacts and bound sources, never the large FASTQ inputs.
"""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import shlex
import sys

import yaml


COUNTS = ("n_reads_read", "n_bases_read", "n_subreads_ingested", "n_bases_ingested", "n_kmers")
PERFORMANCE = {"SRR17535371", "SRR8617500", "SRR25378184", "SRR1057608"}
PLANTS = {"SRR25378184", "ERR14009273"}
CAUSES = {
    "omitted_self_loop_marker_occurrences": "homopolymer_length_unresolved",
    "cyclic_scc_marker_occurrences": "repeat_structure_unresolved",
    "retained_collision_marker_occurrences": "high_coverage_repeat_suspected",
}


def check(condition, message):
    if not condition:
        raise ValueError(message)


class UniqueLoader(yaml.SafeLoader):
    pass


def unique_mapping(loader, node, deep=False):
    pairs = loader.construct_pairs(node, deep=deep)
    result = {}
    for key, value in pairs:
        check(key not in result, f"Duplicate YAML key: {key}")
        result[key] = value
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)


def read_json(path):
    return json.loads(Path(path).read_text())


def read_yaml(path):
    return yaml.load(Path(path).read_text(), Loader=UniqueLoader)


def fingerprint(path):
    path = Path(path).resolve()
    check(path.suffix not in {".fastq", ".fq"}, "Reviewer must not read large FASTQ inputs")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "sha256": digest.hexdigest(), "size_bytes": path.stat().st_size}


def verify(row):
    check(fingerprint(row["path"]) == row, f"Fingerprint mismatch: {row['path']}")


def read_fasta(path):
    records = []
    header, bases = None, []
    with Path(path).open() as stream:
        for raw in stream:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(bases)))
                header, bases = line[1:], []
                check(header and not header[0].isspace(), "Empty or malformed FASTA identifier")
            else:
                check(header is not None and set(line) <= set("ACGT"), "Invalid FASTA bases")
                bases.append(line)
    if header is not None:
        records.append((header, "".join(bases)))
    check(all(sequence for _, sequence in records), "Empty FASTA sequence")
    return records


def annotations(header):
    identifier, *description = header.split(maxsplit=1)
    pairs = [token.split("=", 1) for token in shlex.split(description[0] if description else "")]
    check(all(len(pair) == 2 and pair[0] for pair in pairs), "Malformed annotation token")
    check(len({key for key, _ in pairs}) == len(pairs), "Duplicate annotation key")
    return identifier, dict(pairs)


def inspect_outputs(products, sample, role, command, version):
    prefix = sample["sample_prefix"]
    expected = sample["expected"][role]
    stats_name = prefix + ".stats.yaml"
    stats = read_yaml(products / stats_name)
    check(stats["sample"] == prefix and stats["command"] == " ".join(command), "Stats invocation mismatch")
    check(stats["sharkmer_version"] == version and stats["kmer_length"] == 19
          and stats["chunks"] == 0, "Stats version/settings mismatch")
    observed_counts = {key: stats[key] for key in COUNTS}
    check(observed_counts == expected["counts"], "Raw stats count oracle mismatch")
    check(observed_counts["n_reads_read"] == sample["expected_prefix"]["records"]
          and observed_counts["n_bases_read"] == sample["expected_prefix"]["bases"], "Read-prefix count mismatch")
    check(len(stats["pcr_results"]) == len(expected["targets"]), "Wrong target count")
    manifest_files = {stats_name}
    records, targets, all_ids = [], [], set()
    for entry, old_target in zip(stats["pcr_results"], expected["targets"]):
        gene = entry["gene_name"]
        count = entry["n_products"]
        target = {"gene_name": gene, "status": entry["status"], "n_products": count,
                  "n_standard_products": entry.get("n_standard_products", count),
                  "n_uncertain_products": entry.get("n_uncertain_products", 0),
                  "product_lengths": entry.get("product_lengths", [])}
        check(target == old_target, f"Raw target outcome/order oracle mismatch: {gene}")
        targets.append(target)
        filename = prefix + "_" + gene + ".fasta"
        path = products / filename
        check(path.exists() == bool(count), f"Unexpected or absent FASTA: {filename}")
        raw_records = read_fasta(path) if count else []
        if count:
            manifest_files.add(filename)
        check(len(raw_records) == count and [len(sequence) for _, sequence in raw_records]
              == target["product_lengths"], "FASTA count/length disagrees with stats")
        warning_rows = entry.get("uncertain_products", {}).get("products", [])
        warning_by_id = {row["id"]: row for row in warning_rows}
        check(len(warning_by_id) == len(warning_rows), "Duplicate structured warning IDs")
        classes = []
        for header, sequence in raw_records:
            identifier, fields = annotations(header)
            check(identifier not in all_ids, "Duplicate emitted identifier")
            all_ids.add(identifier)
            check("read_support" not in fields, "Unexpected unevaluated read-support annotation")
            warning = fields.get("assembly_status") == "warning"
            check(fields.get("assembly_status") in (None, "warning"), "Unexpected status spelling")
            check(not warning or role == "candidate", "Warning in released output")
            label = "uncertain" if warning else "standard"
            classes.append(label)
            digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
            record = {"gene": gene, "id": identifier, "length": len(sequence), "sha256": digest,
                      "assembly_status": label}
            if warning:
                check(set(fields) == {"sample", "gene", "product", "length", "assembly_status", "assembly_status_causes"},
                      "Warning header field dictionary differs")
                row = warning_by_id[identifier]
                causes = sorted({name for key, name in CAUSES.items()
                                 if any(source[key] > 0 for source in row["sources"])})
                check(causes and fields["assembly_status_causes"] == ",".join(causes)
                      and row["assembly_status_causes"] == causes, "Warning cause evidence mismatch")
                check(row["oriented_sequence"] == sequence and row["sequence_sha256"] == digest
                      and row["sequence_length"] == len(sequence)
                      and row["description"] == header.split(maxsplit=1)[1], "Warning record/evidence join mismatch")
                check(fields["sample"] == prefix and fields["gene"] == gene
                      and fields["length"] == str(len(sequence)), "Warning identity/length annotation mismatch")
            else:
                check("assembly_status_causes" not in fields and identifier not in warning_by_id,
                      "Ordinary record carries warning evidence")
                record["ordinary_fields"] = fields
            records.append(record)
        check(classes == ["standard"] * target["n_standard_products"]
              + ["uncertain"] * target["n_uncertain_products"], "FASTA class counts/order mismatch")
        check(len(warning_rows) == target["n_uncertain_products"], "Structured warning count mismatch")
    check(records == expected["records"], "Exact sequence/identifier/order/header/class oracle mismatch")
    check({p.name for p in products.glob("*.fasta")} == manifest_files - {stats_name}, "Unreported FASTA output")
    if role == "candidate":
        check(stats["product_output_schema_version"] == 3 and stats["product_output_layout"] == "combined_warnings"
              and stats["run_status"] == "complete", "Wrong current output contract")
        check(stats["input_source"] == {"kind": "local_files", "inputs": [sample["input"]["path"]],
              "paired": False, "max_reads": sample["max_reads"]}, "Input provenance mismatch")
        check(stats["stage_timings"]["read_threading_input_s"] is None, "Unexpected read threading")
        manifest_name = prefix + ".manifest.yaml"
        check(stats["output_manifest"] == manifest_name, "Stats manifest link mismatch")
        manifest = read_yaml(products / manifest_name)
        check(manifest["schema_version"] == 3 and manifest["producer"] == "sharkmer"
              and manifest["sample"] == prefix and manifest["run_id"] == stats["run_id"]
              and manifest["status"] == "complete", "Manifest generation mismatch")
        paths = [row["path"] for row in manifest["files"]]
        check(len(paths) == len(set(paths)) and set(paths) == manifest_files, "Manifest file inventory mismatch")
        for row in manifest["files"]:
            check(Path(row["path"]).name == row["path"] and fingerprint(products / row["path"])["sha256"]
                  == row["sha256"], "Manifest output hash mismatch")
    else:
        check(not any(key in stats for key in ("run_id", "run_status", "product_output_schema_version")),
              "Released stats falsely carries current schema")
    return {"counts": observed_counts, "targets": targets, "records": records,
            "allocator_peak_bytes": stats["peak_memory_bytes"], "stage_timings": stats.get("stage_timings")}


def metrics(path):
    raw = read_json(path)
    check(raw["exit_code"] == 0 and all(isinstance(v, (int, float)) and math.isfinite(v) and v >= 0
          for v in raw.values()), "Invalid raw GNU time metrics")
    return dict(raw, cpu_seconds=raw["user_seconds"] + raw["system_seconds"],
                peak_rss_bytes=raw["max_rss_kib"] * 1024)


def median(values):
    ordered = sorted(values)
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2


def same_numbers(left, right):
    if isinstance(left, dict):
        return isinstance(right, dict) and left.keys() == right.keys() and all(same_numbers(left[k], right[k]) for k in left)
    if isinstance(left, list):
        return isinstance(right, list) and len(left) == len(right) and all(same_numbers(a, b) for a, b in zip(left, right))
    if isinstance(left, (int, float)) and not isinstance(left, bool):
        return isinstance(right, (int, float)) and math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-10)
    return left == right


def analyze(runs, protocol):
    cells = []
    for sample in protocol["performance_samples"]:
        selected = {(run["invocation"]["pair_index"], run["invocation"]["role"]): run
                    for run in runs if run["invocation"]["purpose"] == "performance"
                    and run["invocation"]["sample"] == sample}
        check(set(selected) == {(pair, role) for pair in (1, 2, 3) for role in ("baseline", "candidate")},
              "Missing timing pair")
        resources = {}
        for role in ("baseline", "candidate"):
            reference = selected[(1, role)]["outputs"]
            check(all(selected[(pair, role)]["outputs"] == reference for pair in (2, 3)), "Repeated outputs unstable")
            resources[role] = {}
            for name in ("wall_seconds", "cpu_seconds", "user_seconds", "system_seconds", "peak_rss_bytes"):
                values = [selected[(pair, role)]["metrics"][name] for pair in (1, 2, 3)]
                resources[role][name] = {"values": values, "median": median(values), "min": min(values), "max": max(values)}
        old = {(r["gene"], r["sha256"]) for r in selected[(1, "baseline")]["outputs"]["records"]}
        new = {(r["gene"], r["sha256"]) for r in selected[(1, "candidate")]["outputs"]["records"]}
        check(old <= new and selected[(1, "baseline")]["outputs"]["counts"]
              == selected[(1, "candidate")]["outputs"]["counts"], "Pair count/retention mismatch")
        changes, pairs = {}, []
        for name in ("wall_seconds", "cpu_seconds", "peak_rss_bytes"):
            before, after = resources["baseline"][name]["median"], resources["candidate"][name]["median"]
            percent = (after - before) * 100 / before if before else None
            threshold = {"wall_seconds": 5, "peak_rss_bytes": 2}.get(name)
            changes[name] = {"absolute": after - before, "percent": percent,
                             "investigation_flag": threshold is not None and percent is not None and percent > threshold}
        for pair in (1, 2, 3):
            pairs.append({"pair_index": pair, **{name: selected[(pair, "candidate")]["metrics"][name]
                          - selected[(pair, "baseline")]["metrics"][name]
                          for name in ("wall_seconds", "cpu_seconds", "peak_rss_bytes")}})
        cells.append({"sample": sample, "released_sequences": len(old), "retained_sequences": len(old & new),
                      "resources": resources, "changes": changes, "paired_deltas": pairs})
    candidates = {run["invocation"]["sample"]: run["outputs"] for run in runs if run["invocation"]["role"] == "candidate"}
    check(len(candidates) == 12, "Missing candidate inventory")
    totals = {}
    for label, keys in (("historical_panel_inventory", [key for key in candidates if not key.endswith("_corrected")]),
                        ("corrected_panel_inventory", [key for key in candidates if key not in PLANTS])):
        classes = Counter(r["assembly_status"] for key in keys for r in candidates[key]["records"])
        totals[label] = {"records": sum(classes.values()), "standard": classes["standard"], "warning": classes["uncertain"]}
    check(totals["historical_panel_inventory"] == {"records": 326, "standard": 64, "warning": 262}
          and totals["corrected_panel_inventory"] == {"records": 287, "standard": 66, "warning": 221}, "Inventory totals differ")
    retained, with_warning = 0, 0
    for key in candidates:
        if key.endswith("_corrected"):
            continue
        old = {(r["gene"], r["sha256"]) for r in protocol["samples"][key]["expected"]["baseline"]["records"]}
        current = candidates[key + "_corrected" if key in PLANTS else key]["records"]
        present = {(r["gene"], r["sha256"]) for r in current}
        warning = {(r["gene"], r["sha256"]) for r in current if r["assembly_status"] == "uncertain"}
        check(old <= present, "Current inventory loses released identity")
        retained += len(old)
        with_warning += len(old & warning)
    check((retained, with_warning) == (86, 22), "Released identity/class totals differ")
    return {"cells": cells, **totals, "invocations": len(runs), "released_identities": retained,
            "retained_identities": retained, "retained_with_warning": with_warning, "lost_identities": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    check(not args.output.exists(), "Review destination must be new")
    report = {"status": "failed", "reviewer": fingerprint(__file__), "python": sys.version,
              "scope": "Independent raw-artifact review; FASTQ hash/count checks are attested by frozen runner receipts, not repeated here."}
    try:
        directory = args.results.resolve()
        result = read_json(directory / "result.json")
        check(result["status"] == "pass" and len(result["runs"]) == 32, "Run did not finish all 32 invocations")
        verify(result["protocol"])
        protocol = read_json(result["protocol"]["path"])
        report["protocol"] = result["protocol"]
        report["slurm_job_id"] = result["slurm_job_id"]
        check(set(protocol["performance_samples"]) == PERFORMANCE and set(protocol["corrected_replaces"]) == PLANTS,
              "Unexpected performance or assay scope")
        check(result["bindings_before"] == result["bindings_after"] == protocol["bindings"], "Frozen source bindings changed")
        for row in protocol["bindings"]:
            verify(row)
        report["source_binding_checks"] = len(protocol["bindings"])
        preflight = read_json(directory / "preflight.json")
        check(preflight["bindings_before"] == result["bindings_before"] and preflight["cpus"] == result["cpus"],
              "Preflight bindings/CPU evidence differs")
        check(preflight["input_preflight"] == result["input_preflight"]
              and len(result["input_preflight"]) == 10, "Input prefix preflight evidence differs")
        for sample in protocol["samples"].values():
            check(sample["max_reads"] == 1000000
                  and result["input_preflight"][sample["input"]["path"]] == sample["expected_prefix"],
                  "Frozen sample prefix differs from checked input")
        cpus = result["cpus"]
        check(len(cpus) == len(set(cpus)) == 2 and len({(r["socket"], r["core"]) for r in result["topology"] if r["cpu"] in cpus}) == 2,
              "Expected two distinct physical cores")
        check(result["output_mount"]["filesystems"][0]["fstype"].startswith("nfs"), "Wrong timing output storage")
        check(len(protocol["schedule"]) == len({r["id"] for r in protocol["schedule"]}) == 32, "Invalid schedule")
        check([run["invocation"] for run in result["runs"]] == protocol["schedule"], "Executed schedule differs")
        for pair in (1, 2, 3):
            for sample in PERFORMANCE:
                roles = [row["role"] for row in protocol["schedule"] if row["sample"] == sample and row["pair_index"] == pair]
                check(roles == (["candidate", "baseline"] if pair == 2 else ["baseline", "candidate"]), "Pair order not alternated")
        verified_runs = []
        for original, item in zip(result["runs"], protocol["schedule"]):
            run_dir = directory / item["id"]
            receipt = read_json(run_dir / "receipt.json")
            check(receipt == original and receipt["status"] == "pass", "Result/individual receipt mismatch")
            sample, role = protocol["samples"][item["sample"]], item["role"]
            check(receipt["input_before"] == receipt["input_after"] == sample["input"],
                  "Per-invocation input fingerprint changed")
            products = run_dir / "outputs"
            command = [protocol["binaries"][role]["path"], "-k", "19", "-t", "2", "--chunks", "0",
                       "--max-reads", "1000000", "-o", str(products) + "/", "-s", sample["sample_prefix"],
                       "--pcr-panel-file", sample["panel"]["path"], sample["input"]["path"]]
            execution = receipt["execution"]
            check(execution["command"] == command and execution["returncode"] == 0
                  and not execution.get("timed_out"), "Invocation command/exit differs")
            execution_files = [fingerprint(run_dir / name) for name in
                               ("constraints.json", "stderr.log", "stdout.log", "time.json")]
            check(execution["files"] == execution_files, "Execution artifact receipt set/hash mismatch")
            constraints = read_json(run_dir / "constraints.json")
            check(constraints["status"] == "pass" and constraints["affinity"] == cpus
                  and constraints["rlimit_as_soft_hard"] == [protocol["address_space_bytes"]] * 2,
                  "Process constraints differ")
            check(constraints["command"][0] == "/usr/bin/time" and constraints["command"][-len(command):] == command,
                  "Measured command differs")
            raw_metrics = metrics(run_dir / "time.json")
            check(raw_metrics == execution["metrics"], "Timing receipt differs from raw GNU time")
            observed = inspect_outputs(products, sample, role, command, protocol["versions"][role])
            for key in ("counts", "targets", "allocator_peak_bytes", "stage_timings"):
                check(observed[key] == receipt["outputs"][key], f"Raw/receipt {key} differs")
            normalized = [{key: value for key, value in row.items() if not (key == "ordinary_fields" and value is None)}
                          for row in receipt["outputs"]["records"]]
            check(observed["records"] == normalized, "Raw FASTA differs from receipt inventory")
            actual_files = [fingerprint(path) for path in sorted(products.iterdir()) if path.is_file()]
            check(actual_files == receipt["outputs"]["files"], "Output artifact receipt set/hash mismatch")
            stable = {key: observed[key] for key in ("counts", "targets", "records")}
            verified_runs.append({"invocation": item, "outputs": stable, "metrics": raw_metrics})
        check(Counter(r["invocation"]["purpose"] for r in verified_runs)
              == {"performance": 24, "inventory": 6, "corrected_plant": 2}, "Invocation purpose totals differ")
        analysis = analyze(verified_runs, protocol)
        summary = read_json(directory / "summary.json")
        verify(result["summary"])
        check(summary["status"] == "pass" and all(same_numbers(value, summary[key]) for key, value in analysis.items()),
              "Independent resource/inventory summary differs")
        report.update(analysis=analysis, reviewed_invocations=len(verified_runs), source_unchanged=True,
                      artifacts=[fingerprint(path) for path in sorted(directory.rglob("*")) if path.is_file()
                                 and path.resolve() != args.output.resolve()], status="pass")
    except Exception as error:
        report["error"] = repr(error)
        raise
    finally:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    print(json.dumps({"status": report["status"], "invocations": report["reviewed_invocations"],
                      "artifact_count": len(report["artifacts"]), "output": str(args.output)}))


if __name__ == "__main__":
    main()
