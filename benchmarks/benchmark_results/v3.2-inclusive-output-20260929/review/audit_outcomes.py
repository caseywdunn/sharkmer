#!/usr/bin/env python3
"""Independent recapture of six inclusive-output regressions; no raw-read scan."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shlex
import sys

import yaml

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "scripts"))
from sharkmer_validate import runner


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return {"path": str(path), "sha256": digest.hexdigest(), "size_bytes": path.stat().st_size}


def verify(row):
    require(fingerprint(row["path"]) == row, "Identity changed: " + row["path"])


def records(path):
    result = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            result.append({"header": line[1:], "sequence": ""})
        elif line.strip():
            require(bool(result), "Sequence before header")
            result[-1]["sequence"] += line.strip()
    require(all(row["sequence"] for row in result), "Empty sequence")
    return result


def normalize_assessment(value):
    value = copy.deepcopy(value)
    value.pop("output_file", None)
    for row in value["products"]:
        row.pop("id")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "Output already exists")
    verification_identity = fingerprint(args.verification)
    value = json.loads(args.verification.read_text())
    require(value["status"] == "pass", "Verification is not passing")
    require(value["source_before"] == value["source_after"], "Source bookends differ")
    bound = {str(args.verification.resolve()): verification_identity}
    for group in (value["source_before"], value["extra_sources"], value["compiler"]):
        for row in group.values():
            verify(row)
            bound[row["path"]] = row
    verify(value["binary"])
    verify(value["checks"])
    bound[value["checks"]["path"]] = value["checks"]
    checks = json.loads(Path(value["checks"]["path"]).read_text())
    require(checks["status"] == "pass", "Compiled checks are not passing")
    require(checks["sources_before"] == checks["sources_after"] == value["source_before"], "Compiled source differs")
    names = ["format", "clippy-default", "release-tests-default", "clippy-fxhashmap", "release-tests-fxhashmap", "python-tests", "reference-audit"]
    require([row["name"] for row in checks["steps"]] == names, "Check steps differ")
    test_totals = {}
    for step in checks["steps"]:
        require(step["status"] == "pass" and step["returncode"] == 0, "Failed check")
        verify(step["log"])
        bound[step["log"]["path"]] = step["log"]
        log = Path(step["log"]["path"]).read_text()
        if step["name"].startswith("release-tests"):
            import re
            totals = [int(n) for n in re.findall(r"test result: ok\. (\d+) passed; 0 failed", log)]
            require(totals == [247, 23, 6], "Unexpected Rust test totals: " + str(totals))
            test_totals[step["name"]] = {"unit": totals[0], "integration": sum(totals[1:3])}
        elif step["name"] == "python-tests":
            require("Ran 197 tests" in log and "OK (skipped=1)" in log, "Unexpected Python test result")
            test_totals[step["name"]] = {"run": 197, "skipped": 1}
    for name in ("build_step", "checks_step", "comparison_controls"):
        require(value[name]["status"] == "pass" and value[name]["returncode"] == 0, "Failed execution step")
        verify(value[name]["log"])
        bound[value[name]["log"]["path"]] = value[name]["log"]
    expected = {"SRR27962769": (11, 17), "SRR31887760": (9, 4), "SRR1057608": (6, 13)}
    expected_jobs = {(sample, mode) for sample in expected for mode in ("default", "standard_only")}
    require(len(value["jobs"]) == 6 and {(j["accession"], j["mode"]) for j in value["jobs"]} == expected_jobs, "Wrong jobs")
    output, recaptured = [], {}
    for job in value["jobs"]:
        sample, mode = job["accession"], job["mode"]
        require(job["status"] == "pass" and job["execution"]["returncode"] == 0, "Failed data run")
        require(job["prefix_before"] == job["prefix_after"], "Prefix evidence differs")
        verify(job["receipt"])
        bound[job["receipt"]["path"]] = job["receipt"]
        saved = json.loads(Path(job["receipt"]["path"]).read_text())
        require(saved == {k: v for k, v in job.items() if k != "receipt"}, "Per-job receipt differs")
        verify(job["prior_receipt"])
        bound[job["prior_receipt"]["path"]] = job["prior_receipt"]
        prior = json.loads(Path(job["prior_receipt"]["path"]).read_text())
        require(prior["status"] == "complete", "Prior run incomplete")
        for record in [*job["outputs"]["files"].values(), *prior["outputs"]["files"].values()]:
            verify(record)
            bound[record["path"]] = record
        directory = Path(job["outputs"]["stats"]["path"]).parent
        stats = yaml.safe_load(Path(job["outputs"]["stats"]["path"]).read_text())
        manifest = yaml.safe_load(Path(job["outputs"]["manifest"]["path"]).read_text())
        constraints = json.loads(Path(job["outputs"]["files"]["process_constraints.json"]["path"]).read_text())
        require(constraints["status"] == "pass" and constraints["command"] == job["command"], "Constraint command mismatch")
        require(constraints["affinity"] == value["affinity"] and len(value["affinity"]) == 2, "Wrong CPU affinity")
        require(constraints["rlimit_as_soft_hard"] == [42949672960, 42949672960], "Wrong address-space limit")
        require(stats["product_output_schema_version"] == manifest["schema_version"] == 2, "Wrong schema")
        require(stats["product_output_layout"] == "combined_classified", "Wrong layout")
        require(stats["stage_timings"]["read_threading_input_s"] is None, "Unexpected threading")
        require(("--no-emit-uncertain" in job["command"]) == (mode == "standard_only"), "Wrong delivery option")
        require("--emit-uncertain" not in job["command"], "Default regression used opt-in")
        old_stats = prior["outputs"]["stats"]
        old_directory = Path(prior["outputs"]["files"][stats["sample"] + ".stats.yaml"]["path"]).parent
        require(yaml.safe_load((old_directory / (stats["sample"] + ".stats.yaml")).read_text()) == old_stats, "Prior inline stats differ")
        for field in ("n_reads_read", "n_bases_read", "n_subreads_ingested", "n_bases_ingested", "n_kmers", "n_multi_kmers", "n_singleton_kmers"):
            require(stats.get(field) == old_stats.get(field) and (field in stats) == (field in old_stats), "Count changed: " + field)
        old_genes = {row["gene_name"]: row for row in old_stats["pcr_results"]}
        require(len(stats["pcr_results"]) == len(old_genes) == 21, "Target count differs")
        standard, uncertain, genes, identifiers = [], [], [], []
        for entry in stats["pcr_results"]:
            gene = entry["gene_name"]
            old = old_genes[gene]
            new_records = records(directory / entry["output_file"]) if entry.get("output_file") else []
            count = entry["n_standard_products"]
            std, unc = new_records[:count], new_records[count:]
            old_std = records(old_directory / old["output_file"]) if old.get("output_file") else []
            old_unc = records(old_directory / old["uncertain_products"]["output_file"]) if old["uncertain_products"].get("output_file") else []
            require(entry["n_products"] == len(new_records) == count + entry["n_uncertain_products"], "Inclusive count mismatch")
            require(entry.get("product_lengths", []) == [len(row["sequence"]) for row in new_records], "Inclusive lengths mismatch")
            require(entry["status"] == ("success" if std else "uncertain" if unc else "fail"), "Status mismatch")
            require(entry.get("failure_reason") == old.get("failure_reason"), "Failure reason changed")
            require(entry.get("threshold_diagnostics", []) == old.get("threshold_diagnostics", []), "Threshold evidence changed")
            require(len(std) == old["n_products"] == len(old_std), "Standard count changed")
            for row, old_row in zip(std, old_std):
                require(row["header"].endswith(" assembly_status=standard"), "Standard suffix missing")
                require(row["header"][:-len(" assembly_status=standard")] == old_row["header"] and row["sequence"] == old_row["sequence"], "Literal standard record changed")
            if mode == "default":
                require(len(unc) == len(old_unc), "Uncertain count changed")
                require(normalize_assessment(entry["uncertain_products"]) == normalize_assessment(old["uncertain_products"]), "Assessment changed beyond ID/path")
                for index, (row, old_row) in enumerate(zip(unc, old_unc)):
                    require(row["header"].split(" ", 1)[1] == old_row["header"].split(" ", 1)[1] and row["sequence"] == old_row["sequence"], "Uncertain sequence/description changed")
                    require(row["header"].split()[0] == f'{stats["sample"]}_{gene}_{index}.uncertain', "Uncertain ID differs")
            else:
                require(not unc and "uncertain_products" not in entry, "Opt-out retained uncertainty")
            for index, row in enumerate(new_records):
                labels = [token for token in shlex.split(row["header"]) if token.startswith("assembly_status=")]
                require(labels == ["assembly_status=" + ("standard" if index < count else "uncertain")], "Record classification differs")
                identifiers.append(row["header"].split()[0])
            standard.extend({"gene": gene, **row} for row in std)
            uncertain.extend({"gene": gene, **row} for row in unc)
            genes.append({"gene": gene, "standard": len(std), "uncertain": len(unc), "total": len(new_records), "status": entry["status"]})
        require(len(identifiers) == len(set(identifiers)), "Cross-target ID collision")
        require(standard == job["outputs"]["standard"] and uncertain == job["outputs"]["uncertain"] and genes == job["outputs"]["genes"], "Captured summary differs from raw FASTA/stats")
        require((len(standard), len(uncertain)) == (expected[sample][0], expected[sample][1] if mode == "default" else 0), "Unexpected output totals")
        consumed = runner.parse_fasta_products(stats["sample"], directory)
        require([(g["gene"], r["header"], r["sequence"]) for g in consumed for r in g["products"]] == [(r["gene"], r["header"], r["sequence"]) for r in standard], "Current consumer did not select exactly standard records")
        output.append({"accession": sample, "mode": mode, "standard": len(standard), "uncertain": len(uncertain), "total": len(standard) + len(uncertain), "targets": 21})
        recaptured[(sample, mode)] = uncertain
    historical = []
    for loss in value["historical_losses"]:
        matches = [r["header"] for r in recaptured[(loss["accession"], "default")] if r["gene"] == "insecta_" + loss["gene"] and hashlib.sha256(r["sequence"].encode()).hexdigest() == loss["sequence_sha256"]]
        require(matches == loss["exact_matches"], "Historical match summary differs")
        historical.append({**loss, "exact_matches": matches})
    require(len(historical) == 7 and sum(bool(row["exact_matches"]) for row in historical) == 5, "Historical result differs")
    for row in bound.values():
        verify(row)
    result = {"schema_version": 1, "status": "pass", "verification": verification_identity,
              "checker": fingerprint(__file__), "independent_raw_output_recapture": True,
              "raw_reads_rescanned_by_this_review": False, "run_summary": output,
              "historical_losses": historical, "test_totals": test_totals,
              "reviewed_artifacts": bound,
              "scope": "Actual FASTA/stats recapture and exact old/current comparison; consumer standard subset; checked source, executable, constraints and receipt bindings. Input-prefix identities are inherited from the reviewed producer's before/after scans, not rescanned here."}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": "pass", "runs": 6, "historical_exact": 5, "output": fingerprint(args.output)}))


if __name__ == "__main__":
    main()
