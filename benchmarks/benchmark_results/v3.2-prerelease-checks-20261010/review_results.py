#!/usr/bin/env python3
"""Independently verify saved final checks and plant outputs; launch no tests."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

import yaml

BASE = Path(__file__).resolve().parent
RESULTS = BASE / "results-v1"
ROOT = BASE.parents[2]
receipts = {}


def bind(path, expected=None):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open("rb") as stream:
        while data := stream.read(1024 * 1024):
            h.update(data)
    seen = {"path": str(path), "sha256": h.hexdigest(), "size_bytes": path.stat().st_size}
    if expected:
        assert all(seen[k] == expected[k] for k in expected), str(path)
    receipts[str(path)] = seen
    return seen


def read(path):
    bind(path)
    return json.loads(Path(path).read_text())


def records(directory, sample):
    rows = []
    for path in sorted(directory.glob("*.fasta")):
        bind(path)
        gene = path.stem.removeprefix(sample + "_")
        for block in path.read_text().split(">")[1:]:
            header, *lines = block.splitlines()
            sequence = "".join(lines)
            assert sequence and set(sequence) <= set("ACGT")
            status = re.findall(r"(?:^|\s)assembly_status=(\S+)", header)
            assert len(status) == 1 and status[0] in ("standard", "uncertain")
            rows.append({"gene": gene, "sequence": sequence, "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "length": len(sequence), "assembly_status": status[0], "header": header})
    return rows


def flags(value, path=()):
    found = {}
    if isinstance(value, dict):
        for key, child in value.items():
            if key in {"search_limits_reached", "retention_truncated", "marker_detail_truncated", "marker_detail", "path_search_evaluated", "search_limit_scope"}:
                found[path + (key,)] = child
            found.update(flags(child, path + (key,)))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.update(flags(child, path + (index,)))
    return found


def main():
    result = read(RESULTS / "result.json")
    assert result["status"] == "pass" and result["not_run_gates"] == []
    assert result["slurm_job_id"] == "29147703"
    plan = read(result["plan"]["path"])
    bind(result["plan"]["path"], result["plan"])
    assert result["plan"]["sha256"] == "bda19fa8e6164ed689605713454663a18caa0fdb2f896770aafcc109c4dbf0c0"
    assert result["source_files_unchanged"] and result["bindings_unchanged"]
    assert read(RESULTS / "sources.before.json") == read(RESULTS / "sources.after.json") == plan["source_files"]
    for row in list(plan["source_files"].values()) + plan["bindings"]:
        bind(row["path"], row)
    bind(result["frozen_binary"]["path"], result["frozen_binary"])
    assert read(RESULTS / "frozen_binary.json")["binary"] == result["frozen_binary"]
    assert [g["name"] for g in result["gates"]] == plan["gate_order"]
    assert len({(r["socket"], r["core"]) for r in result["topology"] if r["cpu"] in result["cpus"]}) == 2
    rust_counts = {}
    for gate in result["gates"]:
        directory = RESULTS / "gates" / gate["name"]
        assert read(directory / "receipt.json") == gate
        assert gate["status"] == "pass" and gate["returncode"] == 0
        assert read(directory / "time.json") == gate["metrics"] and gate["metrics"]["exit_code"] == 0
        constraints = read(directory / "constraints.json")
        assert constraints["status"] == "pass" and constraints["affinity"] == result["cpus"]
        assert constraints["rlimit_as_soft_hard"] == [12 * 1024**3] * 2
        assert constraints["command"][-len(gate["command"]):] == gate["command"]
        for row in gate["files"]:
            bind(row["path"], row)
        if gate["name"].startswith("test_"):
            log = (directory / "stdout.log").read_text()
            counts = re.findall(r"test result: ok\. (\d+) passed; 0 failed; 0 ignored; 0 measured; 0 filtered out", log)
            assert counts == ["251", "23", "6"]
            rust_counts[gate["name"]] = list(map(int, counts))
    python = (RESULTS / "gates/python_all/stderr.log").read_text()
    assert "Ran 200 tests" in python and "\nOK\n" in python and "skipped=" not in python
    assert re.search(r"test_real_database_build_and_query_binds_logical_gene_metadata.*\.\.\. ok", python)
    reference = read(RESULTS / "references/verification.json")["audit"]
    assert reference["panel_count"] == 10 and reference["total_references"] == reference["verified_references"] == 158
    assert reference["failure_count"] == reference["excluded_references"] == 0
    assert sum(p["verified_count"] for p in reference["panels"]) == 158
    panel = read(RESULTS / "panels/verification.json")
    assert (panel["schema_checks"], panel["builtin_exports"], panel["loader_dry_runs"]) == (10, 9, 19)
    assert panel["binary"] == result["frozen_binary"]
    for row in panel["panels"]:
        bind(row["panel"]["path"], row["panel"])
        assert row["schema"] == "pass" and all(c["returncode"] == 0 for c in row["commands"])
        export = RESULTS / "panels" / (Path(row["panel"]["path"]).stem + ".export.yaml")
        if export.exists():
            assert yaml.safe_load(export.read_text()) == yaml.safe_load(Path(row["panel"]["path"]).read_text())
    fixture = read(RESULTS / "fixtures/verification.json")
    bind(fixture["fixture"]["path"], fixture["fixture"])
    bind(fixture["oracle"]["path"], fixture["oracle"])
    truth = yaml.safe_load(Path(fixture["oracle"]["path"]).read_text())["samples"][0]
    for k, run in zip((19, 31), fixture["runs"]):
        sample = "fixture_k" + str(k)
        directory = RESULTS / "fixtures" / sample
        stats = yaml.safe_load((directory / (sample + ".stats.yaml")).read_text())
        assert stats["run_status"] == "complete" and stats["command"] == " ".join(run["command"])
        assert {key: stats[key] for key in run["counts"]} == run["counts"]
        raw = records(directory, sample)
        assert [{key: r[key] for key in ("gene", "header", "length", "sha256")} for r in raw] == run["all_records"]
        if k == 19:
            assert stats["n_kmers"] == truth["expected"]["counting-only"]["n_kmers"] == 13197385 and not raw
        else:
            expected = {"cnidaria_" + gene: (rows[0]["length"], rows[0]["sha256"]) for gene, rows in truth["expected"]["end-to-end"]["genes"].items()}
            assert {r["gene"]: (r["length"], r["sha256"]) for r in raw if r["assembly_status"] == "standard"} == expected
    plant = read(RESULTS / "plant/summary.json")
    plant_plan = read(RESULTS / "plant/plan.json")
    assert plant["status"] == "pass" and len(plant["jobs"]) == 8
    assert plant_plan["binary"] == result["frozen_binary"] and plant_plan["cpus"] == [min(result["cpus"])]
    for row in list(plant_plan["bindings"].values()) + list(plant_plan["panels"].values()):
        bind(row["path"], row)
    raw_jobs = {}
    for job in plant["jobs"]:
        directory = Path(job["directory"])
        assert read(directory / "receipt.json") == job and job["status"] == job["execution"]["status"] == "pass"
        bind(job["input_file"]["path"], job["input_file"])
        assert job["input_file"]["sha256"] == job["expected_prefix"]["sha256"]
        for row in list(job["execution"]["files"].values()) + list(job["outputs"]["files"].values()):
            bind(row["path"], row)
        timing = read(directory / "time.json")
        assert all(job["execution"]["metrics"][k] == v for k, v in timing.items()) and timing["exit_code"] == 0
        constraints = read(directory / "process_constraints.json")
        assert constraints["affinity"] == plant_plan["cpus"] and constraints["rlimit_as_soft_hard"] == [job["address_space_bytes"]] * 2
        assert job["address_space_bytes"] == (2 if job["scope"] == "synthetic" else 12) * 1024**3
        output = Path(job["outputs_directory"])
        sample = job["sample_prefix"]
        stats = yaml.safe_load((output / (sample + ".stats.yaml")).read_text())
        assert stats["command"] == " ".join(job["command"]) and stats["run_status"] == "complete"
        assert stats["pcr_results"] == job["outputs"]["pcr_results"] and len(stats["pcr_results"]) == 10
        assert {k: stats[k] for k in job["outputs"]["counts"]} == job["outputs"]["counts"]
        assert stats["n_reads_read"] == job["expected_prefix"]["records"] and stats["n_bases_read"] == job["expected_prefix"]["bases"]
        manifest = yaml.safe_load((output / stats["output_manifest"]).read_text())
        assert manifest["schema_version"] == 2 and manifest["status"] == "complete" and manifest["run_id"] == stats["run_id"]
        for row in manifest["files"]:
            assert bind(output / row["path"])["sha256"] == row["sha256"]
        raw = records(output, sample)
        if raw:
            captured = [{k: r[k] for k in raw[0]} for r in job["outputs"]["records"]]
            order = lambda r: (r["gene"], r["header"])
            assert sorted(raw, key=order) == sorted(captured, key=order)
        else:
            assert not job["outputs"]["records"]
        raw_jobs[(job["scope"], job["sample_id"], job["variant"])] = (raw, stats)
    changed = "angiospermae_trnV-atpE"
    changes = []
    for comparison in plant["comparisons"]:
        old, new = [raw_jobs[(comparison["scope"], comparison["sample_id"], variant)] for variant in ("before", "after")]
        assert all(old[1][k] == new[1][k] for k in comparison["aggregate_counts"]["before"])
        old_targets = {r["gene_name"]: r for r in old[1]["pcr_results"]}
        new_targets = {r["gene_name"]: r for r in new[1]["pcr_results"]}
        for gene in old_targets.keys() - {changed}:
            identity = lambda rows: Counter((r["sha256"], r["assembly_status"]) for r in rows if r["gene"] == gene)
            assert identity(old[0]) == identity(new[0])
            for key in ("status", "n_products", "n_standard_products", "n_uncertain_products", "failure_reason", "threshold_diagnostics"):
                assert old_targets[gene].get(key) == new_targets[gene].get(key)
            assert flags(old_targets[gene].get("uncertain_products")) == flags(new_targets[gene].get("uncertain_products"))
        after = [r for r in new[0] if r["gene"] == changed]
        before = [r for r in old[0] if r["gene"] == changed]
        assert len(after) == 1 and after[0]["assembly_status"] == "standard"
        expected = (909, "8150cf0b4ec2c207ea6f7a45c85ffed36be9c74a12bdf4e3188bec739f1f5379") if comparison["sample_id"] in ("SRR25378184", "NC_008326.1") else (929, "60a67ee9a9fac06519b8aff377abc4f2ab9e252b12121adf28ff37421e9dfc9a")
        assert (after[0]["length"], after[0]["sha256"]) == expected
        assert len(before) == (41 if comparison["sample_id"] == "SRR25378184" else 0)
        assert all(r["assembly_status"] == "uncertain" for r in before)
        changes.append({"scope": comparison["scope"], "sample": comparison["sample_id"], "before_uncertain": len(before), "after_standard": 1, "after_length": expected[0], "after_sha256": expected[1], "unchanged_control_targets": 9})
    for path in sorted(RESULTS.rglob("*")):
        if path.is_file():
            bind(path)
    bind(Path(__file__))
    review = {"status": "pass", "slurm_job_id": "29147703", "gates": 11, "rust_per_backend": rust_counts,
        "python_tests": 200, "python_skips": 0, "references": 158, "schema_panels": 10, "loader_checks": 19,
        "plant_jobs": 8, "plant_comparisons": changes, "frozen_binary": result["frozen_binary"], "source_receipts": receipts,
        "method": "Independent checks of saved raw files, literal sequences, outcomes, counts, diagnostics, constraints and receipts; no test or assembly rerun.",
        "limits": ["Partial-reference assessment of legacy changed-target products is reviewed separately when its post-run report is complete.", "Versioned assay correction changes intended endpoints; old41-to-new1 is not an unchanged-assay retention regression.", "Reference equality does not establish full sample phasing; Acer opalus is a related-species proxy for the Acer sample.", "These finite local gates do not establish held-out or all-panel biological performance, final macOS results, laptop performance, collaborator acceptance or release authorization."]}
    with (BASE / "independent_results_review.json").open("x") as stream:
        json.dump(review, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": "pass", "gates": 11, "plant_jobs": 8, "receipts": len(receipts)}))


if __name__ == "__main__":
    main()
