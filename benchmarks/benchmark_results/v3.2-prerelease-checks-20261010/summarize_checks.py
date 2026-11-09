#!/usr/bin/env python3
"""Summarize completed, frozen prerelease receipts without rerunning any gate."""
import csv
import hashlib
import io
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent


def fingerprint(path):
    raw = Path(path).read_bytes()
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def main():
    sources = {}

    def read(path, expected=None):
        path = Path(path)
        observed = fingerprint(path)
        if expected and observed != {key: expected[key] for key in observed}:
            raise ValueError("Changed receipt or output: " + str(path))
        sources[str(path)] = observed
        return path.read_text()

    run = HERE / "results-v1"
    result = json.loads(read(run / "result.json"))
    assert result["status"] == "pass" and result["source_files_unchanged"] and result["bindings_unchanged"] and not result["not_run_gates"]
    plan = json.loads(read(result["plan"]["path"], result["plan"]))
    before = json.loads(read(run / "sources.before.json"))
    after = json.loads(read(run / "sources.after.json"))
    assert before == after == plan["source_files"]
    for row in [*after.values(), *plan["bindings"], result["frozen_binary"]]:
        observed = fingerprint(row["path"])
        assert observed == {key: row[key] for key in observed}, row["path"]
        sources[row["path"]] = observed
    assert [row["name"] for row in result["gates"]] == plan["gate_order"]
    gate_rows = []
    rust = {}
    for row in result["gates"]:
        directory = run / "gates" / row["name"]
        assert json.loads(read(directory / "receipt.json")) == row
        assert row["status"] == "pass" and row["returncode"] == row["metrics"]["exit_code"] == 0
        for file in row["files"]:
            read(file["path"], file)
        constraints = json.loads((directory / "constraints.json").read_text())
        assert constraints["status"] == "pass" and constraints["affinity"] == result["cpus"]
        assert constraints["rlimit_as_soft_hard"] == [plan["resources"]["child_address_space_bytes"]] * 2
        gate_rows.append({"gate": row["name"], "status": row["status"], "wall_seconds": row["metrics"]["wall_seconds"],
                          "user_seconds": row["metrics"]["user_seconds"], "system_seconds": row["metrics"]["system_seconds"],
                          "peak_rss_mib": row["metrics"]["max_rss_kib"] / 1024,
                          "timeout_seconds": row["timeout_seconds"], "receipt": str(directory / "receipt.json")})
        if row["name"] in ("test_default", "test_fxhash"):
            summaries = re.findall(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out", (directory / "stdout.log").read_text())
            assert len(summaries) == 3 and all(all(int(n) == 0 for n in counts[1:]) for counts in summaries)
            rust[row["name"]] = {"unit": int(summaries[0][0]), "spcr_integration": int(summaries[1][0]),
                                "uncertain_integration": int(summaries[2][0]), "passed": sum(int(row[0]) for row in summaries),
                                "failed": 0, "ignored": 0}
    python_log = (run / "gates/python_all/stderr.log").read_text()
    python_count = int(re.search(r"Ran (\d+) tests", python_log)[1])
    assert python_log.rstrip().endswith("OK") and "skipped=" not in python_log
    assert re.search(r"test_real_database_build_and_query_binds_logical_gene_metadata.*\.\.\. ok", python_log)
    reference = json.loads(read(run / "references/verification.json"))
    panels = json.loads(read(run / "panels/verification.json"))
    fixtures = json.loads(read(run / "fixtures/verification.json"))
    plant = json.loads(read(run / "plant/summary.json"))
    assert all(item["status"] == "pass" for item in (reference, panels, fixtures, plant))
    assert len(plant["jobs"]) == 8 and len(fixtures["runs"]) == 2
    gates = io.StringIO(newline="")
    writer = csv.DictWriter(gates, fieldnames=list(gate_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(gate_rows)
    report = {"status": "pass", "slurm_job_id": result["slurm_job_id"], "machine": result["machine"], "cpus": result["cpus"],
        "resources": plan["resources"], "tool_versions": result["tool_versions"], "dependency_versions": plan["dependency_versions"],
        "frozen_default_binary": result["frozen_binary"], "source_file_count": len(before), "binding_count": len(plan["bindings"]),
        "source_files_and_bindings_unchanged": True, "gates": gate_rows, "rust_tests": rust,
        "python": {"passed": python_count, "skipped": 0, "real_blast_test": "passed"},
        "references": {key: reference["audit"][key] for key in ("panel_count", "total_references", "verified_references", "excluded_references", "failure_count")},
        "panels": {key: panels[key] for key in ("jsonschema_version", "schema_checks", "builtin_exports", "loader_dry_runs")},
        "fixture_oracles": fixtures, "plant": {"status": plant["status"], "invocations": len(plant["jobs"]), "summary": str(run / "plant/summary.json")},
        "source_files": sources, "producer": fingerprint(Path(__file__)),
        "limits": ["Local frozen-build quality and assay checks. No release or performance acceptance is implied.",
                   "Plant runs are single-CPU functionality checks with a changed assay, not a runtime comparison.",
                   "The k31 regression oracle checks exact standard products; every emitted record is preserved in fixture evidence.",
                   "GNU peak RSS, allocator memory, process address-space bounds, and Slurm memory allocation are distinct quantities.",
                   "The Python suite actually used local BLAST tools and has no skipped tests."]}
    output = HERE / "summary-v1"
    output.mkdir(exist_ok=False)
    (output / "gates.tsv").write_text(gates.getvalue())
    (output / "summary.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "pass", "rust": rust, "python_passed": python_count, "source_files": len(before), "bindings": len(plan["bindings"])}))


if __name__ == "__main__":
    main()
