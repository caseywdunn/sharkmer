"""Independent existing-artifact review; no assembler or new biological analysis."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import shlex
import sys
import yaml

ROOT = Path(__file__).resolve().parents[4]
BUNDLE = Path(__file__).resolve().parents[1]
SCRATCH = ROOT / "tmp/inclusive-release-20260929/dros12s"
OBSERVED = SCRATCH / "replay-v1"
CONTROL = ROOT / "tmp/inclusive-output-20260929/attempt-v1/runs/SRR31887760/default"
SAMPLE = "uncertain_SRR31887760"


def fingerprint(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return {"path": str(path), "sha256": digest.hexdigest(), "size_bytes": path.stat().st_size}


def main():
    receipt_path = OBSERVED / "replay.receipt.json"
    receipt = json.loads(receipt_path.read_text())
    assessment_path = OBSERVED / "assessment.json"
    assessment = json.loads(assessment_path.read_text())
    assert receipt["status"] == assessment["status"] == "pass" and receipt["returncode"] == 0
    assert receipt["before"] == receipt["after"] == receipt["after_audit"]
    assert receipt["audit_inputs_before"] == receipt["audit_inputs_after"]
    artifacts = {}
    for rows in (receipt["before"], receipt["outputs"]):
        for row in rows.values():
            assert fingerprint(row["path"]) == row, row["path"]
            artifacts[row["path"]] = row
    for relative, row in receipt["source_freeze"]["files"].items():
        actual = fingerprint(SCRATCH / "source-v1" / relative)
        assert all(actual[key] == row[key] for key in ("sha256", "size_bytes"))
        artifacts[actual["path"]] = actual
    assert receipt["source_freeze"] == json.loads((SCRATCH / "source-freeze-v1.json").read_text())
    assert receipt["observed_constraints"]["command"] == receipt["argv"]
    assert receipt["observed_constraints"]["affinity"] == receipt["declared_limits"]["cpu_affinity"]
    assert receipt["observed_constraints"]["rlimit_as_soft_hard"] == [40 * 1024**3] * 2
    assert receipt["assessment"] == fingerprint(assessment_path)
    stats = yaml.safe_load((OBSERVED / f"{SAMPLE}.stats.yaml").read_text())
    control = yaml.safe_load((CONTROL / f"{SAMPLE}.stats.yaml").read_text())
    keys = ("n_reads_read", "n_bases_read", "n_subreads_ingested", "n_bases_ingested", "n_kmers",
            "kmer_length", "chunks", "input_source", "sample", "pcr_results", "product_output_schema_version", "product_output_layout")
    assert all(stats[key] == control[key] for key in keys)
    assert len(stats["pcr_results"]) == 21
    assert stats["n_reads_read"] == 1000000 and stats["n_bases_read"] == 101000000
    assert stats["stage_timings"]["read_threading_input_s"] is None
    expected_ledgers, all_rows, marker_budget = set(), {}, None
    output_counts = Counter()
    cause_keys = ("omitted_self_loop_occurrences", "cyclic_scc_occurrences", "retained_collision_occurrences")
    source_keys = ("omitted_self_loop_marker_occurrences", "cyclic_scc_marker_occurrences", "retained_collision_marker_occurrences")
    for gene in stats["pcr_results"]:
        name = gene["gene_name"]
        if gene.get("output_file"):
            raw = (OBSERVED / gene["output_file"]).read_bytes()
            assert raw == (CONTROL / gene["output_file"]).read_bytes()
            for line in raw.decode().splitlines():
                if line.startswith(">"):
                    labels = [part.split("=", 1)[1] for part in shlex.split(line) if part.startswith("assembly_status=")]
                    assert len(labels) == 1 and labels[0] in ("standard", "uncertain")
                    output_counts[labels[0]] += 1
        uncertainty = gene["uncertain_products"]
        retained = Counter()
        all_rows[name] = []
        for search in uncertainty["thresholds"]:
            if not search["path_search_evaluated"]:
                assert search["withheld_candidate_paths"] == 0
                continue
            path = OBSERVED / f"completion_ledger_{hashlib.sha256(name.encode()).hexdigest()}_{search['threshold']}.yaml"
            expected_ledgers.add(path.name)
            doc = yaml.safe_load(path.read_text())
            ledger = doc["ledger"]
            assert (doc["sample"], doc["gene"], doc["threshold"]) == (SAMPLE, name, search["threshold"])
            assert doc["complete"] and ledger["omitted_rows"] == 0
            rows = ledger["rows"]
            assert len(rows) == ledger["observed_completions"] == search["withheld_candidate_paths"]
            assert [row["ordinal"] for row in rows] == list(range(1, len(rows) + 1))
            counts = Counter(row["retention_disposition"] for row in rows)
            assert counts == ledger["disposition_counts"]
            retention = search["retention"]
            assert retention["observed_withheld_paths"] == len(rows)
            assert retention["retained_paths"] == counts["retained"]
            assert retention["observed_not_retained_paths"] == len(rows) - counts["retained"]
            for reason in ("oversize_sequence", "oversize_marker_set", "path_cap", "sequence_base_cap", "marker_cap"):
                assert retention["dropped_" + reason] == counts[reason]
            for aggregate, cause in zip(("withheld_by_self_loop_node_paths", "withheld_by_scc_node_paths", "withheld_by_collision_edge_paths"), cause_keys):
                assert search[aggregate] == sum(row[cause] > 0 for row in rows)
            for row in rows:
                all_rows[name].append(dict(row, threshold=search["threshold"]))
                if row["retention_disposition"] == "retained":
                    retained[(search["threshold"], row["sequence_sha256"], row["sequence_length"], *(row[key] for key in cause_keys))] += 1
            if name == "insecta_12S":
                limits = uncertainty["limits"]
                allocated = limits["allocated_run_share_for_gene"]["marker_occurrences"]
                retained_rows = [row for row in rows if row["retention_disposition"] == "retained"]
                used = sum(sum(row[key] for key in cause_keys) for row in retained_rows)
                assert allocated == 1561 and used == 1274 == retention["retained_marker_occurrences"]
                assert counts == {"retained": 3, "marker_cap": 23} and len({row["sequence_sha256"] for row in rows}) == 26
                assert [sum(row[key] for key in cause_keys) for row in retained_rows] == [414, 427, 433]
                assert all(used + sum(row[key] for key in cause_keys) > allocated for row in rows if row["retention_disposition"] == "marker_cap")
                marker_budget = {"allocated_run_share": allocated, "retained_marker_occurrences": used, "remaining": allocated - used,
                                 "threshold_cap": limits["threshold"]["marker_occurrences"], "gene_cap": limits["gene"]["marker_occurrences"],
                                 "retained_paths": retention["retained_paths"], "allocated_path_share": limits["allocated_run_share_for_gene"]["paths"],
                                 "retained_bases": retention["retained_sequence_bases"], "allocated_base_share": limits["allocated_run_share_for_gene"]["sequence_bases"]}
        expected_retained = Counter()
        for product in uncertainty["products"]:
            sequence = product["oriented_sequence"]
            assert hashlib.sha256(sequence.encode()).hexdigest() == product["sequence_sha256"]
            assert len(sequence) == product["sequence_length"]
            for source in product["sources"]:
                expected_retained[(source["threshold"], product["sequence_sha256"], len(sequence), *(source[key] for key in source_keys))] += 1
        assert uncertainty["excluded_standard_equivalent_path_records"] == 0
        assert retained == expected_retained
    assert {p.name for p in OBSERVED.glob("completion_ledger_*")} == expected_ledgers
    assert len(expected_ledgers) == 16 and output_counts == {"standard": 9, "uncertain": 4}
    history_path = ROOT / "benchmarks/benchmark_results/v3.2-uncertain-products-20260927/inputs/historical_losses.json"
    history = json.loads(history_path.read_text())["rows"]
    findings = []
    for historical in history:
        if historical["accession"] != "SRR31887760":
            continue
        sequence = historical["sequence"]
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        assert digest == historical["sequence_sha256"] and len(sequence) == historical["length"]
        matches = [row for row in all_rows["insecta_" + historical["gene"]] if row["sequence_sha256"] == digest and row["sequence_length"] == len(sequence)]
        assert len(matches) == 1
        row = matches[0]
        reverse = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        assert hashlib.sha256(reverse.encode()).hexdigest() == row["reverse_complement_sha256"]
        findings.append({"gene": "insecta_" + historical["gene"], "length": len(sequence), "sha256": digest,
                         "threshold": row["threshold"], "ordinal": row["ordinal"], "disposition": row["retention_disposition"],
                         "marker_occurrences": sum(row[key] for key in cause_keys)})
    missing = {row["length"]: row for row in findings if row["gene"] == "insecta_12S" and row["length"] in (487, 499)}
    assert missing[499]["ordinal"] == 13 and missing[499]["marker_occurrences"] == 439
    assert missing[487]["ordinal"] == 19 and missing[487]["marker_occurrences"] == 426
    assert all(row["disposition"] == "marker_cap" for row in missing.values())
    for p in (Path(__file__), receipt_path, assessment_path, history_path):
        row = fingerprint(p); artifacts[row["path"]] = row
    result = {"schema_version": 1, "status": "pass", "date": "2026-09-30", "purpose": "independent_scratch_completion_outcome_review",
              "reviewed_artifacts": artifacts, "normal_output_parity": True, "output_counts": dict(output_counts), "complete_threshold_ledgers": 16,
              "historical_exact_matches": findings, "marker_budget": marker_budget,
              "conclusion": "Both 499-bp and 487-bp historical identities were enumerated exactly but omitted from retained uncertainty payloads by the static run marker share allocated to 12S.",
              "limitations": ["Scratch observer changed no ordinary products, search counters or retention budgets.",
                              "The two sequences remain absent from production FASTA output; this diagnostic does not restore them.",
                              "Exact historical identity does not establish biological correctness, repeat length or independent read support.",
                              "The enumeration remains bounded by graph, DFS, length and visit constraints; observer timing is not a release benchmark."]}
    output = BUNDLE / "review/dros12s_outcome_review.json"
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True); stream.write("\n")
    print(json.dumps({"status": "pass", "review": fingerprint(output), "historical": findings, "budget": marker_budget}))


if __name__ == "__main__":
    main()
