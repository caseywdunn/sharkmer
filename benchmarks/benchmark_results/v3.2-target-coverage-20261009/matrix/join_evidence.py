#!/usr/bin/env python3
"""Join available, identity-bound evidence into the frozen full target matrix.

Additional evidence JSON files contain sequence_evidence and/or target_evidence
arrays. Sequence keys are accession, gene and sha256. Each may provide
read_evidence_status/links/details or reference_evidence_status/links/details.
Evidence is appended with its source; older assessments are never discarded.
"""

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from build_matrix import HERE, RESULTS, ROOT, Sources, identity, render_tsv, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = Sources()

    def read_json(path):
        return json.loads(sources.read(path))

    base = read_json(HERE / "validation.json")
    if base["status"] != "pass":
        raise ValueError("Base matrix not validated")
    for name, receipt in base["outputs"].items():
        sources.read(HERE / name, receipt)
    for path, receipt in base["source_files"].items():
        sources.read(path, receipt)
    rows = read_json(HERE / "target_matrix.json")
    current = read_json(HERE / "sequence_inventory.json")
    lookup = {identity(r): r for r in current}
    target_lookup = {(r["accession"], r["gene"]): r for r in rows}
    evidence = {key: {"read": [], "reference": []} for key in lookup}
    target_notes = defaultdict(lambda: {"read": [], "reference": []})

    def add(key, kind, status, links, details, source):
        if key not in lookup:
            raise ValueError(f"Evidence sequence not in current frozen inventory: {key}")
        evidence[key][kind].append({"status": status, "links": links,
                                    "details": details, "evidence_source": source,
                                    "whole_product_truth": "not_established"})

    def archive_member(bundle, member):
        directory = RESULTS / bundle
        expected = sources.read(directory / "SHA256SUMS").decode().split()[0]
        raw_archive = sources.read(directory / "evidence.tar.gz")
        if sha(raw_archive) != expected:
            raise ValueError("Evidence archive checksum differs")
        inventory = read_json(directory / "ARCHIVE_CONTENTS.json")
        members = inventory["members"] if isinstance(inventory, dict) else inventory
        expected_member = next(r for r in members if r["path"] == member)
        return json.loads(sources.member(directory / "evidence.tar.gz", member, expected_member))

    # Recover exact existing reference classifications, including the 51
    # supported identities. These criteria assess gene relationship only.
    reference_path = "v3.2-public-replacements-20260913"
    member = "reclassification/saved_output_reclassification.json"
    reference = archive_member(reference_path, member)
    historical_reference_count = historical_supported_count = 0
    for product in reference["results"]["products"]:
        if product["invocation"]["version"] != "baseline":
            continue
        command, match = product["invocation"], product["new_reference_match"]
        key = command["input"], command["panel"] + "_" + product["gene"], product["sha256"]
        details = {k: match.get(k) for k in ("target_support", "taxon_support", "matched_accession", "matched_gene",
                   "matched_taxon", "pct_identity", "query_coverage_pct", "alignment_structure_status",
                   "sequence_relationship", "primer_region_support", "haplotype_truth", "thresholds")}
        add(key, "reference", match["status"], [str(RESULTS / reference_path / "README.md"),
            str(ROOT / "dev_docs/PUBLIC_REFERENCE_REPLACEMENTS.md")], details, reference_path + "::" + member)
        historical_reference_count += 1
        historical_supported_count += match["target_support"] == "supported"
    if historical_reference_count != 86 or historical_supported_count != 51:
        raise ValueError("Historical reference classification scope changed")

    gain_path = RESULTS / "v3.2-sequence-first-retention-20261001/blast-gains-v3/per_sequence.json"
    gain_assessments = read_json(gain_path)
    for product in gain_assessments:
        status = ("expected_target_hsp_present_gene_validation_not_established" if product["best_expected_target_hsp"]
                  else "other_public_sequence_hsp_only" if product["best_public_hsp"]
                  else "no_reported_hsp_in_frozen_search")
        add(identity(product), "reference", status, [str(gain_path)],
            {k: product[k] for k in ("query_id", "best_expected_target_hsp", "best_panel_hsp", "best_public_hsp",
                                     "best_same_taxon_hsp", "exact_full_query_same_taxon_sources")}, str(gain_path))
    if len(gain_assessments) != 200:
        raise ValueError("Saved gain reference assessment scope changed")

    # Reuse the historical exact read audit, binding each tested candidate by
    # sequence SHA. Marginal tiles and complete structural events stay distinct.
    read_bundle = "v3.2-read-backed-audit-20260913"
    historical_read_candidates = 0
    for accession in ("SRR27962769", "SRR31887760", "SRR1057608"):
        old_manifest = archive_member(read_bundle, f"scan_specs/{accession}.scanner_manifest.json")
        old_raw = archive_member(read_bundle, f"execution_replication_outputs/{accession}.exact.json")
        for candidate in old_manifest["candidates"]:
            possible = [key for key in lookup if key[0] == accession and key[2] == candidate["sequence_sha256"]]
            if len(possible) != 1:
                raise ValueError("Historical read candidate identity ambiguous or absent")
            key = possible[0]
            event_rows = [e for e in old_raw["events"] if e["candidate_id"] == candidate["id"]]
            marginal = [e for e in event_rows if e["role"] == "marginal_diagnostic"]
            structural = [e for e in event_rows if e["role"] != "marginal_diagnostic"]
            observed = sum(e["records"] > 0 for e in marginal)
            status = "historical_local_strings_observed" if observed else "historical_local_strings_not_observed"
            if any(e["local_bridge_status"] == "corroborated_local_bridge" and not e["longer_than_every_read"] for e in structural):
                status += ";complete_local_bridge_corroborated"
            add(key, "read", status, [str(RESULTS / read_bundle / "README.md"), str(ROOT / "dev_docs/READ_BACKED_AUDIT_RESULTS.md")],
                {"candidate_id": candidate["id"], "marginal_windows": len(marginal), "marginal_windows_observed": observed,
                 "structural_events": [{"id": e["id"], "records": e["records"],
                                        "callability": "longer_than_every_read" if e["longer_than_every_read"] else e["callability"],
                                        "window_length": e["window_end"] - e["window_start"],
                                        "local_bridge_status": e["local_bridge_status"]} for e in structural],
                 "scope": "Exact Q20 single-read windows from original first-million R1; no full-product or repeat-copy validation",
                 "input_subset": old_raw["input_subset"]}, read_bundle)
            historical_read_candidates += 1

    dros_path = RESULTS / "v3.2-dros12s-read-pilot-20261009/summary.json"
    dros = read_json(dros_path)
    k19_path = RESULTS / "v3.2-dros12s-read-pilot-20261009/k19-check/run-v1/results/report.json"
    k19 = read_json(k19_path)
    short_by_sha = {r["sequence_sha256"]: r for r in k19["candidates"] if r["role"] == "emitted_12S"}
    for product in dros["per_candidate"]:
        short = short_by_sha[product["sha256"]]
        if short["absent_19mers"] or short["minimum_raw_occurrences"] < 2:
            raise ValueError("Drosophila short-word evidence differs")
        add(("SRR31887760", "insecta_12S", product["sha256"]), "read",
            "all_canonical_19mers_observed_longer_assayed_strings_unobserved", [str(dros_path), str(k19_path)],
            {"bounded_local_assay": product, "raw_19mer_assay": short,
             "scope": "Independent short-word occurrence cannot establish their linkage; absence of longer exact matches is not biological absence"}, str(dros_path))

    its_path = RESULTS / "v3.2-its-linkage-20261009/run-v1/report.json"
    its = read_json(its_path)
    for candidate in its["candidates"]:
        key = "SRR27962769", "insecta_ITS_2", candidate["sequence_sha256"]
        add(key, "read", "saved_ledger_" + candidate["disposition"], [str(its_path)],
            {k: candidate[k] for k in ("candidate_id", "role", "length", "blocks", "segment_read_ordinals", "full_candidate_phase")}
            | {"scope": "Analysis of previously recruited 25 ITS reads; an unobserved combination is unobserved in that ledger, not a full-prefix negative"}, str(its_path))

    # Optional completed new assays append evidence at exact identity level.
    for overlay_path in args.evidence:
        overlay = read_json(overlay_path)
        for item in overlay.get("sequence_evidence", []):
            key = identity(item)
            for kind in ("read", "reference"):
                field = kind + "_evidence_status"
                if field in item:
                    links = item.get(kind + "_evidence_links", [str(overlay_path)])
                    if not isinstance(links, list) or not isinstance(item[field], str):
                        raise ValueError("Evidence status/links have wrong types")
                    add(key, kind, item[field], links, item.get(kind + "_evidence_details", {}), str(overlay_path))
        for item in overlay.get("target_evidence", []):
            key = item["accession"], item["gene"]
            if key not in target_lookup:
                raise ValueError("Evidence target outside frozen matrix")
            for kind in ("read", "reference"):
                if kind + "_evidence_status" in item:
                    target_notes[key][kind].append(item)

    # Preserve the target's actual assembly outcome separately from CLI success.
    target_results = {}
    for run in base["current_runs"]:
        receipt = read_json(run["receipt"])
        outcome = receipt["outputs"]["pcr_results"]
        if {x["gene_name"] for x in outcome} != {r["gene"] for r in rows if r["accession"] == run["accession"]}:
            raise ValueError("Raw target outcomes do not cover the complete panel")
        for result in outcome:
            target_results[(run["accession"], result["gene_name"])] = result
    for record in current:
        key = identity(record)
        for kind in ("read", "reference"):
            records = evidence[key][kind]
            record[kind + "_evidence_status"] = ";".join(sorted({r["status"] for r in records})) or "unassessed"
            record[kind + "_evidence_links"] = sorted({link for r in records for link in r["links"]})
            record[kind + "_evidence_details"] = records
    for row in rows:
        key = row["accession"], row["gene"]
        result = target_results[key]
        if (result["n_products"], result["n_standard_products"], result["n_uncertain_products"]) != (row["current_total"], row["current_standard"], row["current_uncertain"]):
            raise ValueError("Matrix product counts differ from actual target outcome")
        row["target_assembly_status"] = result["status"]
        row["target_failure_reason"] = result.get("failure_reason")
        diagnostics = result.get("threshold_diagnostics", [])
        row["search_limit_flags"] = sorted({flag for step in diagnostics for flag in
            ("dfs_limit_reached", "node_budget_reached", "path_quota_reached", "max_length_reached") if step.get(flag)})
        row["target_diagnostic_receipt"] = next(r["receipt"] for r in base["current_runs"] if r["accession"] == row["accession"])
        selected = [r for r in current if (r["accession"], r["gene"]) == key]
        for kind in ("read", "reference"):
            assessed = [r for r in selected if evidence[identity(r)][kind]]
            statuses = Counter(r[kind + "_evidence_status"] for r in assessed)
            row[kind + "_assessed_candidates"] = len(assessed)
            row[kind + "_unassessed_candidates"] = len(selected) - len(assessed)
            row[kind + "_evidence_status"] = ("not_applicable_no_emitted_product" if not selected else
                "unassessed" if not assessed else "all_emitted_candidates_have_scoped_assessments" if len(assessed) == len(selected)
                else "subset_of_emitted_candidates_has_scoped_assessments")
            row[kind + "_status_counts"] = dict(sorted(statuses.items()))
            row[kind + "_evidence_links"] = sorted({link for r in assessed for link in r[kind + "_evidence_links"]})
            row[kind + "_target_notes"] = target_notes[key][kind]
        row["historical_gene_supported_products"] = sum(any(e["details"].get("target_support") == "supported" for e in evidence[identity(r)]["reference"]) for r in selected)
    if sum(r["historical_gene_supported_products"] for r in rows) != 51:
        raise ValueError("Historical 51 supported identities not preserved")
    for before, after in zip(read_json(HERE / "target_matrix.json"), rows):
        for key in before:
            if key not in ("read_evidence_status", "read_evidence_links", "reference_evidence_status", "reference_evidence_links") and before[key] != after[key]:
                raise ValueError("Evidence join changed frozen sequence accounting")
    novel = [r for r in current if r["newly_represented_target"]]
    payloads = {}
    for name, values in (("target_matrix", rows), ("sequence_inventory", current), ("novel_target_inventory", novel)):
        payloads[name + ".json"] = (json.dumps(values, indent=2, sort_keys=True) + "\n").encode()
        payloads[name + ".tsv"] = render_tsv(values)
    summary = {"status": "pass", "base_counts": base["checks"], "historical_reference_identities_assessed": historical_reference_count,
               "historical_gene_supported_identities": historical_supported_count, "saved_gain_reference_identities_assessed": len(gain_assessments),
               "historical_read_candidates_assessed": historical_read_candidates,
               "read_assessed_candidates": sum(r["read_assessed_candidates"] for r in rows),
               "reference_assessed_candidates": sum(r["reference_assessed_candidates"] for r in rows),
               "actual_target_status_counts": dict(Counter(r["target_assembly_status"] for r in rows)),
               "source_files": sources.files, "archive_members": sources.members,
               "producer_sha256": sha(Path(__file__).read_bytes()),
               "outputs": {n: {"sha256": sha(v), "size_bytes": len(v)} for n, v in payloads.items()},
               "limitations": ["The complete matrix reports emitted sequence availability and scoped evidence, not validated target recovery.",
                   "No emitted product is not a biological absence observation. CLI success and target assembly failure are separate columns.",
                   "Reference assessments concern gene relationship or local HSPs; primer endpoints, read linkage and whole-product truth remain separate.",
                   "Old and new assays are appended by full sequence identity; aggregate target support never implies every current sequence was tested.",
                   "Whole-product validation remains not_established for every emitted record."]}
    payloads["validation.json"] = (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode()
    args.output.mkdir(parents=True, exist_ok=True)
    for name in payloads:
        if (args.output / name).exists():
            raise ValueError("Refusing evidence output overwrite")
    for name, value in payloads.items():
        (args.output / name).write_bytes(value)
    print(json.dumps({k: summary[k] for k in ("status", "historical_gene_supported_identities", "read_assessed_candidates", "reference_assessed_candidates", "actual_target_status_counts")}))


if __name__ == "__main__":
    main()
