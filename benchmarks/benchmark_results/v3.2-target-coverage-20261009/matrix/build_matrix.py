#!/usr/bin/env python3
"""Reconstruct the full frozen sample-target comparison without running assembly."""

import argparse
import csv
import hashlib
import io
import json
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RESULTS = ROOT / "benchmarks/benchmark_results"
PRIOR = RESULTS / "v3.2-inclusive-release-20260929/inputs"
RETENTION = RESULTS / "v3.2-sequence-first-retention-20261001"
RUNS = ROOT / "tmp/sequence-retention-20261001/full-v3"
EXPECTED_NEW_TARGETS = {
    ("SRR25378184", "angiospermae_psbA-trnH"): 49,
    ("SRR25378184", "angiospermae_trnV-atpE"): 41,
    ("SRR25378184", "angiospermae_atpB-rbcL"): 1,
    ("ERR14009273", "angiospermae_trnL-F"): 1,
    ("SRR1057608", "insecta_ND4"): 44,
    ("SRR27962769", "insecta_18S_1"): 18,
}
SUPPORT_DEFAULTS = {
    "read_evidence_status": "unassessed", "read_evidence_links": [],
    "reference_evidence_status": "unassessed", "reference_evidence_links": [],
    "whole_product_validation": "not_established",
}


def sha(value):
    return hashlib.sha256(value).hexdigest()


def parse_fasta(raw):
    result = []
    for line in raw.decode("ascii").splitlines():
        if line.startswith(">"):
            result.append({"header": line[1:], "sequence": ""})
        elif line:
            if not result or set(line.upper()) - set("ACGTN"):
                raise ValueError("Malformed FASTA")
            result[-1]["sequence"] += line
    if any(not r["sequence"] for r in result):
        raise ValueError("Empty sequence")
    return result


def identity(row):
    return row["accession"], row["gene"], row["sha256"]


def render_tsv(rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows({k: json.dumps(v, sort_keys=True, separators=(",", ":")) if isinstance(v, (list, dict)) else v
                     for k, v in row.items()} for row in rows)
    return stream.getvalue().encode()


class Sources:
    def __init__(self):
        self.files = {}
        self.members = {}
        self.archives = {}

    def read(self, path, expected=None):
        path = Path(path)
        value = path.read_bytes()
        receipt = {"sha256": sha(value), "size_bytes": len(value)}
        if expected and (receipt["sha256"] != expected["sha256"] or
                         receipt["size_bytes"] != expected.get("size_bytes", expected.get("bytes"))):
            raise ValueError(f"Source fingerprint differs: {path}")
        self.files[str(path)] = receipt
        return value

    def member(self, archive_path, name, expected=None):
        archive_path = Path(archive_path)
        if archive_path not in self.archives:
            self.archives[archive_path] = tarfile.open(fileobj=io.BytesIO(self.read(archive_path)), mode="r:gz")
        value = self.archives[archive_path].extractfile(name).read()
        if expected and (sha(value) != expected["sha256"] or len(value) != expected.get("size_bytes", expected.get("bytes"))):
            raise ValueError(f"Archive member differs: {name}")
        self.members[f"{archive_path}::{name}"] = {"sha256": sha(value), "size_bytes": len(value)}
        return value

    def current(self, path, expected=None):
        path = Path(path)
        if path.is_file():
            return self.read(path, expected)
        receipt = json.loads(self.read(RETENTION / "archive_receipt.json"))
        self.read(RETENTION / "evidence.tar.gz", receipt)
        name = "scratch/full-v3/" + str(path.relative_to(RUNS))
        return self.member(RETENTION / "evidence.tar.gz", name, expected)


def apply_evidence(rows, overlays, key_fields):
    lookup = {tuple(row[k] for k in key_fields): row for row in rows}
    seen = set()
    for evidence in overlays:
        key = tuple(evidence[k] for k in key_fields)
        if key not in lookup or key in seen:
            raise ValueError("Unknown or duplicate evidence identity")
        seen.add(key)
        extras = set(evidence) - set(key_fields) - set(SUPPORT_DEFAULTS)
        if extras:
            raise ValueError(f"Unsupported evidence overlay fields: {extras}")
        for field in SUPPORT_DEFAULTS:
            if field in evidence:
                lookup[key][field] = evidence[field]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--evidence", type=Path, help="Optional JSON with target_evidence and sequence_evidence arrays; identity keys accession,gene[,sha256]")
    args = parser.parse_args()
    sources = Sources()
    inventory = json.loads(sources.read(PRIOR / "inventory.json"))
    historical = json.loads(sources.read(PRIOR / "historical_identities.json"))
    if historical["status"] != "pass" or historical["sources_before"] != historical["sources_after"]:
        raise ValueError("Historical provenance failed")
    for item in historical["sources_before"].values():
        sources.read(item["path"], item)
    samples = {r["accession"]: r for r in inventory["samples"]}
    panels = {}
    for entry in inventory["panels"]:
        panel = yaml.safe_load(sources.read(entry["assembly_panel"]["path"], entry["assembly_panel"]))
        targets = [str(p["gene"]) + ("-" + str(p["region"]) if "region" in p else "")
                   + ("_" + str(p["index"]) if "index" in p else "") for p in panel["primers"]]
        if len(targets) != entry["target_count"] or len(set(targets)) != len(targets):
            raise ValueError("Panel target scope differs")
        panels[entry["name"]] = targets
    query_map = defaultdict(list)
    blast_path = RETENTION / "blast-gains-v3/per_sequence.json"
    for row in json.loads(sources.read(blast_path)):
        query_map[identity(row)].append({"query_id": row["query_id"], "namespace": "sequence_retention_blast_gains_v3", "source": str(blast_path)})
    old_members = {}
    for invocation in historical["invocation_evidence"]:
        if invocation["role"] == "baseline" and invocation["pair_index"] == 1:
            old_members.update({f["path"]: f for f in invocation["authenticated_fasta_members"]})
    old_parsed = {}
    old_archive = historical["sources_before"]["final_archive"]["path"]
    for name, fingerprint in old_members.items():
        old_parsed[name] = parse_fasta(sources.member(old_archive, name, fingerprint))
    old, current, run_sources = [], [], []
    for expectation in historical["sample_expectations"]:
        accession, panel = expectation["accession"], expectation["panel"]
        sample = samples[accession]
        baseline = expectation["roles"]["baseline"]
        if not baseline["replicate_sequence_parity"]:
            raise ValueError("Historical replicate identity unverified")
        for target, records in baseline["by_gene"].items():
            for record in records:
                if {"header": record["description"], "sequence": record["sequence"]} not in old_parsed[record["archived_member"]]:
                    raise ValueError("Historical sequence differs from authenticated FASTA")
                value = {"accession": accession, "taxon": sample["taxon"], "panel": panel,
                         "target": target, "gene": panel + "_" + target, "length": record["length"],
                         "sha256": record["sequence_sha256"], "sequence": record["sequence"],
                         "header": record["description"], "assembly_status": "historical_emitted",
                         "source": str(old_archive) + "::" + record["archived_member"],
                         "scope": record["scope"]}
                if sha(value["sequence"].encode()) != value["sha256"] or len(value["sequence"]) != value["length"]:
                    raise ValueError("Historical sequence digest differs")
                old.append(value)
    old_keys = {identity(r) for r in old}
    if len(old_keys) != len(old):
        raise ValueError("Duplicate historical identity")
    for index, sample in enumerate(inventory["samples"]):
        accession, panel = sample["accession"], sample["panel"]
        stem = f"{index:03d}_{panel}_{accession}_1000000_pair1_candidate"
        path = RUNS / (stem + ".json")
        job = json.loads(sources.current(path))
        if job["status"] != "pass" or job["execution"]["returncode"] != 0 or job["invocation"]["input"] != accession:
            raise ValueError("Current candidate run failed or misbound")
        if job["prefix_before"] != sample["expected_prefix"]:
            raise ValueError("Current and historical consumed-prefix identities differ")
        parsed = []
        for filename, fingerprint in job["outputs"]["files"].items():
            if filename.endswith(".fasta"):
                gene = filename[len(job["invocation"]["sample_prefix"]) + 1:-6]
                for record in parse_fasta(sources.current(fingerprint["path"], fingerprint)):
                    parsed.append((gene, record["header"], record["sequence"]))
        expected_records = [(r["gene"], r["header"], r["sequence"]) for r in job["outputs"]["records"]]
        if sorted(parsed) != sorted(expected_records):
            raise ValueError("Current receipt differs from authenticated FASTAs")
        for record in job["outputs"]["records"]:
            if sha(record["sequence"].encode()) != record["sha256"] or len(record["sequence"]) != record["length"]:
                raise ValueError("Current sequence digest differs")
            key = accession, record["gene"], record["sha256"]
            if record["assembly_status"] not in ("standard", "uncertain"):
                raise ValueError("Unexpected current assembly status")
            queries = query_map[key]
            current.append({"accession": accession, "taxon": sample["taxon"], "panel": panel,
                            "target": record["gene"][len(panel) + 1:], "gene": record["gene"],
                            "length": record["length"], "sha256": record["sha256"],
                            "sequence": record["sequence"], "header": record["header"],
                            "assembly_status": record["assembly_status"], "historical_retained": key in old_keys,
                            "new_vs_release": key not in old_keys, "existing_blast_query_id": queries[0]["query_id"] if queries else None,
                            "blast_queries": queries, "scope": record["scope"], "source": str(path), **SUPPORT_DEFAULTS})
        run_sources.append({"accession": accession, "panel": panel, "receipt": str(path),
                            "status": job["status"], "prefix": job["prefix_before"],
                            "products": len(job["outputs"]["records"]), "invocation_id": stem})
    current.sort(key=identity)
    old.sort(key=identity)
    current_keys = {identity(r) for r in current}
    if len(current_keys) != len(current):
        raise ValueError("Duplicate current identity")
    rows = []
    for sample in inventory["samples"]:
        accession, panel = sample["accession"], sample["panel"]
        for target in panels[panel]:
            gene = panel + "_" + target
            before = [r for r in old if (r["accession"], r["gene"]) == (accession, gene)]
            after = [r for r in current if (r["accession"], r["gene"]) == (accession, gene)]
            before_hashes, after_hashes = {r["sha256"] for r in before}, {r["sha256"] for r in after}
            rows.append({"accession": accession, "taxon": sample["taxon"], "panel": panel, "target": target, "gene": gene,
                         "attempted": True, "run_status": "pass", "v3_1_products": len(before),
                         "current_standard": sum(r["assembly_status"] == "standard" for r in after),
                         "current_uncertain": sum(r["assembly_status"] == "uncertain" for r in after),
                         "current_total": len(after), "retained_count": len(before_hashes & after_hashes),
                         "lost_count": len(before_hashes - after_hashes), "new_count": len(after_hashes - before_hashes),
                         "v3_1_identities": sorted(before_hashes), "current_identities": sorted(after_hashes),
                         "retained_identities": sorted(before_hashes & after_hashes),
                         "lost_identities": sorted(before_hashes - after_hashes), "new_identities": sorted(after_hashes - before_hashes),
                         "newly_represented_target": not before and bool(after),
                         "zero_output_current": not after, **SUPPORT_DEFAULTS})
    valid_targets = {(r["accession"], r["gene"]) for r in rows}
    if any((r["accession"], r["gene"]) not in valid_targets for r in old + current):
        raise ValueError("Output outside frozen panel scope")
    new_targets = {(r["accession"], r["gene"]): r["current_total"] for r in rows if r["newly_represented_target"]}
    if new_targets != EXPECTED_NEW_TARGETS:
        raise ValueError(f"Newly represented target identities differ: {new_targets}")
    for record in current:
        record["newly_represented_target"] = (record["accession"], record["gene"]) in new_targets
    for record in old:
        record["retained_in_current"] = identity(record) in current_keys
    if args.evidence:
        overlay = json.loads(sources.read(args.evidence))
        apply_evidence(rows, overlay.get("target_evidence", []), ("accession", "gene"))
        apply_evidence(current, overlay.get("sequence_evidence", []), ("accession", "gene", "sha256"))
    novel = [r for r in current if r["newly_represented_target"]]
    checks = {"target_rows": len(rows), "historical_products": len(old), "current_products": len(current),
              "standard_products": sum(r["current_standard"] for r in rows), "uncertain_products": sum(r["current_uncertain"] for r in rows),
              "retained_exact_products": len(old_keys & current_keys), "lost_exact_products": len(old_keys - current_keys),
              "new_exact_products": len(current_keys - old_keys), "historical_represented_targets": sum(r["v3_1_products"] > 0 for r in rows),
              "current_represented_targets": sum(r["current_total"] > 0 for r in rows), "newly_represented_targets": len(new_targets),
              "novel_target_products": len(novel), "zero_output_current_targets": sum(r["zero_output_current"] for r in rows)}
    expected_checks = {"target_rows": 118, "historical_products": 86, "current_products": 326,
                       "standard_products": 64, "uncertain_products": 262, "retained_exact_products": 86, "lost_exact_products": 0,
                       "new_exact_products": 240, "historical_represented_targets": 68, "current_represented_targets": 74,
                       "newly_represented_targets": 6, "novel_target_products": 154, "zero_output_current_targets": 44}
    if checks != expected_checks or sum(r["retained_count"] + r["lost_count"] for r in rows) != len(old) or sum(r["retained_count"] + r["new_count"] for r in rows) != len(current):
        raise ValueError(f"Set accounting differs: {checks}")
    payloads = {}
    for name, records in (("target_matrix", rows), ("sequence_inventory", current), ("historical_inventory", old), ("novel_target_inventory", novel)):
        payloads[name + ".json"] = (json.dumps(records, indent=2, sort_keys=True) + "\n").encode()
        payloads[name + ".tsv"] = render_tsv(records)
        if name != "target_matrix":
            payloads[name + ".fasta"] = "".join(f">{r['accession']}|{r['gene']}|{r['sha256']} assembly_status={r['assembly_status']}\n{r['sequence']}\n" for r in records).encode()
    validation = {"status": "pass", "checks": checks, "newly_represented_targets": [r for r in rows if r["newly_represented_target"]],
                  "source_files": sources.files, "archive_members": sources.members, "current_runs": run_sources,
                  "producer_sha256": sha(Path(__file__).read_bytes()), "evidence_overlay": str(args.evidence) if args.evidence else None,
                  "outputs": {name: {"sha256": sha(value), "size_bytes": len(value)} for name, value in payloads.items()},
                  "limitations": ["Exact identity is keyed by accession, target and oriented full-sequence SHA-256; no biological correctness follows from retention.",
                      "Historical v3.1.0 products are observations, not references. Current standard and uncertain categories retain their original software meaning.",
                      "All 118 attempted sample-targets from complete frozen panels are represented, including 44 current zero-output targets.",
                      "Read/reference support defaults to unassessed; only explicit evidence overlays alter those fields. Whole-product truth is separate.",
                      "Existing BLAST query IDs are bound only where verified mappings exist and include their namespace; missing IDs are null.",
                      "Replicate stability is inherited from previously completed reviews; this matrix authenticates first-replicate receipts and raw FASTAs without rerunning assembly."]}
    payloads["validation.json"] = (json.dumps(validation, indent=2, sort_keys=True) + "\n").encode()
    args.output.mkdir(parents=True, exist_ok=True)
    for name in payloads:
        if (args.output / name).exists():
            raise ValueError(f"Refusing to overwrite {name}; use a fresh --output directory")
    for name, value in payloads.items():
        (args.output / name).write_bytes(value)
    print(json.dumps({"status": "pass", "checks": checks, "output": str(args.output)}))


if __name__ == "__main__":
    main()
