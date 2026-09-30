#!/usr/bin/env python3
"""Derive literal historical product identities from authenticated archived FASTAs.

No biological classification is copied or inferred. Scope labels are the frozen
high-copy/deferred target grouping. All three historical replicates must have
identical ordered sequences per target before a single expectation is emitted.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import tarfile

from restore_inputs import fingerprint, write_json


HERE = Path(__file__).resolve().parent


def parse_fasta(raw):
    records = []
    for line in raw.decode("ascii").splitlines():
        if line.startswith(">"):
            if not line[1:].strip():
                raise ValueError("Empty archived FASTA header")
            records.append({"id": line[1:].split()[0], "description": line[1:], "sequence": ""})
        elif line:
            if not records or not re.fullmatch("[ACGTNacgtn]+", line):
                raise ValueError("Unexpected archived FASTA sequence")
            records[-1]["sequence"] += line
    for record in records:
        sequence = record["sequence"]
        if not sequence:
            raise ValueError("Empty archived FASTA record")
        reverse = sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
        record.update(length=len(sequence),
                      sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                      canonical_sha256=hashlib.sha256(min(sequence.upper(), reverse).encode()).hexdigest())
    return records


def extract(inventory_path, output):
    inventory = json.loads(inventory_path.read_text())
    sources = {key: inventory["historical_sources"][key] for key in
               ("final_provenance", "final_archive", "final_inventory")}
    sources["inventory"] = fingerprint(inventory_path)
    sources["extractor"] = fingerprint(__file__)
    for row in sources.values():
        if fingerprint(row["path"]) != row:
            raise ValueError("Historical source differs from inventory")
    provenance = json.loads(Path(sources["final_provenance"]["path"]).read_text())
    content_rows = json.loads(Path(sources["final_inventory"]["path"]).read_text())
    contents = {row["path"]: row for row in content_rows}
    if len(contents) != len(content_rows):
        raise ValueError("Duplicate archive inventory member")
    scope = {panel: {row["gene"]: row["scope"] for row in rows}
             for panel, rows in provenance["protocol"]["target_metadata"].items()}
    invocations = []
    with tarfile.open(sources["final_archive"]["path"], "r:gz") as archive:
        members = archive.getmembers()
        if len({member.name for member in members}) != len(members):
            raise ValueError("Duplicate archived member path")
        for command in provenance["schedule"]:
            accession, panel = command["input"], command["panel"]
            prefix = f"final-execution/attempts/{command['invocation_id']}/"
            fasta_members = sorted((member for member in members
                                    if member.name.startswith(prefix) and member.name.endswith(".fasta")),
                                   key=lambda member: member.name)
            by_gene = {gene: [] for gene in scope[panel]}
            authenticated = []
            seen_genes = set()
            for member in fasta_members:
                if not member.isfile() or member.size > 16 * 1024 * 1024:
                    raise ValueError("Invalid or excessive archived FASTA member")
                raw = archive.extractfile(member).read()
                expected = contents[member.name]
                if len(raw) != expected["bytes"] or hashlib.sha256(raw).hexdigest() != expected["sha256"]:
                    raise ValueError("Archived FASTA member differs from manifest")
                name_prefix = f"{panel}_{accession}_{command['depth']}_{panel}_"
                filename = Path(member.name).name
                if not filename.startswith(name_prefix):
                    raise ValueError("Unexpected FASTA filename namespace")
                gene = filename[len(name_prefix):-len(".fasta")]
                if gene not in by_gene or gene in seen_genes:
                    raise ValueError("Unexpected or duplicated archived gene output")
                seen_genes.add(gene)
                by_gene[gene] = parse_fasta(raw)
                for record in by_gene[gene]:
                    record["archived_member"] = member.name
                    record["scope"] = scope[panel][gene]
                authenticated.append(expected)
            invocations.append({"invocation_id": command["invocation_id"], "accession": accession,
                                "panel": panel, "role": command["version"], "pair_index": command["pair_index"],
                                "by_gene": by_gene, "authenticated_fasta_members": authenticated})
    if len(invocations) != 60:
        raise ValueError("Expected exactly 60 historical invocations")
    expectations = []
    totals = {role: {"all_products": 0, "high_copy_products": 0, "deferred_products": 0}
              for role in ("baseline", "candidate")}
    for sample in inventory["samples"]:
        accession, panel = sample["accession"], sample["panel"]
        roles = {}
        for role in ("baseline", "candidate"):
            replicates = sorted((row for row in invocations
                                 if row["accession"] == accession and row["role"] == role),
                                key=lambda row: row["pair_index"])
            if [row["pair_index"] for row in replicates] != [1, 2, 3]:
                raise ValueError("Missing or repeated historical replicate")
            literal = [{gene: [record["sequence"] for record in records]
                        for gene, records in row["by_gene"].items()} for row in replicates]
            if not literal[0] == literal[1] == literal[2]:
                raise ValueError(f"Historical sequence order differs across replicates: {accession}/{role}")
            by_gene = replicates[0]["by_gene"]
            records = [record for rows in by_gene.values() for record in rows]
            count = len(records)
            high_copy = sum(record["scope"] == "high_copy_candidate" for record in records)
            roles[role] = {"all_products": count, "high_copy_products": high_copy,
                           "deferred_products": count - high_copy, "by_gene": by_gene,
                           "replicate_sequence_parity": True,
                           "source_invocations": [row["invocation_id"] for row in replicates]}
            for key in totals[role]:
                totals[role][key] += roles[role][key]
        expectations.append({"accession": accession, "panel": panel, "roles": roles})
    if totals != {"baseline": {"all_products": 86, "high_copy_products": 71, "deferred_products": 15},
                  "candidate": {"all_products": 64, "high_copy_products": 64, "deferred_products": 0}}:
        raise ValueError(f"Historical totals disagree: {totals}")
    after = {key: fingerprint(row["path"]) for key, row in sources.items()}
    if after != sources:
        raise ValueError("Historical sources changed during extraction")
    evidence = [{key: value for key, value in row.items() if key != "by_gene"} for row in invocations]
    result = {"schema_version": 1, "status": "pass", "sources_before": sources, "sources_after": after,
              "scope": "Literal oriented and reverse-complement-canonical identity expectations; no biological validity claims",
              "baseline_commit": "5a664680c91ad59f32b8c2a847b8fef37f34a0ae",
              "candidate_role": "Historical September 11 pre-inclusive candidate; not the new integrated candidate",
              "totals": totals, "sample_expectations": expectations, "invocation_evidence": evidence}
    write_json(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "inventory.json")
    parser.add_argument("--output", type=Path, default=HERE / "historical_identities.json")
    args = parser.parse_args()
    result = extract(args.inventory, args.output)
    print(json.dumps({"status": result["status"], "totals": result["totals"],
                      "output": fingerprint(args.output)}, sort_keys=True))
