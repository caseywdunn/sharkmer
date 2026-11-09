#!/usr/bin/env python3
"""Bind full sequence-set deltas to the previously verified runtime receipts."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    path = HERE / "saved_analysis.json"
    raw = path.read_bytes()
    analysis = json.loads(raw)
    sources = {str(path): {"sha256": sha(raw), "size_bytes": len(raw)}}
    rows, prefix_bindings = [], []
    for accession in analysis["samples"]:
        runs = [r for r in analysis["runs"] if r["accession"] == accession]
        hashes = {r["input_prefix_sha256"] for r in runs}
        if len(hashes) != 1:
            raise ValueError("Historical input-prefix digests differ")
        prefix_bindings.append({"accession": accession, "raw_receipt_input_sha256": next(iter(hashes)),
                                "invocations_checked": len(runs), "scope": "Archived pre-run prefix validation; this supplement does not reread the four FASTQ files."})
    for phase in ("release_to_inclusive", "inclusive_to_retention"):
        for accession, taxon in analysis["samples"].items():
            products = {}
            for role in ("baseline", "candidate"):
                row = next(r for r in analysis["runs"] if r["phase"] == phase and r["accession"] == accession and r["role"] == role and r["pair"] == 1)
                raw = Path(row["receipt"]).read_bytes()
                seen = {"sha256": sha(raw), "size_bytes": len(raw)}
                if seen != analysis["source_files"][row["receipt"]]:
                    raise ValueError("Verified receipt changed")
                sources[row["receipt"]] = seen
                records = json.loads(raw)["outputs"]["records"]
                products[role] = {(r["gene"], r["sha256"]): {k: r[k] for k in ("gene", "sha256", "length", "assembly_status")} for r in records}
                if len(products[role]) != len(records):
                    raise ValueError("Duplicate literal identities")
            before, after = products["baseline"], products["candidate"]
            retained, lost, gained = before.keys() & after.keys(), before.keys() - after.keys(), after.keys() - before.keys()
            rows.append({"phase": phase, "accession": accession, "taxon": taxon,
                         "counts": {"baseline": len(before), "candidate": len(after), "retained": len(retained), "lost": len(lost), "gained": len(gained)},
                         "retained": [{**after[k], "previous_assembly_status": before[k]["assembly_status"]} for k in sorted(retained)],
                         "lost": [before[k] for k in sorted(lost)], "gained": [after[k] for k in sorted(gained)]})
    result = {"status": "pass", "rows": rows, "prefix_bindings": prefix_bindings, "source_files": sources,
              "producer_sha256": sha(Path(__file__).read_bytes()), "identity_rule": "Same target and full literal sequence SHA-256, independent of assembly-status reclassification."}
    with (HERE / "saved_identity_deltas.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps([{k: r[k] for k in ("phase", "accession", "counts")} for r in rows]))


if __name__ == "__main__":
    main()
