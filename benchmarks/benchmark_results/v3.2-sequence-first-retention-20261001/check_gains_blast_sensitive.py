#!/usr/bin/env python3
"""Sensitive BLAST of gained hypotheses against a pinned public-record catalog."""

import argparse
import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.sharkmer_validate import reference_provenance
from scripts.sharkmer_validate.reference_targets import logical_gene_name


def run(command):
    subprocess.run(command, check=True)


def receipt(path):
    return {"path": str(path), "size_bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def load_base():
    path = Path(__file__).with_name("check_gains_blast.py")
    spec = importlib.util.spec_from_file_location("check_gains_blast", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def query_rows(path):
    with path.open(newline="") as stream:
        for fields in csv.reader(stream, delimiter="\t"):
            if not fields:
                continue
            if len(fields) != 14:
                raise ValueError(f"Unexpected BLAST row width: {fields}")
            (qid, sid, pident, length, qlen, qstart, qend, slen, sstart, send,
             evalue, bitscore, gaps, qcovs) = fields
            yield {"query_id": qid, "subject_id": sid,
                   "identity_pct": float(pident), "alignment_length": int(length),
                   "query_length": int(qlen), "qstart": int(qstart), "qend": int(qend),
                   "subject_length": int(slen), "sstart": int(sstart), "send": int(send),
                   "evalue": float(evalue), "bitscore": float(bitscore),
                   "gaps": int(gaps), "subject_query_coverage_pct": float(qcovs),
                   "hsp_query_coverage_pct": round(100 * (abs(int(qend) - int(qstart)) + 1) / int(qlen), 3)}


def best(hits):
    return max(hits, key=lambda h: (h["hsp_query_coverage_pct"], h["identity_pct"],
                                     h["bitscore"]), default=None)


def compact(hit):
    if hit is None:
        return None
    return {key: value for key, value in hit.items() if key != "query_id"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    base = load_base()
    rows = base.load_sequences(base.load_gains())
    catalog, catalog_receipt = reference_provenance.load_catalog()
    panels = {panel: yaml.safe_load((ROOT / "panels" / f"{panel}.yaml").read_text())
              for panel in {r["panel"] for r in rows}}
    query_path = out / "gained_queries.fasta"
    with query_path.open("w") as stream:
        for i, row in enumerate(rows):
            stream.write(f">q{i:03d}\n{row['sequence']}\n")
    source_path = out / "public_sources.fasta"
    source_meta = {}
    with source_path.open("w") as stream:
        for i, (accession, record) in enumerate(sorted(catalog["records"].items())):
            sid = f"source_{i:03d}"
            source_meta[sid] = {"accession": accession, "taxon": record["organism"],
                                "length": len(record["sequence"]),
                                "sha256": record["sequence_sha256"]}
            stream.write(f">{sid}\n{record['sequence']}\n")
    (out / "public_sources_metadata.json").write_text(
        json.dumps(source_meta, indent=2, sort_keys=True) + "\n")
    source_db = out / "public_sources_db"
    run(["makeblastdb", "-in", str(source_path), "-dbtype", "nucl", "-out", str(source_db)])
    fields = "qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore gaps qcovs"
    common = ["-task", "blastn", "-dust", "no", "-evalue", "1e-10",
              "-max_target_seqs", "121", "-max_hsps", "20", "-outfmt", f"6 {fields}"]
    sources_tsv = out / "public_sources_hits.tsv"
    run(["blastn", "-query", str(query_path), "-db", str(source_db),
         *common, "-out", str(sources_tsv)])
    all_source_hits = {f"q{i:03d}": [] for i in range(len(rows))}
    for hit in query_rows(sources_tsv):
        if hit["query_id"] not in all_source_hits or hit["subject_id"] not in source_meta:
            raise ValueError("Unrecognized query or public source in BLAST hits")
        hit.update(source_meta[hit["subject_id"]])
        all_source_hits[hit["query_id"]].append(hit)
    panel_hits = {f"q{i:03d}": [] for i in range(len(rows))}
    panel_receipts = {}
    for panel_name, panel in panels.items():
        refs = reference_provenance.verified_references(panel)
        path = out / f"{panel_name}_references.fasta"
        meta = {}
        with path.open("w") as stream:
            for i, ref in enumerate(refs):
                sid = f"reference_{i:03d}"
                meta[sid] = {"gene": ref["gene_name"],
                             "logical_gene": logical_gene_name(ref["gene_name"], panel),
                             "taxon": ref["taxon"], "accession": ref["accession"],
                             "length": len(ref["sequence"])}
                stream.write(f">{sid}\n{ref['sequence']}\n")
        db = out / f"{panel_name}_references_db"
        run(["makeblastdb", "-in", str(path), "-dbtype", "nucl", "-out", str(db)])
        qpath = out / f"{panel_name}_queries.fasta"
        with qpath.open("w") as stream:
            for i, row in enumerate(rows):
                if row["panel"] == panel_name:
                    stream.write(f">q{i:03d}\n{row['sequence']}\n")
        hits_path = out / f"{panel_name}_reference_hits.tsv"
        run(["blastn", "-query", str(qpath), "-db", str(db),
             *common, "-out", str(hits_path)])
        for hit in query_rows(hits_path):
            if hit["query_id"] not in panel_hits or hit["subject_id"] not in meta:
                raise ValueError("Unrecognized query or panel reference in BLAST hits")
            hit.update(meta[hit["subject_id"]])
            panel_hits[hit["query_id"]].append(hit)
        panel_receipts[panel_name] = {"references": receipt(path), "hits": receipt(hits_path),
                                      "reference_count": len(refs)}
    results = []
    for i, row in enumerate(rows):
        qid = f"q{i:03d}"
        gene = row["gene"].removeprefix(f"{row['panel']}_")
        logical = logical_gene_name(gene, panels[row["panel"]])
        exact = [meta for sid, meta in source_meta.items()
                 if meta["taxon"] == row["taxon"] and
                 (row["sequence"] in catalog["records"][meta["accession"]]["sequence"] or
                  reference_provenance.reverse_complement(row["sequence"]) in
                  catalog["records"][meta["accession"]]["sequence"])]
        result = {key: value for key, value in row.items() if key not in {"sequence", "header"}}
        result.update({
            "query_id": qid,
            "target": gene,
            "best_public_hsp": compact(best(all_source_hits[qid])),
            "best_same_taxon_hsp": compact(best(
                [h for h in all_source_hits[qid] if h["taxon"] == row["taxon"]])),
            "best_panel_hsp": compact(best(panel_hits[qid])),
            "best_expected_target_hsp": compact(best(
                [h for h in panel_hits[qid] if h["logical_gene"] == logical])),
            "exact_full_query_same_taxon_sources": sorted(x["accession"] for x in exact),
            "public_hit_count": len(all_source_hits[qid]),
            "panel_hit_count": len(panel_hits[qid]),
        })
        results.append(result)
    strong = lambda h: h is not None and h["hsp_query_coverage_pct"] >= 90 and h["identity_pct"] >= 90
    summary = {"schema_version": 1, "total": len(results),
               "blastn_version": subprocess.check_output(["blastn", "-version"], text=True).splitlines()[0],
               "makeblastdb_version": subprocess.check_output(["makeblastdb", "-version"], text=True).splitlines()[0],
               "catalog": catalog_receipt,
               "query_fasta": receipt(query_path), "source_fasta": receipt(source_path),
               "source_hits": receipt(sources_tsv), "panel": panel_receipts,
               "settings": {"task": "blastn", "dust": "no", "evalue": "1e-10",
                            "max_target_seqs": 121, "max_hsps": 20,
                            "strong_single_hsp": ">=90% query coverage and >=90% identity"},
               "counts": {"strong_public_single_hsp": sum(strong(x["best_public_hsp"]) for x in results),
                          "strong_same_taxon_single_hsp": sum(strong(x["best_same_taxon_hsp"]) for x in results),
                          "strong_expected_target_single_hsp": sum(strong(x["best_expected_target_hsp"]) for x in results),
                          "exact_full_query_same_taxon": sum(bool(x["exact_full_query_same_taxon_sources"]) for x in results)},
               "by_target": {f"{accession}:{gene}": {
                   "count": len(group),
                   "strong_public_single_hsp": sum(strong(x["best_public_hsp"]) for x in group),
                   "strong_same_taxon_single_hsp": sum(strong(x["best_same_taxon_hsp"]) for x in group),
                   "strong_expected_target_single_hsp": sum(strong(x["best_expected_target_hsp"]) for x in group),
                   "exact_full_query_same_taxon": sum(bool(x["exact_full_query_same_taxon_sources"]) for x in group),
                   "public_hit_any": sum(x["public_hit_count"] > 0 for x in group),
                   "panel_hit_any": sum(x["panel_hit_count"] > 0 for x in group),
               } for accession, gene in sorted({(x["accession"], x["gene"]) for x in results})
                   for group in [[x for x in results if x["accession"] == accession and x["gene"] == gene]]}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (out / "per_sequence.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary["by_target"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
