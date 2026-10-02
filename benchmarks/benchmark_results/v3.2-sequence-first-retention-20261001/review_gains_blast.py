#!/usr/bin/env python3
"""Independently check the saved BLAST review of newly retained sequences."""

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BUNDLE = Path(__file__).resolve().parent
RESULTS = BUNDLE / "blast-gains-v3"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def fasta(path):
    records = {}
    name = None
    pieces = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                assert name not in records, (path, name)
                records[name] = "".join(pieces)
            name = line[1:].split()[0]
            pieces = []
        elif line.strip():
            assert name is not None, path
            pieces.append(line.strip())
    if name is not None:
        assert name not in records, (path, name)
        records[name] = "".join(pieces)
    return records


def check_receipt(receipt):
    path = ROOT / receipt["path"]
    data = path.read_bytes()
    assert len(data) == receipt["size_bytes"], path
    assert digest(data) == receipt["sha256"], path


def read_hits(path, queries, subjects):
    hit_count = Counter()
    high_coverage = []
    high_identity = 0
    strong = []
    with path.open(newline="") as stream:
        for fields in csv.reader(stream, delimiter="\t"):
            assert len(fields) == 14, (path, fields)
            qid, sid, identity, _, qlen, qstart, qend, slen, *_rest = fields
            assert qid in queries and sid in subjects, (path, qid, sid)
            assert int(qlen) == len(queries[qid]), (path, qid)
            assert int(slen) == len(subjects[sid]), (path, sid)
            assert 1 <= int(qstart) <= int(qlen) and 1 <= int(qend) <= int(qlen)
            pct_coverage = 100 * (abs(int(qend) - int(qstart)) + 1) / int(qlen)
            pct_identity = float(identity)
            assert 0 <= pct_identity <= 100
            hit_count[qid] += 1
            if pct_coverage >= 90:
                high_coverage.append((qid, sid, round(pct_coverage, 3), pct_identity))
            if pct_identity >= 90:
                high_identity += 1
            if pct_coverage >= 90 and pct_identity >= 90:
                strong.append((qid, sid))
    return hit_count, high_coverage, high_identity, strong


def main():
    if not __debug__:
        raise RuntimeError("Run this evidence audit without Python optimization")
    summary = json.loads((RESULTS / "summary.json").read_text())
    gains = list(csv.DictReader((BUNDLE / "report-v3/uncertain_gains.tsv").open(), delimiter="\t"))
    assert len(gains) == summary["total"] == 200
    assert len({(r["accession"], r["gene"], r["sha256"]) for r in gains}) == 200
    queries = fasta(RESULTS / "gained_queries.fasta")
    assert len(queries) == 200
    rows = json.loads((RESULTS / "per_sequence.json").read_text())
    assert len(rows) == 200
    for i, (gain, row) in enumerate(zip(gains, rows)):
        qid = f"q{i:03d}"
        assert row["query_id"] == qid and qid in queries
        assert all(row[k] == gain[k] for k in ("accession", "gene", "sha256"))
        assert len(queries[qid]) == int(gain["length"]) == row["length"]
        assert digest(queries[qid].encode("ascii")) == gain["sha256"]

    catalog_path = ROOT / "panels/reference_sources.json.gz"
    assert digest(catalog_path.read_bytes()) == summary["catalog"]["sha256"]
    catalog = json.loads(gzip.decompress(catalog_path.read_bytes()))["records"]
    sources = fasta(RESULTS / "public_sources.fasta")
    metadata = json.loads((RESULTS / "public_sources_metadata.json").read_text())
    assert len(sources) == len(metadata) == len(catalog) == 121
    for sid, entry in metadata.items():
        original = catalog[entry["accession"]]
        assert sources[sid] == original["sequence"]
        assert entry["taxon"] == original["organism"]
        assert entry["length"] == len(sources[sid])
        assert entry["sha256"] == digest(sources[sid].encode("ascii"))

    for key in ("query_fasta", "source_fasta", "source_hits"):
        check_receipt(summary[key])
    public_counts, public_high_cov, public_high_id, public_strong = read_hits(
        RESULTS / "public_sources_hits.tsv", queries, sources
    )
    assert not public_strong
    assert sum(bool(public_counts[q]) for q in queries) == 168
    assert all(row["public_hit_count"] == public_counts[row["query_id"]] for row in rows)

    panel_results = {}
    panel_counts = Counter()
    for panel, receipts in sorted(summary["panel"].items()):
        check_receipt(receipts["references"])
        check_receipt(receipts["hits"])
        subjects = fasta(RESULTS / f"{panel}_references.fasta")
        assert len(subjects) == receipts["reference_count"]
        subset = {row["query_id"]: queries[row["query_id"]] for row in rows if row["panel"] == panel}
        counts, high_cov, high_id, strong = read_hits(
            RESULTS / f"{panel}_reference_hits.tsv", subset, subjects
        )
        assert not strong
        panel_counts.update(counts)
        panel_results[panel] = {"references": len(subjects), "hsp_count": sum(counts.values()),
                                "queries_with_hits": len(counts),
                                "hsps_at_least_90pct_coverage": len(high_cov),
                                "hsps_at_least_90pct_identity": high_id,
                                "hsps_meeting_both": len(strong)}
    assert all(row["panel_hit_count"] == panel_counts[row["query_id"]] for row in rows)

    exact = 0
    same_taxon_sources = defaultdict(list)
    for sid, source in sources.items():
        same_taxon_sources[metadata[sid]["taxon"]].append(source)
    complement = str.maketrans("ACGT", "TGCA")
    by_target = defaultdict(lambda: {"queries": 0, "public_hits": 0, "panel_hits": 0})
    for row in rows:
        sequence = queries[row["query_id"]]
        revcomp = sequence.translate(complement)[::-1]
        matches = [s for s in same_taxon_sources[row["taxon"]]
                   if sequence in s or revcomp in s]
        assert bool(matches) == bool(row["exact_full_query_same_taxon_sources"])
        exact += bool(matches)
        group = by_target[f"{row['accession']}:{row['gene']}"]
        group["queries"] += 1
        group["public_hits"] += bool(public_counts[row["query_id"]])
        group["panel_hits"] += bool(panel_counts[row["query_id"]])
    assert exact == 0
    assert summary["counts"] == {"strong_public_single_hsp": 0,
                                  "strong_same_taxon_single_hsp": 0,
                                  "strong_expected_target_single_hsp": 0,
                                  "exact_full_query_same_taxon": 0}
    assert len(public_high_cov) == 5
    review = {"status": "passed", "scope": "saved v3 local BLAST outputs and pinned source catalog",
              "queries": len(queries), "sources": len(sources),
              "public_hsps": sum(public_counts.values()),
              "public_queries_with_hits": len(public_counts),
              "public_hsps_at_least_90pct_coverage": len(public_high_cov),
              "public_hsps_at_least_90pct_identity": public_high_id,
              "public_hsps_meeting_both": len(public_strong),
              "exact_full_query_same_taxon": exact,
              "high_coverage_hsps": [{"query_id": q, "source_accession": metadata[s]["accession"],
                                      "query_taxon": rows[int(q[1:])]["taxon"],
                                      "source_taxon": metadata[s]["taxon"],
                                      "query_coverage_pct": cov, "identity_pct": identity}
                                     for q, s, cov, identity in public_high_cov],
              "panels": panel_results, "by_target": dict(sorted(by_target.items())),
              "interpretation": "No qualifying alignment in this bounded catalog; biological validity and read support remain unassessed."}
    (RESULTS / "review.json").write_text(json.dumps(review, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: review[k] for k in ("status", "queries", "public_hsps",
                                              "public_hsps_meeting_both", "exact_full_query_same_taxon")},
                     indent=2))


if __name__ == "__main__":
    main()
