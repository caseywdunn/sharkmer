#!/usr/bin/env python3
"""Describe saved uncertain gains without assembly, references, or read scans."""

import argparse
import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path


K_VALUES = (19, 31, 51)
LIMITATIONS = [
    "Only the 200 saved gained sequences are compared; the nearest neighbor is restricted to gains in the same accession and target.",
    "Canonical k-mer sets ignore multiplicity, order, orientation, and long-range linkage; identical sets do not establish identical biological products.",
    "Differing windows use the longest exact common prefix, then a common suffix that does not overlap that prefix. This is a single-window representation, not an alignment or edit distance.",
    "Nearest means minimizing maximum allele-window length, then total allele-window length, then neighbor SHA-256. Original emitted sequence orientation is used.",
    "Flanks are up to 19 exact common bases per side. Their uniqueness and suitability for read recruitment have not been evaluated; repeat-related indels may have other equally valid placements.",
    "No read support, biological correctness, repeat copy number, full-product phase, or target specificity is established. No production behavior changes.",
]


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical_kmers(sequence, k):
    return frozenset(
        min(sequence[start:start + k], reverse_complement(sequence[start:start + k]))
        for start in range(len(sequence) - k + 1)
    )


def differing_window(first, second):
    """Return exact shared flanks and half-open differing windows in each input."""
    prefix = 0
    while prefix < min(len(first), len(second)) and first[prefix] == second[prefix]:
        prefix += 1
    suffix = 0
    while (suffix < min(len(first), len(second)) - prefix
           and first[len(first) - suffix - 1] == second[len(second) - suffix - 1]):
        suffix += 1
    first_end, second_end = len(first) - suffix, len(second) - suffix
    return {
        "common_prefix_length": prefix,
        "common_suffix_length": suffix,
        "query_start": prefix,
        "query_end": first_end,
        "neighbor_start": prefix,
        "neighbor_end": second_end,
        "query_allele": first[prefix:first_end],
        "neighbor_allele": second[prefix:second_end],
        "left_flank": first[max(0, prefix - 19):prefix],
        "right_flank": first[first_end:min(len(first), first_end + 19)],
        "max_allele_length": max(first_end - prefix, second_end - prefix),
        "total_allele_length": first_end + second_end - 2 * prefix,
    }


def parse_fasta(text):
    records = {}
    name, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                if name in records:
                    raise ValueError("duplicate FASTA ID")
                records[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        elif line.strip():
            if name is None:
                raise ValueError("sequence before FASTA header")
            chunks.append(line.strip())
    if name is not None:
        if name in records:
            raise ValueError("duplicate FASTA ID")
        records[name] = "".join(chunks)
    return records


def bind_inputs(gains, mappings, sequences, expected_count=200):
    def key(row):
        return row["accession"], row["gene"], int(row["length"]), row["sha256"]

    gain_keys = [key(row) for row in gains]
    map_keys = [key(row) for row in mappings]
    query_ids = [row["query_id"] for row in mappings]
    if (len(gain_keys) != expected_count or len(set(gain_keys)) != expected_count
            or len(map_keys) != expected_count or len(set(map_keys)) != expected_count
            or set(gain_keys) != set(map_keys) or len(set(query_ids)) != expected_count
            or set(query_ids) != set(sequences)):
        raise ValueError("gain table, query mapping, and FASTA must bind one-to-one")
    result = []
    for row in sorted(mappings, key=lambda row: row["query_id"]):
        sequence = sequences[row["query_id"]]
        if not sequence or set(sequence) - set("ACGT"):
            raise ValueError("query contains empty or non-ACGT sequence")
        if len(sequence) != int(row["length"]) or sha256_bytes(sequence.encode()) != row["sha256"]:
            raise ValueError("query sequence length/hash mismatch")
        result.append({name: row[name] for name in ("query_id", "accession", "gene", "sha256")}
                      | {"length": len(sequence), "sequence": sequence})
    return result


def analyze(records):
    by_target = defaultdict(list)
    for record in records:
        by_target[(record["accession"], record["gene"])].append(record)
    targets, windows = [], []
    for (accession, gene), group in sorted(by_target.items()):
        lengths = Counter(record["length"] for record in group)
        target = {"accession": accession, "gene": gene, "gains": len(group),
                  "min_length": min(lengths), "max_length": max(lengths),
                  "length_counts": dict(sorted(lengths.items())), "kmer_sets": []}
        for k in K_VALUES:
            set_groups = defaultdict(list)
            for record in group:
                set_groups[canonical_kmers(record["sequence"], k)].append(record)
            groups = []
            for kmers, members in set_groups.items():
                if len(members) < 2:
                    continue
                groups.append({
                    "set_sha256": sha256_bytes(("\n".join(sorted(kmers)) + "\n").encode()),
                    "distinct_kmers": len(kmers),
                    "query_ids": sorted(record["query_id"] for record in members),
                    "lengths": sorted({record["length"] for record in members}),
                })
            target["kmer_sets"].append({"k": k, "distinct_sets": len(set_groups),
                                        "shared_groups": sorted(groups, key=lambda row: row["query_ids"])})
        target_windows = []
        for query in group:
            choices = []
            for neighbor in group:
                if neighbor["query_id"] == query["query_id"]:
                    continue
                window = differing_window(query["sequence"], neighbor["sequence"])
                choices.append(((window["max_allele_length"], window["total_allele_length"],
                                 neighbor["sha256"]), neighbor, window))
            if not choices:
                continue
            score, neighbor, window = min(choices, key=lambda choice: choice[0])
            row = {"query_id": query["query_id"], "accession": accession, "gene": gene,
                   "query_sha256": query["sha256"], "query_length": query["length"],
                   "neighbor_query_id": neighbor["query_id"],
                   "neighbor_sha256": neighbor["sha256"], "neighbor_length": neighbor["length"],
                   "equal_score_neighbors_before_hash_tiebreak": sum(
                       choice[0][:2] == score[:2] for choice in choices)} | window
            windows.append(row)
            target_windows.append(row)
        target["nearest_window_max_allele_length_counts"] = dict(sorted(Counter(
            row["max_allele_length"] for row in target_windows).items()))
        target["queries_with_nearest_window_at_most_3_bases"] = sum(
            row["max_allele_length"] <= 3 for row in target_windows)
        target["queries_with_nearest_window_at_most_5_bases"] = sum(
            row["max_allele_length"] <= 5 for row in target_windows)
        target["largest_nearest_window"] = max(
            (row["max_allele_length"] for row in target_windows), default=None)
        target["windows_with_two_full_19bp_flanks"] = sum(
            len(row["left_flank"]) == len(row["right_flank"]) == 19 for row in target_windows)
        target["windows_with_edge_limited_flanks"] = (
            len(target_windows) - target["windows_with_two_full_19bp_flanks"])
        targets.append(target)
    return targets, windows


def tsv_bytes(rows, fields):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode()


def write_results(output, payloads):
    """Refuse existing result paths before creating any output file."""
    existing = [name for name in payloads if (output / name).exists()]
    if existing:
        raise FileExistsError("refusing to overwrite existing results: " + ", ".join(existing))
    output.mkdir(parents=True, exist_ok=True)
    for name, data in payloads.items():
        with (output / name).open("xb") as handle:
            handle.write(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1]
                        / "v3.2-sequence-first-retention-20261001")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    files = {
        "gains": args.source / "report-v3/uncertain_gains.tsv",
        "mapping": args.source / "blast-gains-v3/per_sequence.json",
        "fasta": args.source / "blast-gains-v3/gained_queries.fasta",
    }
    raw = {name: path.read_bytes() for name, path in files.items()}
    records = bind_inputs(list(csv.DictReader(io.StringIO(raw["gains"].decode()), delimiter="\t")),
                          json.loads(raw["mapping"]), parse_fasta(raw["fasta"].decode()))
    targets, windows = analyze(records)
    window_tsv = tsv_bytes(windows, list(windows[0]))
    fields = ["accession", "gene", "gains", "min_length", "max_length",
              "k19_distinct_sets", "k31_distinct_sets", "k51_distinct_sets",
              "queries_with_nearest_window_at_most_3_bases",
              "queries_with_nearest_window_at_most_5_bases", "largest_nearest_window",
              "windows_with_two_full_19bp_flanks", "windows_with_edge_limited_flanks"]
    target_rows = []
    for target in targets:
        row = {key: target[key] for key in fields if key in target}
        row.update({f"k{entry['k']}_distinct_sets": entry["distinct_sets"]
                    for entry in target["kmer_sets"]})
        target_rows.append(row)
    target_tsv = tsv_bytes(target_rows, fields)
    report = {
        "schema_version": 1,
        "purpose": "saved_sequence_structure_triage_without_read_support_evaluation",
        "coordinate_system": "zero_based_half_open_in_emitted_oriented_sequences",
        "source_files": {name: {"path": str(files[name]), "bytes": len(value),
                                "sha256": sha256_bytes(value)} for name, value in raw.items()},
        "producer_sha256": sha256_bytes(Path(__file__).read_bytes()),
        "verified_gain_queries": len(records), "target_groups": len(targets),
        "limitations": LIMITATIONS, "targets": targets,
        "nearest_window_count": len(windows),
        "largest_nearest_window": max(row["max_allele_length"] for row in windows),
        "queries_with_nearest_window_at_most_5_bases": sum(
            row["max_allele_length"] <= 5 for row in windows),
        "windows_with_two_full_19bp_flanks": sum(
            len(row["left_flank"]) == len(row["right_flank"]) == 19 for row in windows),
        "windows_with_edge_limited_flanks": sum(
            len(row["left_flank"]) < 19 or len(row["right_flank"]) < 19 for row in windows),
        "result_files": {
            name: {"path": str(args.output / name), "bytes": len(data), "sha256": sha256_bytes(data)}
            for name, data in {"nearest_windows.tsv": window_tsv, "target_summary.tsv": target_tsv}.items()
        },
    }
    write_results(args.output, {
        "nearest_windows.tsv": window_tsv,
        "target_summary.tsv": target_tsv,
        "report.json": (json.dumps(report, indent=2, sort_keys=True) + "\n").encode(),
    })
    print(json.dumps({"verified_gain_queries": len(records), "target_groups": len(targets),
                      "nearest_windows": len(windows), "output": str(args.output)}))


if __name__ == "__main__":
    main()
