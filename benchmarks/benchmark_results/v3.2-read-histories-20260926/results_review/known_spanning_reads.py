#!/usr/bin/env python3
"""Describe known alignment records in completed exact-history reports.

This post-outcome helper does not rescan FASTQ, change a sensitivity, or infer
repeat count. Earlier gapped alignment intervals are contextual annotations.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


MAX_REPORT_BYTES = 1024 * 1024 * 1024
PRIOR_REPORT = "dev_docs/ITS_READ_RECRUITMENT_RESULTS.md"
PRIOR_ALIGNMENTS = {
    "SRR27962769.700435": {
        "class_id": "region_1_class_2", "class_read_strand": "-",
        "alignments": [{"layer": "unfiltered", "interval": [17, 141],
                        "minimum_quality": 11, "cigar": "2=1X121="}],
        "note": "No Q20 alignment under the preceding frozen edit limits.",
    },
    "SRR27962769.443018": {
        "class_id": "region_2_class_2", "class_read_strand": "-",
        "alignments": [{"layer": "unfiltered_and_q20", "interval": [1, 148],
                        "minimum_quality": 25, "cigar": "71=1D3=1I10=1X61="}],
        "note": "The preceding minimum-edit traceback is not unique.",
    },
    "SRR27962769.442978": {
        "class_id": "region_2_class_2", "class_read_strand": "-",
        "alignments": [{"layer": "unfiltered", "interval": [1, 148],
                        "minimum_quality": 11, "cigar": "71=1D3=1I10=1X61="},
                       {"layer": "q20", "interval": [1, 147],
                        "minimum_quality": 25, "cigar": "71=1D3=1I10=1X60=1D"}],
        "note": "Earlier Q20 end deletion excluded a Q11 observed endpoint; it is not a biological deletion claim.",
    },
}


def digest_file(path):
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def raw_interval(interval, strand, read_length):
    start, end = interval
    return [start, end] if strand == "+" else [read_length - end, read_length - start]


def overlap(left, right):
    start, end = max(left[0], right[0]), min(left[1], right[1])
    return [start, end] if start < end else None


def prior_span_context(record, classes, k):
    """Project earlier observed spans; coverage does not mean a class match."""
    prior = PRIOR_ALIGNMENTS[record["read_id"]]
    item = classes[prior["class_id"]]
    length = len(record["sequence"])
    contexts = []
    for alignment in prior["alignments"]:
        for representation in item["representations"]:
            if not representation["represented"]:
                continue
            graph_strand = prior["class_read_strand"]
            interval = alignment["interval"]
            if representation["strand"] == "-":
                graph_strand = "+" if graph_strand == "-" else "-"
                interval = [length - interval[1], length - interval[0]]
            # A k-mer is wholly within the aligned observed span iff its start
            # lies in this half-open interval. No gapped hypothesis projection
            # is attempted, because the earlier traceback can be ambiguous.
            starts = [interval[0], max(interval[0], interval[1] - k + 1)]
            layers = {}
            for layer, history in record["histories"][graph_strand].items():
                layers[layer] = {
                    "covering_run_indices": [index for index, run in enumerate(history["runs"])
                                             if run["start"] <= interval[0] and run["end"] >= interval[1]],
                    "overlapping_runs": [{"run_index": index, "run_interval": [run["start"], run["end"]],
                                          "overlap": intersect}
                                         for index, run in enumerate(history["runs"])
                                         if (intersect := overlap(interval, [run["start"], run["end"]])) is not None],
                    "broken_kmer_starts_within_observed_span": [
                        {"interval": intersect, "reasons": interruption["reasons"]}
                        for interruption in history["breaks"]
                        if (intersect := overlap(starts, [interruption["start"], interruption["end"]])) is not None],
                }
            contexts.append({"earlier_alignment": alignment,
                             "class_in_graph_strand": representation["strand"],
                             "graph_read_strand": graph_strand,
                             "observed_interval_in_graph_read_orientation": interval,
                             "raw_observed_interval": raw_interval(interval, graph_strand, length),
                             "whole_span_kmer_start_interval": starts, "layers": layers})
    return {"earlier_report": PRIOR_REPORT, **prior, "projections": contexts}


def summarize_record(record, classes, k):
    length = len(record["sequence"])
    low = [{"position": index, "quality": ord(score) - 33}
           for index, score in enumerate(record["quality"]) if ord(score) - 33 < 20]
    histories = {}
    for strand, layers in record["histories"].items():
        histories[strand] = {}
        for layer, history in layers.items():
            runs = []
            for index, run in enumerate(history["runs"]):
                runs.append({"run_index": index, "oriented_interval": [run["start"], run["end"]],
                             "raw_interval": raw_interval([run["start"], run["end"]], strand, length),
                             "sequence": run["sequence"], "edge_visits": len(run["edges"]),
                             "edges": run["edges"], "branch_events": run["branch_events"],
                             "cyclic_node_visits": run["cyclic_node_visits"],
                             "cyclic_scc_ids": run["cyclic_scc_ids"],
                             "repeated_edge_occurrences": run["repeated_edge_occurrences"]})
            interruptions = []
            for interruption in history["breaks"]:
                interval = [interruption["start"], interruption["end"]]
                # Union of the bases covered by these consecutive invalid
                # windows; this is not a precise causal-base localization.
                covered = [interval[0], interval[1] + k - 1]
                interruptions.append({"oriented_kmer_start_interval": interval,
                                      "window_base_coverage_union": covered,
                                      "raw_window_base_coverage_union": raw_interval(covered, strand, length),
                                      "reasons": interruption["reasons"]})
            histories[strand][layer] = {"runs": runs, "breaks": interruptions}
    return {"read_id": record["read_id"], "ordinal": record["ordinal"],
            "sequence": record["sequence"], "quality": record["quality"],
            "sequence_sha256": record["sequence_sha256"],
            "canonical_sequence_sha256": record["canonical_sequence_sha256"],
            "raw_low_quality_bases": low,
            "raw_invalid_bases": [{"position": index, "base": base}
                                  for index, base in enumerate(record["sequence"]) if base not in "ACGT"],
            "histories": histories, "observations": record["observations"],
            "earlier_alignment_context": prior_span_context(record, classes, k)}


def load_summary(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("Report must be a regular nonsymlink file")
    before = path.stat()
    if before.st_size > MAX_REPORT_BYTES:
        raise ValueError("Report exceeds scanner publication cap")
    source_hash = digest_file(path)
    with path.open() as stream:
        report = json.load(stream)
    if (report.get("status") != "complete" or report.get("schema_version") != 1 or
            report.get("production_admission_changed") is not False):
        raise ValueError("Expected a complete experimental schema-1 report")
    classes = {item["id"]: item for item in report["classes"]}
    records = {}
    for record in report["records"]:
        if record["read_id"] not in PRIOR_ALIGNMENTS:
            continue
        if record["read_id"] in records:
            raise ValueError("Known read ID appears more than once in the retained ledger")
        records[record["read_id"]] = summarize_record(record, classes, report["target"]["kmer_length"])
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
        raise ValueError("Report changed during analysis")
    return {"path": str(path), "sha256": source_hash, "graph_phase": report["graph_phase"],
            "consumed_prefix": report["consumed_prefix"], "sources_before": report["sources_before"],
            "counts": report["counts"], "layer_counts": report["layer_counts"],
            "classes": report["classes"], "records": records,
            "requested_ids_missing_from_retained_ledger": sorted(set(PRIOR_ALIGNMENTS) - set(records))}


def compare(pre, post):
    result = {}
    for read_id in PRIOR_ALIGNMENTS:
        left, right = pre["records"].get(read_id), post["records"].get(read_id)
        if left is None or right is None:
            result[read_id] = {"status": "not_in_both_retained_ledgers"}
            continue
        # Full histories include graph-specific SCC IDs; report their equality
        # separately from the actual ordered sequences and discontinuities.
        result[read_id] = {"status": "compared", "sequence_equal": left["sequence"] == right["sequence"],
                           "quality_equal": left["quality"] == right["quality"],
                           "ordinal_equal": left["ordinal"] == right["ordinal"],
                           "full_histories_equal": left["histories"] == right["histories"], "layers": {}}
        for strand in ("+", "-"):
            result[read_id]["layers"][strand] = {}
            for layer in ("unfiltered", "q20"):
                a, b = left["histories"][strand][layer], right["histories"][strand][layer]
                runs_a = [(run["oriented_interval"], run["edges"]) for run in a["runs"]]
                runs_b = [(run["oriented_interval"], run["edges"]) for run in b["runs"]]
                result[read_id]["layers"][strand][layer] = {
                    "positioned_edge_runs_equal": runs_a == runs_b,
                    "breaks_equal": a["breaks"] == b["breaks"],
                    "pre_run_intervals": [run["oriented_interval"] for run in a["runs"]],
                    "post_run_intervals": [run["oriented_interval"] for run in b["runs"]],
                }
    return result


def duplicate_comparison(summary):
    left = summary["records"].get("SRR27962769.442978")
    right = summary["records"].get("SRR27962769.443018")
    if left is None or right is None:
        return {"status": "not_both_in_retained_ledger"}
    return {"status": "compared", "sequence_equal": left["sequence"] == right["sequence"],
            "canonical_sequence_group_equal": left["canonical_sequence_sha256"] == right["canonical_sequence_sha256"],
            "quality_equal": left["quality"] == right["quality"],
            "unfiltered_histories_equal": all(left["histories"][strand]["unfiltered"] == right["histories"][strand]["unfiltered"] for strand in ("+", "-")),
            "q20_histories_equal": all(left["histories"][strand]["q20"] == right["histories"][strand]["q20"] for strand in ("+", "-")),
            "interpretation_limit": "Quality-dependent exact run boundaries do not create independent sequence groups or positioned molecules."}


def publish(path, document):
    path = Path(path)
    descriptor, temporary = tempfile.mkstemp(prefix=".known-spanning-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as stream:
            json.dump(document, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("pre", "post", "output"):
        parser.add_argument("--" + option, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise ValueError("Output exists; no overwrite permitted")
    pre, post = load_summary(args.pre), load_summary(args.post)
    if pre["consumed_prefix"] != post["consumed_prefix"]:
        raise ValueError("Reports consumed different FASTQ prefixes")
    for identity in ("scanner", "manifest", "expected_prefix"):
        if pre["sources_before"][identity]["sha256"] != post["sources_before"][identity]["sha256"]:
            raise ValueError("Reports used different " + identity + " contents")
    document = {"schema_version": 1, "status": "complete", "helper_sha256": digest_file(Path(__file__)),
                "coordinate_convention": "zero-based half-open; base spans and kmer-start spans explicitly distinguished",
                "pre": pre, "post": post, "known_record_comparison": compare(pre, post),
                "known_duplicate_comparison": {"pre": duplicate_comparison(pre), "post": duplicate_comparison(post)},
                "limitations": ["Reads are selected from retained exact graph/class ledgers, not a new FASTQ scan",
                                "Missing from a retained ledger does not mean missing from the input prefix",
                                "Prior gapped alignments are contextual annotations, not graph-walk evidence",
                                "Break window unions identify covered bases, not unique causal substitutions or indels",
                                "No copy-number, independent-molecule, whole-product, or admission claim"]}
    publish(args.output, document)
    print(json.dumps({"status": "complete", "output": str(args.output),
                      "pre_found_ids": sorted(pre["records"]), "post_found_ids": sorted(post["records"])}, sort_keys=True))


if __name__ == "__main__":
    main()
