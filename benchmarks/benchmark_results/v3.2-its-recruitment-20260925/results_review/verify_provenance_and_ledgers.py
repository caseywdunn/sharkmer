#!/usr/bin/env python3
"""Independently bind frozen provenance and ledgers; never import a scanner.

The FASTQ pass verifies identity and ledger membership only. Control alignments
are exhaustively checked within recorded ledger reads, not across all reads.
"""

import datetime
import hashlib
import json
from pathlib import Path
import sys


BASE = Path(__file__).resolve().parent.parent
OUTPUT = Path(__file__).with_name("provenance_ledger_review.json")


def check(condition, message):
    if not condition:
        raise ValueError(message)


def stat_tuple(path):
    state = path.stat()
    return (state.st_dev, state.st_ino, state.st_size, state.st_mtime_ns, state.st_ctime_ns)


def fingerprint(path):
    path = Path(path)
    check(path.is_file() and not path.is_symlink(), "Nonregular input: " + str(path))
    before = stat_tuple(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    check(before == stat_tuple(path), "Input changed during hashing: " + str(path))
    return {"size_bytes": before[2], "sha256": digest.hexdigest()}


def read_json(path):
    return json.loads(Path(path).read_text())


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def brute_hits(sequence, pattern, limit):
    hits = []
    for start in range(len(sequence) - len(pattern) + 1):
        observed = sequence[start:start + len(pattern)]
        if not set(observed) <= set("ACGT"):
            continue
        differences = sum(a != b for a, b in zip(observed, pattern))
        if differences <= limit:
            hits.append((start, differences))
    return hits


def placements(pattern, candidates, limit):
    result = []
    for candidate_id, sequence in candidates.items():
        for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
            for position, differences in brute_hits(oriented, pattern, limit):
                result.append({"candidate_id": candidate_id, "strand": strand,
                               "start": position, "substitutions": differences})
    return result


def verify_controls(payload, manifest, expected_retained):
    candidates = {item["id"]: item["sequence"] for item in manifest["candidates"]}
    for item in manifest["candidates"]:
        check(hashlib.sha256(item["sequence"].encode("ascii")).hexdigest() == item["sequence_sha256"],
              "Control candidate hash differs")
    limit = payload["max_substitutions"]
    check(payload["minimum_quality"] == 20 and payload["flank_bases"] == 21 and not payload["paired"],
          "Unexpected control quality/flank/pairing parameters")
    windows = {}
    for event in manifest["events"]:
        sequence = candidates[event["candidate_id"]]
        left, right = event["start"] - 21, event["end"] + 21
        check(0 <= left < right <= len(sequence), "Invalid control window")
        windows[event["id"]] = (event, left, right, sequence[left:right])
    by_event = {event_id: [] for event_id in windows}
    placement_checks = 0
    for row in payload["read_evidence"]:
        reconstructed = []
        for event_id, (event, left, right, window) in windows.items():
            matches = []
            for strand, oriented in (("+", window), ("-", reverse_complement(window))):
                for start, differences in brute_hits(row["sequence"], oriented, limit):
                    end = start + len(window)
                    if any(ord(char) - 33 < 20 for char in row["quality"][start:end]):
                        continue
                    observed = row["sequence"][start:end]
                    offset = start if strand == "+" else len(row["sequence"]) - end
                    if strand == "-":
                        observed = reverse_complement(observed)
                    matches.append({"strand": strand, "read_start": start, "read_end": end,
                                    "substitutions": differences,
                                    "projected_footprint": [left - offset, left - offset + len(row["sequence"])],
                                    "observed_template_placements": placements(observed, candidates, limit)})
                    placement_checks += 1
            if matches:
                reconstructed.append({"event_id": event_id, "matches": matches,
                    "candidate_id": event["candidate_id"],
                    "comparison_group": event.get("comparison_group", event["candidate_id"]),
                    "unique": len(matches) == 1 and len(matches[0]["observed_template_placements"]) == 1})
        groups = {}
        for item in reconstructed:
            groups.setdefault(item["comparison_group"], set()).add(item["candidate_id"])
        for item in reconstructed:
            if len(groups[item["comparison_group"]]) > 1:
                item["unique"] = False
        check(reconstructed == row["observations"], "Control ledger placements differ at ordinal " + str(row["ordinal"]))
        for item in reconstructed:
            by_event[item["event_id"]].append((row, item))
    report = []
    for summary in payload["events"]:
        event, left, right, window = windows[summary["id"]]
        rows = by_event[summary["id"]]
        check(summary["window_sequence"] == window and summary["window_start"] == left and summary["window_end"] == right,
              "Manifested control window differs")
        check(summary["template_placements"] == placements(window, candidates, limit), "Control template placements differ")
        counters = {
            "records": len(rows), "unique_records": sum(item["unique"] for _, item in rows),
            "ambiguous_records": sum(not item["unique"] for _, item in rows),
            "exact_records": sum(any(match["substitutions"] == 0 for match in item["matches"]) for _, item in rows),
            "one_substitution_records": sum(any(match["substitutions"] == 1 for match in item["matches"]) for _, item in rows),
        }
        unique_rows = [(row, item) for row, item in rows if item["unique"]]
        ids = {row["header"].split()[0] for row, _ in unique_rows}
        sequences = {min(row["sequence"], reverse_complement(row["sequence"])) for row, _ in unique_rows}
        footprints = {tuple(item["matches"][0]["projected_footprint"]) for _, item in unique_rows}
        strands = {item["matches"][0]["strand"] for _, item in unique_rows}
        counters.update(distinct_ids=len(ids), distinct_sequences=len(sequences), distinct_footprints=len(footprints))
        for name, value in counters.items():
            check(summary[name] == value, "Control summary mismatch: " + name)
        check(summary["footprints"] == [list(item) for item in sorted(footprints)] and summary["strands"] == sorted(strands),
              "Control footprint/strand summary mismatch")
        witness = summary["corroboration_witness"]
        links = {(row["header"].split()[0], tuple(item["matches"][0]["projected_footprint"])) for row, item in unique_rows}
        check(all((item["read_id"], tuple(item["footprint"])) in links for item in witness), "Invalid corroboration witness")
        check(len({item["read_id"] for item in witness}) == len(witness) == len({tuple(item["footprint"]) for item in witness}),
              "Corroboration witness lacks distinct IDs/footprints")
        retained = event["role"] == "supplementary_retained_corresponding_interval"
        check(counters["records"] == (expected_retained if retained else 0), "Historical control count differs")
        check(counters["exact_records"] == (27 if retained else 0), "Historical exact count differs")
        check(counters["one_substitution_records"] == (limit if retained else 0), "Historical substitution count differs")
        check(summary["local_bridge_status"] == ("corroborated_local_bridge" if retained else "unresolved"), "Unexpected control status")
        check(len(witness) >= 3 if retained else not witness, "Unexpected corroboration witness count")
        report.append({"event_id": summary["id"], "role": event["role"], **counters})
    check(len(report) == 2, "Expected both historical control events")
    return {"max_substitutions": limit, "ledger_records": len(payload["read_evidence"]),
            "independently_recomputed_matches": placement_checks, "events": report}


def main():
    freeze_path = BASE / "freeze.json"
    frozen_identity = fingerprint(freeze_path)
    plan = read_json(freeze_path)
    root, run = Path(plan["root"]), Path(plan["run_dir"])
    check(len(plan["sources"]) == 36, "Unexpected frozen source count")
    for relative, expected in plan["sources"].items():
        check(fingerprint(root / relative) == expected, "Frozen source mismatch: " + relative)
        check(fingerprint(run / "source_snapshot" / relative) == expected, "Frozen snapshot mismatch: " + relative)
    execution = read_json(run / "execution.json")
    started = read_json(run / "started.json")
    check(execution["status"] == "complete" and execution["freeze"] == frozen_identity == started["freeze"],
          "Execution/freeze identity mismatch")
    check(started["python"] == plan["python"], "Execution Python differs from freeze")
    check(plan["frozen_utc"] <= started["started_utc"] <= execution["finished_utc"], "Invalid freeze/execution ordering")
    check(execution["commands_planned"] == [item["id"] for item in plan["commands"]], "Planned commands differ")
    check(len(execution["commands"]) == len(plan["commands"]) == 3, "Expected three completed commands")
    payloads = {}
    artifact_hashes = {"freeze.json": frozen_identity, "execution.json": fingerprint(run / "execution.json"),
                       "started.json": fingerprint(run / "started.json")}
    expected_prefix_path = BASE / "inputs/originals/inputs/expected_prefix.json"
    prefix = read_json(expected_prefix_path)
    check(prefix["records"] == 1000000, "Unexpected prefix size")
    prior_finish = started["started_utc"]
    for command, outcome in zip(plan["commands"], execution["commands"]):
        identifier = command["id"]
        receipt_path = run / (identifier + ".receipt.json")
        check(read_json(receipt_path) == outcome, "Inline/standalone command receipt mismatch")
        artifact_hashes[receipt_path.name] = fingerprint(receipt_path)
        check(outcome["id"] == identifier and outcome["status"] == "complete" and outcome["exit_code"] == 0
              and outcome["launched"] and outcome["source_input_bookends_match"], "Command did not complete with verified bookends")
        check(outcome["freeze"] == frozen_identity and outcome["input_before"] == outcome["input_after"] == plan["input_identity"],
              "Command provenance/bookend identity mismatch")
        check(outcome["argv"] == ["/usr/bin/time", "-v", "-o", str(run / (identifier + ".time.txt")), *command["argv"]],
              "Executed command differs from freeze")
        check(prior_finish <= outcome["started_utc"] <= outcome["finished_utc"] <= execution["finished_utc"],
              "Command timestamps inconsistent")
        prior_finish = outcome["finished_utc"]
        check(len(outcome["outputs"]) == 4, "Unexpected command output count")
        for name, expected in outcome["outputs"].items():
            actual = fingerprint(run / name)
            check(actual == expected, "Command output hash mismatch: " + name)
            artifact_hashes[name] = actual
        payload = read_json(run / (identifier + ".json"))
        payloads[identifier] = payload
        check(payload["input_subset"] == prefix, "Scanner prefix identity differs")
        check(payload["provenance"]["scanner_sha256"] == fingerprint(Path(command["argv"][1]))["sha256"], "Scanner provenance hash differs")
        check(payload["provenance"]["expected_prefix_manifest_sha256"] == fingerprint(expected_prefix_path)["sha256"], "Prefix manifest hash differs")
    recruitment = payloads["recruitment"]
    check(recruitment["provenance"]["argv"] == plan["commands"][0]["argv"][1:], "Recruitment argv differs")
    check(recruitment["provenance"]["manifest_sha256"] == fingerprint(BASE / "inputs/manifest.json")["sha256"], "Recruitment manifest hash differs")
    check(recruitment["status"] == "completed", "Recruitment did not complete")
    check(recruitment["staged_file_before"] == recruitment["staged_file_after"] == recruitment["provenance"]["staged_file_before"],
          "Recruitment stat/hash bookends differ")
    for key, value in plan["input_identity"].items():
        check(recruitment["staged_file_after"][key] == value, "Recruitment bookend identity differs from freeze")
    ledgers = {}
    for identifier, payload in payloads.items():
        ledgers[identifier] = {row["ordinal"]: row for row in payload["read_evidence"]}
        check(len(ledgers[identifier]) == len(payload["read_evidence"]), "Duplicate ledger ordinal")
        check(list(ledgers[identifier]) == sorted(ledgers[identifier]), "Unsorted ledger ordinals")
    wanted = set().union(*(set(rows) for rows in ledgers.values()))
    matched = {identifier: 0 for identifier in ledgers}
    fastq = Path(plan["fastq"])
    check(not fastq.is_symlink(), "FASTQ unexpectedly became a symlink")
    before = stat_tuple(fastq)
    digest = hashlib.sha256()
    bases = size = maximum = 0
    with fastq.open("rb") as stream:
        for ordinal in range(1, prefix["records"] + 1):
            lines = [stream.readline() for _ in range(4)]
            check(all(lines), "Premature FASTQ end")
            raw = b"".join(lines)
            digest.update(raw)
            size += len(raw)
            header, sequence, separator, quality = [line.rstrip(b"\r\n") for line in lines]
            check(header.startswith(b"@") and separator.startswith(b"+") and len(sequence) == len(quality) > 0,
                  "Malformed FASTQ record")
            bases += len(sequence)
            maximum = max(maximum, len(sequence))
            if ordinal not in wanted:
                continue
            for identifier, rows in ledgers.items():
                if ordinal not in rows:
                    continue
                row = rows[ordinal]
                check(row["header"] == header.decode("ascii") and row["sequence"] == sequence.decode("ascii")
                      and row["quality"] == quality.decode("ascii"), "Ledger record differs from original FASTQ")
                if "read_id" in row:
                    check(row["read_id"] == row["header"].split()[0], "Ledger read ID differs")
                matched[identifier] += 1
        check(stream.read(1) == b"", "Staged FASTQ contains reads beyond the frozen prefix")
    observed = {"records": prefix["records"], "bases": bases, "size_bytes": size, "sha256": digest.hexdigest()}
    check(observed == prefix, "Streamed FASTQ identity differs from expected prefix")
    check(before == stat_tuple(fastq), "FASTQ changed during ledger verification")
    check({key: observed[key] for key in ("sha256", "size_bytes")} == plan["input_identity"], "Input differs from freeze")
    bookend = recruitment["staged_file_after"]
    check(before == (bookend["device"], bookend["inode"], bookend["size_bytes"], bookend["mtime_ns"], bookend["ctime_ns"]),
          "Current input metadata differs from production bookend")
    for identifier, rows in ledgers.items():
        check(matched[identifier] == len(rows), "Ledger record not found in prefix")
        check(payloads[identifier]["maximum_read_length"] == maximum, "Maximum read length differs")
    controls = []
    manifest_path = BASE / "inputs/originals/inputs/co1_control_manifest.json"
    manifest = read_json(manifest_path)
    for substitutions in (0, 1):
        payload = payloads["co1_substitutions_" + str(substitutions)]
        check(payload["provenance"]["event_manifest_sha256"] == fingerprint(manifest_path)["sha256"], "Control manifest hash differs")
        check(payload["max_substitutions"] == substitutions, "Control substitution limit differs")
        controls.append(verify_controls(payload, manifest, 27 + substitutions))
    check(set(ledgers["co1_substitutions_0"]) < set(ledgers["co1_substitutions_1"]), "Exact-control reads not preserved by tolerance")
    # End-of-review recheck protects the identities reported in this review.
    check(fingerprint(freeze_path) == frozen_identity, "Freeze changed during review")
    for relative, expected in plan["sources"].items():
        check(fingerprint(root / relative) == fingerprint(run / "source_snapshot" / relative) == expected,
              "Source or snapshot changed during review")
    for name, expected in artifact_hashes.items():
        check(fingerprint(freeze_path if name == "freeze.json" else run / name) == expected, "Artifact changed during review")
    return {"status": "pass", "reviewed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "review_source": fingerprint(Path(__file__)), "frozen_sources_verified": len(plan["sources"]),
            "snapshots_verified": len(plan["sources"]), "commands_verified": 3,
            "artifacts": artifact_hashes, "streamed_prefix": observed, "maximum_read_length": maximum,
            "ledger_records_bound": matched, "unique_ordinals_bound": len(wanted), "controls": controls,
            "checks": ["source and snapshot hashes before and after review", "freeze, argv, command receipts, output hashes and input bookends",
                       "each ledger ordinal/header/sequence/quality bound to the original frozen prefix",
                       "exhaustive control placements on each recorded control read, Q20, strands, footprints and competitor uniqueness",
                       "historical control counts retained 27 exact / 28 with one substitution; lost 0 in both"],
            "limitations": ["Ledger binding verifies every reported record; this review does not independently establish absence of unreported matches among non-ledger reads.",
                            "No production scanner was imported or rerun; only new results_review files were written."]}


if __name__ == "__main__":
    try:
        result = main()
    except Exception as error:
        result = {"status": "fail", "error_type": type(error).__name__, "error": str(error),
                  "reviewed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in ("status", "frozen_sources_verified", "snapshots_verified", "ledger_records_bound", "unique_ordinals_bound", "error") if key in result}, sort_keys=True))
    sys.exit(0 if result["status"] == "pass" else 1)
