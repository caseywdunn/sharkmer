#!/usr/bin/env python3
"""Independent post-run provenance, raw-record, graph-history and class audit.

No prototype code is imported. Exhaustive k-mer windows, literal placements,
and a residual-network matching check independently reconstruct the evidence.
Run only after the separate graph experiment has completed.
"""

import argparse
import ast
from collections import Counter, defaultdict, deque
import datetime
import hashlib
from itertools import groupby
import json
import os
from pathlib import Path
import tarfile
import tempfile
import time


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent
ROOT = EVIDENCE.parents[2]
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
EXPECTED_PREFIX = {"records": 1000000, "bases": 150000000, "size_bytes": 365356898,
                   "sha256": "6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794"}
ARCHIVE_SHA = "14799e1b8da1c1efc8675647a8c5a26076edde9894e55107f37db0ec0ba6abc3"
GRAPH_SHA = {
    "pre": "0810b40673c4e4bdc7dcad0573bbe7f3f49e8a0b4542ff5789bd7ce2ec545bf9",
    "post": "b23e89f11d6ed653c7d8113a3f2c07506932049c8ed2397d4083a2f5fd5dfd40",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def identity(path):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), f"Nonregular input: {path}")
    before = path.stat()
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    after = path.stat()
    require(all(getattr(before, key) == getattr(after, key) for key in
                ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")), f"Input changed while hashing: {path}")
    return {"size_bytes": after.st_size, "sha256": hasher.hexdigest()}


def basic(value):
    return {key: value[key] for key in ("size_bytes", "sha256")}


def read_json(path):
    return json.loads(Path(path).read_text())


def canonical(values):
    return sorted(values, key=lambda item: json.dumps(item, sort_keys=True))


def publish(path, value):
    raw = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    descriptor, temporary = tempfile.mkstemp(prefix=".history-review-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(raw)
            output.flush()
            os.fsync(output.fileno())
        os.link(temporary, path)
    finally:
        os.unlink(temporary)


def frozen_limits(path):
    def number(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return node.value
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
            return number(node.left) * number(node.right)
        raise ValueError("Unsupported frozen integer limit expression")
    values = {}
    for node in ast.parse(Path(path).read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name.startswith("MAX_"):
                values[name] = number(node.value)
    return values


def graph_model(document):
    k = document["target"]["kmer_length"]
    nodes = {item["node_id"]: item for item in document["graph"]["nodes"]}
    edges = {item["edge_id"]: item for item in document["graph"]["edges"]}
    require(len(nodes) == len(document["graph"]["nodes"]) and len(edges) == len(document["graph"]["edges"]), "Duplicate graph identities")
    require(len(nodes) == document["graph"]["node_count"] and len(edges) == document["graph"]["edge_count"], "Graph identity/count mismatch")
    incoming, outgoing = Counter(), Counter()
    for sequence, edge in edges.items():
        require(sequence == edge["oriented_sequence"] and len(sequence) == k, "Invalid graph edge spelling")
        require(edge["source_node_id"] == sequence[:-1] and edge["target_node_id"] == sequence[1:], "Invalid graph edge endpoints")
        require(sequence[:-1] in nodes and sequence[1:] in nodes, "Missing graph node")
        outgoing[sequence[:-1]] += 1
        incoming[sequence[1:]] += 1
    for sequence, node in nodes.items():
        require(sequence == node["oriented_sequence"] and len(sequence) == k - 1, "Invalid graph node spelling")
        require(incoming[sequence] == node["incoming_degree"] and outgoing[sequence] == node["outgoing_degree"], "Graph degree mismatch")
    return {"k": k, "nodes": nodes, "edges": edges, "incoming": incoming, "outgoing": outgoing}


def history_oracle(sequence, quality, graph, filtered):
    """Partition all window starts by exact rejection reason, then spell runs."""
    k = graph["k"]
    windows = []
    for start in range(max(0, len(sequence) - k + 1)):
        word = sequence[start:start + k]
        reasons = []
        if any(base not in "ACGT" for base in word):
            reasons.append("invalid_base")
        if filtered and any(ord(score) < 53 for score in quality[start:start + k]):
            reasons.append("low_quality")
        if word not in graph["edges"]:
            reasons.append("missing_graph_kmer")
        windows.append((start, tuple(reasons)))
    result = {"runs": [], "breaks": []}
    for reason, group in groupby(windows, key=lambda item: item[1]):
        offsets = [item[0] for item in group]
        start, stop = offsets[0], offsets[-1] + 1
        if reason:
            result["breaks"].append({"start": start, "end": stop, "reasons": list(reason)})
            continue
        end = stop + k - 1
        edges = [sequence[offset:offset + k] for offset in offsets]
        require(all(first[1:] == second[:-1] for first, second in zip(edges, edges[1:])), "Oracle run disconnected")
        spelled = edges[0] + "".join(edge[-1] for edge in edges[1:])
        require(spelled == sequence[start:end], "Oracle run spelling differs")
        branch = []
        for index in range(len(edges) - 1):
            node = edges[index][1:]
            if graph["incoming"][node] > 1 or graph["outgoing"][node] > 1:
                branch.append({"edge_offset": index, "read_kmer_start": start + index,
                               "incoming_edge": edges[index], "outgoing_edge": edges[index + 1]})
        cyclic = []
        run_nodes = [edges[0][:-1]] + [edge[1:] for edge in edges]
        for index, node in enumerate(run_nodes):
            component = graph["nodes"][node].get("cyclic_scc_id")
            if component is not None:
                cyclic.append({"node_offset": index, "read_node_start": start + index,
                               "node_id": node, "cyclic_scc_id": component})
        result["runs"].append({"start": start, "end": end, "sequence": spelled, "edges": edges,
            "branch_events": branch, "cyclic_node_visits": cyclic,
            "repeated_edge_occurrences": len(edges) - len(set(edges)),
            "cyclic_scc_ids": sorted({item["cyclic_scc_id"] for item in cyclic}, key=str)})
    return result


def representations(sequence, graph):
    k = graph["k"]
    result = []
    for strand, spelling in (("+", sequence), ("-", rc(sequence))):
        missing = [start for start in range(len(spelling) - k + 1) if spelling[start:start + k] not in graph["edges"]]
        anchors = [spelling[:k - 1], spelling[-(k - 1):]]
        result.append({"strand": strand, "represented": not missing, "missing_edge_starts": missing,
                       "anchor_nodes": anchors, "noncyclic_anchors": all(node in graph["nodes"] and graph["nodes"][node].get("cyclic_scc_id") is None for node in anchors)})
    return result


def placements(sequence, quality, classes, histories):
    result = []
    length = len(sequence)
    for item in classes:
        literal = item["sequence"]
        for strand, oriented, scores in (("+", sequence, quality), ("-", rc(sequence), quality[::-1])):
            for start in range(len(oriented) - len(literal) + 1):
                if not oriented.startswith(literal, start):
                    continue
                end = start + len(literal)
                physical = [start, end] if strand == "+" else [length - end, length - start]
                paths = {"unfiltered": [], "q20": []}
                anchors = False
                for representation in item["representations"]:
                    graph_strand = strand if representation["strand"] == "+" else ("+" if strand == "-" else "-")
                    interval = physical if graph_strand == "+" else [length - physical[1], length - physical[0]]
                    for layer in paths:
                        for index, run in enumerate(histories[graph_strand][layer]["runs"]):
                            if run["start"] <= interval[0] and interval[1] <= run["end"]:
                                paths[layer].append({"graph_strand": graph_strand, "run_index": index,
                                    "class_in_graph_strand": representation["strand"], "interval": interval})
                                anchors |= layer == "q20" and representation["noncyclic_anchors"]
                result.append({"region_id": item["region_id"], "class_id": item["id"],
                    "class_read_strand": strand, "class_oriented_start": start, "physical_interval": physical,
                    "footprint": [-start, length - start], "minimum_quality": min(map(ord, scores[start:end])) - 33,
                    "paths": paths, "noncyclic_q20_anchors": anchors})
    counts = Counter(item["region_id"] for item in result)
    for item in result:
        item["unique_unfiltered_placement"] = counts[item["region_id"]] == 1
        item["qualifying_local_witness"] = bool(counts[item["region_id"]] == 1 and item["minimum_quality"] >= 20 and item["noncyclic_q20_anchors"] and item["paths"]["q20"])
    return result


def matching_size_three(links):
    """Independent residual-network breadth-first matching, stop at three."""
    neighbors, residual = defaultdict(set), {}
    source, sink = ("source",), ("sink",)
    def arc(left, right):
        neighbors[left].add(right)
        neighbors[right].add(left)
        residual[left, right], residual[right, left] = 1, 0
    for read_id, footprints in links.items():
        arc(source, ("id", read_id))
        for footprint in footprints:
            arc(("id", read_id), ("footprint", footprint))
    for footprint in set().union(*links.values()) if links else set():
        arc(("footprint", footprint), sink)
    flow = 0
    while flow < 3:
        parent, queue = {source: None}, deque([source])
        while queue and sink not in parent:
            left = queue.popleft()
            for right in neighbors[left]:
                if right not in parent and residual[left, right]:
                    parent[right] = left
                    queue.append(right)
        if sink not in parent:
            break
        cursor = sink
        while cursor != source:
            previous = parent[cursor]
            residual[previous, cursor] -= 1
            residual[cursor, previous] += 1
            cursor = previous
        flow += 1
    return flow


def audit_class_summary(item, summary, records):
    require(all(summary[key] == item.get(key) for key in ("id", "region_id", "sequence", "kind")), "Class summary metadata differs")
    require(summary["representations"] == item["representations"], "Class graph representability mismatch")
    exact, q20, ambiguous = set(), set(), set()
    links, groups, qualifying = defaultdict(set), defaultdict(set), []
    for record in records:
        for observation in record["observations"]:
            if observation["class_id"] != item["id"]:
                continue
            exact.add(record["ordinal"])
            if observation["minimum_quality"] >= 20:
                q20.add(record["ordinal"])
            if not observation["unique_unfiltered_placement"]:
                ambiguous.add(record["ordinal"])
            if observation["qualifying_local_witness"]:
                footprint = tuple(observation["footprint"])
                links[record["read_id"]].add(footprint)
                groups[record["canonical_sequence_sha256"]].add(footprint)
                qualifying.append({"ordinal": record["ordinal"], "read_id": record["read_id"],
                    "canonical_sequence_sha256": record["canonical_sequence_sha256"],
                    "footprint": list(footprint), "class_read_strand": observation["class_read_strand"]})
    require(all(len(values) == 1 for values in groups.values()), "Canonical sequence has inconsistent unique footprints")
    fields = {"exact_records": len(exact), "q20_exact_records": len(q20), "ambiguous_exact_records": len(ambiguous),
              "qualifying_records": len(qualifying), "distinct_read_ids": len(links),
              "distinct_footprints": len(set().union(*links.values())) if links else 0}
    require(all(summary[key] == value for key, value in fields.items()), f"Class summary counts differ: {item['id']}")
    expected_groups = [{"sha256": group, "footprints": [list(value) for value in sorted(values)]} for group, values in sorted(groups.items())]
    require(summary["canonical_sequence_groups"] == expected_groups, "Sequence group summary differs")
    witness = summary["witness"]
    require(len({member["read_id"] for member in witness}) == len(witness), "Witness repeats read ID")
    require(len({tuple(member["footprint"]) for member in witness}) == len(witness), "Witness repeats footprint")
    for member in witness:
        require(canonical(member["records"]) == canonical([record for record in qualifying if record["read_id"] == member["read_id"] and record["footprint"] == member["footprint"]]), "Witness raw record binding differs")
        require(member["records"], "Empty witness assignment")
    maximum = matching_size_three(links)
    require(len(witness) == maximum, "Witness does not attain bounded coupled matching size")
    eligible = maximum == 3 and len(groups) >= 3
    require(summary["experimental_local_eligibility"] == eligible, "Experimental eligibility differs")
    require(summary["full_product_phase"] == "not_evaluated", "Unexpected full-product phase claim")
    return {"id": item["id"], "region_id": item["region_id"], **fields,
            "represented_orientations": [entry["strand"] for entry in item["representations"] if entry["represented"]],
            "canonical_sequence_groups": len(groups), "witness_size": maximum,
            "experimental_local_eligibility": eligible}


def validate_report(report, graph, classes):
    require(report["status"] == "complete" and report["experimental_only"] and not report["production_admission_changed"], "Not a completed standalone diagnostic")
    require(report["consumed_prefix"] == EXPECTED_PREFIX, "Wrong consumed prefix")
    records = report["records"]
    require(all(type(record["ordinal"]) is int and 1 <= record["ordinal"] <= EXPECTED_PREFIX["records"] for record in records), "Retained ordinal outside exact prefix")
    require(len({record["ordinal"] for record in records}) == len(records), "Duplicate retained ordinal")
    require([record["ordinal"] for record in records] == sorted(record["ordinal"] for record in records), "Ledger ordinals out of order")
    layer_counts = {strand: {layer: Counter() for layer in ("unfiltered", "q20")} for strand in ("+", "-")}
    visits = 0
    for record in records:
        sequence, quality = record["sequence"], record["quality"]
        require(set(record["histories"]) == {"+", "-"} and all(set(layers) == {"unfiltered", "q20"} for layers in record["histories"].values()), "Unexpected history strand/layer")
        require(len(sequence) == len(quality), "Read/quality length mismatch")
        require(record["sequence_sha256"] == sha(sequence.encode()), "Read sequence hash differs")
        require(record["canonical_sequence_sha256"] == sha(min(sequence, rc(sequence)).encode()), "Canonical read sequence hash differs")
        for strand, oriented, scores in (("+", sequence, quality), ("-", rc(sequence), quality[::-1])):
            for layer in ("unfiltered", "q20"):
                expected = history_oracle(oriented, scores, graph, layer == "q20")
                require(record["histories"][strand][layer] == expected, f"History differs at ordinal {record['ordinal']} {strand}/{layer}")
                counter = layer_counts[strand][layer]
                counter.update(records_with_runs=bool(expected["runs"]), runs=len(expected["runs"]), break_intervals=len(expected["breaks"]))
                for interruption in expected["breaks"]:
                    for reason in interruption["reasons"]:
                        counter["broken_windows_" + reason] += interruption["end"] - interruption["start"]
                for run in expected["runs"]:
                    counter.update(edge_visits=len(run["edges"]), branch_events=len(run["branch_events"]),
                        cyclic_node_visits=len(run["cyclic_node_visits"]), repeated_edge_occurrences=run["repeated_edge_occurrences"],
                        runs_touching_cyclic_scc=bool(run["cyclic_scc_ids"]))
                    visits += len(run["edges"])
        require(canonical(record["observations"]) == canonical(placements(sequence, quality, classes, record["histories"])), f"Class placements differ at ordinal {record['ordinal']}")
    for strand in layer_counts:
        for layer in layer_counts[strand]:
            require(layer_counts[strand][layer] == Counter(report["layer_counts"][strand][layer]), "Layer summary counters differ")
    require(report["counts"]["retained_records"] == len(records) and report["counts"]["retained_edge_visits"] == visits, "Retention counters differ")
    limits = report["limits"]
    require(len(records) <= limits["MAX_HISTORY_RECORDS"] and visits <= limits["MAX_TOTAL_EDGE_VISITS"], "Evidence bound exceeded")
    require(all(len(record["sequence"]) <= limits["MAX_READ_LENGTH"] and len(record["observations"]) <= limits["MAX_PLACEMENTS_PER_READ"] for record in records), "Per-read bound exceeded")
    summaries = {item["id"]: item for item in report["classes"]}
    require(len(summaries) == len(report["classes"]), "Duplicate class summary")
    require(set(summaries) == {item["id"] for item in classes}, "Class summary identities differ")
    return {"retained_records_validated": len(records), "edge_visits_validated": visits,
            "classes": [audit_class_summary(item, summaries[item["id"]], records) for item in classes],
            "layer_counts": layer_counts, "limits": limits}


def prefix_binding(path, reports, graphs, classes):
    """One independent prefix pass checks ledger completeness and raw binding."""
    ledger = [{record["ordinal"]: record for record in report["records"]} for report in reports]
    k = graphs[0]["k"]
    require(all(graph["k"] == k for graph in graphs), "Graph k differs")
    edge_mask = defaultdict(int)
    for index, graph in enumerate(graphs):
        for edge in graph["edges"]:
            edge_mask[edge] |= 1 << index
            edge_mask[rc(edge)] |= 1 << index
    all_graphs = (1 << len(graphs)) - 1
    masks = dict(edge_mask)
    totals = [Counter() for _ in reports]
    class_counts = {item["id"]: Counter() for item in classes}
    hasher, bases, size = hashlib.sha256(), 0, 0
    before = identity(path)
    with Path(path).open("rb") as stream:
        for ordinal in range(1, EXPECTED_PREFIX["records"] + 1):
            lines = [stream.readline() for _ in range(4)]
            require(all(line.endswith(b"\n") for line in lines), f"Malformed prefix at ordinal {ordinal}")
            raw = b"".join(lines)
            header, sequence, separator, quality = [line.rstrip(b"\r\n").decode("ascii") for line in lines]
            require(header.startswith("@") and separator.startswith("+") and len(sequence) == len(quality), "Malformed raw FASTQ")
            hasher.update(raw)
            size += len(raw)
            bases += len(sequence)
            mask = 0
            for start in range(len(sequence) - k + 1):
                mask |= masks.get(sequence[start:start + k], 0)
                if mask == all_graphs:
                    break
            reverse, class_hit = rc(sequence), False
            for item in classes:
                literal = item["sequence"]
                found, qualified = False, False
                for oriented, scores in ((sequence, quality), (reverse, quality[::-1])):
                    start = oriented.find(literal)
                    while start >= 0:
                        found = True
                        qualified |= min(map(ord, scores[start:start + len(literal)])) >= 53
                        start = oriented.find(literal, start + 1)
                class_counts[item["id"]].update(exact_records=found, q20_exact_records=qualified)
                class_hit |= found
            for index in range(len(reports)):
                graph_hit = bool(mask & (1 << index))
                required = graph_hit or class_hit
                require((ordinal in ledger[index]) == required, f"Incomplete or surplus graph {index} ledger at ordinal {ordinal}")
                totals[index]["records"] += 1
                if not required:
                    totals[index]["records_without_graph_edge"] += 1
                    continue
                totals[index]["retained_records"] += 1
                if not graph_hit:
                    totals[index]["class_only_records_without_graph_edge"] += 1
                record = ledger[index][ordinal]
                require(record["raw_fastq"].encode() == raw and record["sequence"] == sequence and record["quality"] == quality and record["read_id"] == header[1:].split()[0], f"Raw prefix binding differs at ordinal {ordinal}")
    observed = {"records": EXPECTED_PREFIX["records"], "bases": bases, "size_bytes": size, "sha256": hasher.hexdigest()}
    require(observed == EXPECTED_PREFIX, "Independent prefix identity differs")
    require(identity(path) == before, "Raw input changed during independent audit")
    for index, report in enumerate(reports):
        for key in ("records", "records_without_graph_edge", "class_only_records_without_graph_edge"):
            require(totals[index][key] == report["counts"].get(key, 0), f"Prefix retention count differs: {key}")
        require(totals[index]["retained_records"] == len(report["records"]) == report["counts"]["retained_records"], "Independently required ledger size differs")
        for summary in report["classes"]:
            require(all(summary[key] == class_counts[summary["id"]][key] for key in ("exact_records", "q20_exact_records")), "Literal class count differs from full independent prefix pass")
    return {"consumed_prefix": observed, "full_staged_file": before,
            "all_prefix_records_checked_for_ledger_completeness": EXPECTED_PREFIX["records"],
            "literal_class_counts": class_counts}


def archived_exact_counts(classes):
    archive = ROOT / "benchmarks/benchmark_results/v3.2-its-scc-audit-20260913/evidence.tar.gz"
    require(identity(archive)["sha256"] == ARCHIVE_SHA, "Archived SCC evidence changed")
    inventory = {item["path"]: item for item in read_json(EVIDENCE / "inputs/archive_contents.json")["members"]}
    with tarfile.open(archive, "r:gz") as stream:
        raw = stream.extractfile("read_summary.json").read()
    require({"size_bytes": len(raw), "sha256": sha(raw)} == basic(inventory["read_summary.json"]), "Archived read summary differs from inventory")
    summary = json.loads(raw)
    require(summary["analysis_sha256"] == identity(EVIDENCE / "inputs/scc_analysis.json")["sha256"], "Archived summary/analysis mismatch")
    rows = {row["class_id"]: row for row in summary["rows"] if row["max_substitutions"] == 0 and row["region"] in ("region_1", "region_2")}
    require(set(rows) == {item["id"] for item in classes}, "Archived exact class identities differ")
    for item in classes:
        require(rows[item["id"]]["window_sha256"] == sha(item["sequence"].encode()), "Archived class sequence mismatch")
    return {"source": "archived read_summary.json, max_substitutions=0", "identity": basic(inventory["read_summary.json"]),
            "q20_exact_records": {key: row["records"] for key, row in rows.items()}}


def audit(freeze_path):
    freeze_identity = identity(freeze_path)
    plan = read_json(freeze_path)
    run_dir = Path(plan["run_dir"])
    artifacts = {Path(freeze_path): freeze_identity,
                 run_dir / "execution.json": identity(run_dir / "execution.json"),
                 run_dir / "started.json": identity(run_dir / "started.json")}
    require(plan["input_identity"] == basic(EXPECTED_PREFIX), "Frozen input identity differs")
    require(identity(plan["fastq"]) == plan["input_identity"], "Staged input differs from freeze")
    for relative, expected in plan["sources"].items():
        require(identity(ROOT / relative) == expected, f"Frozen live source differs: {relative}")
        require(identity(Path(plan["source_snapshot_root"]) / relative) == expected, f"Frozen snapshot differs: {relative}")
    execution = read_json(run_dir / "execution.json")
    started = read_json(run_dir / "started.json")
    require(execution["status"] == "complete" and execution["freeze"] == freeze_identity and started["freeze"] == freeze_identity, "Execution/freeze receipt mismatch")
    require(execution["commands_planned"] == [item["id"] for item in plan["commands"]], "Planned command IDs differ")
    require(len(execution["commands"]) == len(plan["commands"]) == 2, "Incomplete graph command set")
    require([item["id"] for item in plan["commands"]] == ["histories_pre", "histories_post"], "Expected one pre and one post graph command")
    manifest = read_json(EVIDENCE / "inputs/its_manifest.json")
    class_templates = [{**item, "region_id": region["id"]} for region in manifest["regions"] for item in region["classes"]]
    limits = frozen_limits(ROOT / "scripts/audit_read_histories.py")
    reports, graphs, results = [], [], []
    for command, receipt in zip(plan["commands"], execution["commands"]):
        argv = command["argv"]
        options = {argv[index]: argv[index + 1] for index in range(2, len(argv), 2)}
        phase = command["id"].removeprefix("histories_")
        require(Path(argv[1]) == ROOT / "scripts/audit_read_histories.py" and options["--fastq"] == plan["fastq"], "Unexpected prototype or read input")
        require(Path(options["--graph"]) == EVIDENCE / f"inputs/graph_{phase}.json" and identity(options["--graph"])["sha256"] == GRAPH_SHA[phase], "Unexpected graph input identity")
        require(Path(options["--manifest"]) == EVIDENCE / "inputs/its_manifest.json" and Path(options["--expected-prefix"]) == EVIDENCE / "inputs/expected_prefix.json", "Unexpected class/prefix inputs")
        require(receipt["id"] == command["id"] and receipt["status"] == "complete" and receipt["exit_code"] == 0 and receipt["source_input_bookends_match"], "Command failed or bookends differ")
        require(receipt["argv"] == ["/usr/bin/time", "-v", "-o", str(run_dir / f"{command['id']}.time.txt"), *argv] and receipt["freeze"] == freeze_identity, "Actual command/freeze differs")
        require(receipt["input_before"] == receipt["input_after"] == plan["input_identity"], "Wrapper input bookends differ")
        receipt_path = run_dir / f"{command['id']}.receipt.json"
        artifacts[receipt_path] = identity(receipt_path)
        require(read_json(receipt_path) == receipt, "Inline/standalone command receipts differ")
        for filename, expected in receipt["outputs"].items():
            require(identity(run_dir / filename) == expected, f"Output hash differs: {filename}")
            artifacts[run_dir / filename] = expected
        report_path = Path(options["--output"])
        require(report_path == run_dir / f"{command['id']}.json" and report_path.name in receipt["outputs"], "Prototype report absent from output hash receipt")
        require(report_path.stat().st_size <= limits["MAX_OUTPUT_BYTES"], "Serialized output exceeds frozen bound")
        report = read_json(report_path)
        require(report["sources_before"] == report["sources_after"], "Prototype source bookends differ")
        paths = {"graph": options["--graph"], "manifest": options["--manifest"], "expected_prefix": options["--expected-prefix"], "scanner": argv[1]}
        for key, path in paths.items():
            require(basic(report["sources_before"][key]) == identity(path), f"Prototype {key} identity differs")
        require(report["staged_fastq_before"] == report["staged_fastq_after"] and basic(report["staged_fastq_before"]) == plan["input_identity"], "Prototype FASTQ bookends differ")
        require(report["limits"] == limits and report_path.stat().st_size <= limits["MAX_OUTPUT_BYTES"], "Frozen limits/output bound differs")
        graph_document = read_json(options["--graph"])
        graph = graph_model(graph_document)
        require(report["target"] == graph_document["target"] and report["graph_phase"] == graph_document["graph"]["phase"], "Graph target/phase differs")
        require(len(graph["nodes"]) <= limits["MAX_NODES"] and len(graph["edges"]) <= limits["MAX_EDGES"], "Graph bounds exceeded")
        classes = [{**item, "representations": representations(item["sequence"], graph)} for item in class_templates]
        results.append({"command_id": command["id"], "output": identity(report_path), **validate_report(report, graph, classes)})
        reports.append(report)
        graphs.append(graph)
    binding = prefix_binding(plan["fastq"], reports, graphs, class_templates)
    archived = archived_exact_counts(class_templates)
    for result in results:
        for item in result["classes"]:
            item["archived_q20_exact_records"] = archived["q20_exact_records"][item["id"]]
            item["archived_q20_compatible"] = item["q20_exact_records"] == item["archived_q20_exact_records"]
            require(item["archived_q20_compatible"], "Archived Q20 exact control mismatch")
    require(identity(freeze_path) == freeze_identity, "Freeze changed during independent review")
    for path, expected in artifacts.items():
        require(identity(path) == expected, f"Reviewed artifact changed during audit: {path}")
    for relative, expected in plan["sources"].items():
        require(identity(ROOT / relative) == expected, f"Source changed during independent review: {relative}")
        require(identity(Path(plan["source_snapshot_root"]) / relative) == expected, f"Snapshot changed during independent review: {relative}")
    return {"status": "pass", "freeze": freeze_identity, "execution": artifacts[run_dir / "execution.json"],
            "reviewed_artifacts": {str(path.relative_to(ROOT)): expected for path, expected in artifacts.items()},
            "prefix_binding": binding, "archived_comparison": archived, "graphs": results,
            "scope": "All retained histories and all raw prefix ordinals independently audited; no prototype helpers imported",
            "limitations": ["Exact local compatibility is not biological absence, repeat copy number, independent molecules, or full-product phase", "Archived counts are Q20 exact comparisons; unfiltered exact counts are reported separately", "No production admission or reference change"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "provenance_ledger_review.json")
    args = parser.parse_args()
    failure = args.output.with_name(args.output.name + ".failure.json")
    require(not any(path.exists() or path.is_symlink() for path in (args.output, failure)), "Review output/failure exists")
    source = identity(Path(__file__))
    started = time.monotonic()
    try:
        result = audit(args.freeze)
        require(identity(Path(__file__)) == source, "Checker changed during review")
        result.update(checker=source, reviewed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), elapsed_seconds=time.monotonic() - started)
        publish(args.output, result)
    except Exception as error:
        publish(failure, {"status": "failed", "checker": source, "error_type": type(error).__name__, "error": str(error), "elapsed_seconds": time.monotonic() - started})
        raise
    print(json.dumps({"status": "pass", "output": str(args.output), "graphs": len(result["graphs"]), "elapsed_seconds": result["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
