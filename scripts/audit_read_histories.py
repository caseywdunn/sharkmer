#!/usr/bin/env python3
"""Experimental exact read histories; this tool does not admit assembly paths."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import tempfile


SCHEMA_VERSION = 1
MIN_QUALITY = 20
FLANK_BASES = 21
MAX_GRAPH_BYTES = 64 * 1024 * 1024
MAX_MANIFEST_BYTES = 1024 * 1024
MAX_NODES = 100000
MAX_EDGES = 200000
MAX_CLASSES = 64
MAX_REGIONS = 16
MAX_READ_LENGTH = 1000
MAX_LINE_BYTES = 4096
MAX_HISTORY_RECORDS = 100000
MAX_TOTAL_EDGE_VISITS = 10000000
MAX_PLACEMENTS_PER_READ = 4096
MAX_OUTPUT_BYTES = 1024 * 1024 * 1024
MAX_PREFIX_RECORDS = 1000000
BASES = frozenset(b"ACGT")
COMPLEMENT = bytes.maketrans(b"ACGTRYSWKMBDHVN", b"TGCAYRSWMKVHDBN")


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def positions(sequence, pattern):
    start = sequence.find(pattern)
    while start >= 0:
        yield start
        start = sequence.find(pattern, start + 1)


def bounded_json(path, maximum):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("JSON input must be a regular nonsymlink file")
    with path.open("rb") as stream:
        raw = stream.read(maximum + 1)
    if len(raw) > maximum:
        raise ValueError("JSON input byte cap exceeded")
    return json.loads(raw)


def file_identity(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("Input must be a regular nonsymlink file")
    before = path.stat()
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(chunk)
    after = path.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(before, field) != getattr(after, field) for field in fields):
        raise ValueError("Input changed during hashing")
    return {"sha256": hasher.hexdigest(), "size_bytes": after.st_size,
            "device": after.st_dev, "inode": after.st_ino,
            "mtime_ns": after.st_mtime_ns, "ctime_ns": after.st_ctime_ns}


def literal(value, length=None):
    if not isinstance(value, str):
        raise ValueError("Expected literal DNA string")
    sequence = value.encode("ascii")
    if not sequence or set(sequence) - BASES or (length is not None and len(sequence) != length):
        raise ValueError("Invalid literal DNA or length")
    return sequence


def identifier(value):
    if not isinstance(value, str) or not value or len(value.encode("utf8")) > 256:
        raise ValueError("Invalid identifier")
    return value


def cyclic_components(nodes, outgoing, incoming):
    """Iterative Kosaraju, independently check the graph's SCC annotations."""
    seen, order = set(), []
    for root in nodes:
        if root in seen:
            continue
        seen.add(root)
        stack = [(root, iter(outgoing[root]))]
        while stack:
            node, neighbors = stack[-1]
            target = next(neighbors, None)
            if target is None:
                order.append(node)
                stack.pop()
            elif target not in seen:
                seen.add(target)
                stack.append((target, iter(outgoing[target])))
    seen, result = set(), []
    for root in reversed(order):
        if root in seen:
            continue
        seen.add(root)
        stack, component = [root], set()
        while stack:
            node = stack.pop()
            component.add(node)
            for source in incoming[node]:
                if source not in seen:
                    seen.add(source)
                    stack.append(source)
        if len(component) > 1 or root in outgoing[root]:
            result.append(component)
    return result


def prepare_graph(document):
    if document.get("schema_version") != 1:
        raise ValueError("Expected graph schema version 1")
    k = document["target"]["kmer_length"]
    if type(k) is not int or not 2 <= k <= FLANK_BASES:
        raise ValueError("Unsupported graph k")
    source = document["graph"]
    if not isinstance(source["nodes"], list) or not isinstance(source["edges"], list):
        raise ValueError("Graph nodes and edges must be lists")
    if len(source["nodes"]) > MAX_NODES or len(source["edges"]) > MAX_EDGES:
        raise ValueError("Graph node/edge cap exceeded")
    nodes, edges = {}, {}
    for node in source["nodes"]:
        sequence = literal(node["oriented_sequence"], k - 1)
        if node["node_id"] != sequence.decode() or sequence in nodes:
            raise ValueError("Duplicate or inconsistent oriented node")
        nodes[sequence] = node
    outgoing = {node: [] for node in nodes}
    incoming = {node: [] for node in nodes}
    for edge in source["edges"]:
        sequence = literal(edge["oriented_sequence"], k)
        left, right = sequence[:-1], sequence[1:]
        if (edge["edge_id"] != sequence.decode() or sequence in edges or
                edge["source_node_id"] != left.decode() or edge["target_node_id"] != right.decode() or
                left not in nodes or right not in nodes):
            raise ValueError("Duplicate or inconsistent oriented edge/endpoints")
        if type(edge["count"]) is not int or edge["count"] < 0:
            raise ValueError("Invalid edge count")
        for side, node in (("source", left), ("target", right)):
            if edge.get(side + "_cyclic_scc_id") != nodes[node].get("cyclic_scc_id"):
                raise ValueError("Edge/node SCC annotation differs")
        edges[sequence] = edge
        outgoing[left].append(right)
        incoming[right].append(left)
    if source["node_count"] != len(nodes) or source["edge_count"] != len(edges):
        raise ValueError("Graph counts differ")
    declared = defaultdict(set)
    for sequence, node in nodes.items():
        if node["incoming_degree"] != len(incoming[sequence]) or node["outgoing_degree"] != len(outgoing[sequence]):
            raise ValueError("Graph degrees differ")
        scc = node.get("cyclic_scc_id")
        if scc is not None:
            if not isinstance(scc, (str, int)) or isinstance(scc, bool):
                raise ValueError("Invalid SCC identifier")
            declared[scc].add(sequence)
    actual = cyclic_components(nodes, outgoing, incoming)
    if {frozenset(group) for group in declared.values()} != {frozenset(group) for group in actual}:
        raise ValueError("Graph cyclic SCC labels differ from topology")
    return {"k": k, "nodes": nodes, "edges": edges, "outgoing": outgoing,
            "incoming": incoming, "phase": source.get("phase"), "target": document["target"]}


def prepare_manifest(document, graph):
    if document.get("schema_version") != 1 or not isinstance(document.get("regions"), list):
        raise ValueError("Expected class manifest schema version 1")
    if not 1 <= len(document["regions"]) <= MAX_REGIONS:
        raise ValueError("Region count outside bound")
    classes, seen_regions, seen_classes = [], set(), set()
    k = graph["k"]
    for region in document["regions"]:
        region_id = identifier(region["id"])
        if region_id in seen_regions or not region["classes"]:
            raise ValueError("Duplicate or empty region")
        seen_regions.add(region_id)
        for item in region["classes"]:
            class_id = identifier(item["id"])
            sequence = literal(item["sequence"])
            if class_id in seen_classes or not 2 * FLANK_BASES <= len(sequence) <= MAX_READ_LENGTH:
                raise ValueError("Duplicate class or incomplete/bounded flanks")
            if digest(sequence) != item["sequence_sha256"]:
                raise ValueError("Class sequence checksum differs")
            seen_classes.add(class_id)
            representations = []
            for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
                visits = [oriented[start:start + k] for start in range(len(oriented) - k + 1)]
                missing = [index for index, edge in enumerate(visits) if edge not in graph["edges"]]
                anchors = [oriented[:k - 1], oriented[-(k - 1):]]
                noncyclic = all(node in graph["nodes"] and graph["nodes"][node].get("cyclic_scc_id") is None for node in anchors)
                representations.append({"strand": strand, "represented": not missing,
                                        "missing_edge_starts": missing,
                                        "anchor_nodes": [node.decode() for node in anchors],
                                        "noncyclic_anchors": noncyclic})
            classes.append({"id": class_id, "region_id": region_id, "sequence": sequence,
                            "kind": item.get("kind"), "representations": representations})
            if len(classes) > MAX_CLASSES:
                raise ValueError("Class count cap exceeded")
    return classes


def map_oriented(sequence, quality, graph, quality_min=None):
    """Keep every exact edge occurrence and every window-level break reason."""
    if len(sequence) != len(quality) or len(sequence) > MAX_READ_LENGTH:
        raise ValueError("Read sequence/quality length mismatch or cap exceeded")
    k = graph["k"]
    runs, breaks, visits, start = [], [], [], None

    def flush():
        nonlocal visits, start
        if not visits:
            return
        spelled = visits[0] + bytes(edge[-1] for edge in visits[1:])
        end = start + len(visits) + k - 1
        if spelled != sequence[start:end]:
            raise ValueError("Run spelling differs from read interval")
        branches = []
        for offset, (left, right) in enumerate(zip(visits, visits[1:])):
            node = left[1:]
            if node != right[:-1]:
                raise ValueError("Run endpoints are disconnected")
            if len(graph["incoming"][node]) > 1 or len(graph["outgoing"][node]) > 1:
                branches.append({"edge_offset": offset, "read_kmer_start": start + offset,
                                 "incoming_edge": left.decode(), "outgoing_edge": right.decode()})
        counts = Counter(visits)
        cyclic_visits = []
        for offset, node in enumerate([visits[0][:-1]] + [edge[1:] for edge in visits]):
            scc = graph["nodes"][node].get("cyclic_scc_id")
            if scc is not None:
                cyclic_visits.append({"node_offset": offset, "read_node_start": start + offset,
                                      "node_id": node.decode(), "cyclic_scc_id": scc})
        sccs = {visit["cyclic_scc_id"] for visit in cyclic_visits}
        runs.append({"start": start, "end": end, "sequence": spelled.decode(),
                     "edges": [edge.decode() for edge in visits], "branch_events": branches,
                     "cyclic_node_visits": cyclic_visits,
                     "repeated_edge_occurrences": sum(count - 1 for count in counts.values()),
                     "cyclic_scc_ids": sorted((value for value in sccs if value is not None), key=str)})
        visits, start = [], None

    invalid = [0]
    low = [0]
    for base, score in zip(sequence, quality):
        invalid.append(invalid[-1] + (base not in BASES))
        low.append(low[-1] + (quality_min is not None and score < quality_min + 33))
    for offset in range(max(0, len(sequence) - k + 1)):
        edge = sequence[offset:offset + k]
        reasons = []
        if invalid[offset + k] != invalid[offset]:
            reasons.append("invalid_base")
        if low[offset + k] != low[offset]:
            reasons.append("low_quality")
        if edge not in graph["edges"]:
            reasons.append("missing_graph_kmer")
        if reasons:
            flush()
            if breaks and breaks[-1]["end"] == offset and breaks[-1]["reasons"] == reasons:
                breaks[-1]["end"] = offset + 1
            else:
                breaks.append({"start": offset, "end": offset + 1, "reasons": reasons})
            continue
        if visits and visits[-1][1:] != edge[:-1]:
            raise ValueError("Consecutive graph kmers have inconsistent endpoints")
        if start is None:
            start = offset
        visits.append(edge)
    flush()
    return {"runs": runs, "breaks": breaks}


def read_histories(sequence, quality, graph):
    result = {}
    for strand, oriented, scores in (("+", sequence, quality), ("-", reverse_complement(sequence), quality[::-1])):
        result[strand] = {"unfiltered": map_oriented(oriented, scores, graph),
                          "q20": map_oriented(oriented, scores, graph, MIN_QUALITY)}
    return result


def observations_for_read(sequence, quality, classes, histories):
    """Placement ambiguity is computed before applying quality or graph gates."""
    observations, per_region = [], defaultdict(list)
    for item in classes:
        for strand, oriented, scores in (("+", sequence, quality), ("-", reverse_complement(sequence), quality[::-1])):
            for start in positions(oriented, item["sequence"]):
                end = start + len(item["sequence"])
                physical = [start, end] if strand == "+" else [len(sequence) - end, len(sequence) - start]
                paths = {"unfiltered": [], "q20": []}
                anchors_ok = False
                for representation in item["representations"]:
                    graph_strand = strand if representation["strand"] == "+" else ("-" if strand == "+" else "+")
                    interval = [start, end] if representation["strand"] == "+" else [len(sequence) - end, len(sequence) - start]
                    for layer in paths:
                        for run_index, run in enumerate(histories[graph_strand][layer]["runs"]):
                            if run["start"] <= interval[0] and run["end"] >= interval[1]:
                                paths[layer].append({"graph_strand": graph_strand, "run_index": run_index,
                                                     "class_in_graph_strand": representation["strand"],
                                                     "interval": interval})
                                if layer == "q20" and representation["noncyclic_anchors"]:
                                    anchors_ok = True
                observation = {"region_id": item["region_id"], "class_id": item["id"],
                               "class_read_strand": strand, "class_oriented_start": start,
                               "physical_interval": physical, "footprint": [-start, len(sequence) - start],
                               "minimum_quality": min(scores[start:end]) - 33,
                               "paths": paths, "noncyclic_q20_anchors": anchors_ok}
                observations.append(observation)
                per_region[item["region_id"]].append(observation)
                if len(observations) > MAX_PLACEMENTS_PER_READ:
                    raise ValueError("Read class-placement cap exceeded")
    for observation in observations:
        observation["unique_unfiltered_placement"] = len(per_region[observation["region_id"]]) == 1
        observation["qualifying_local_witness"] = (observation["unique_unfiltered_placement"] and
                                                   observation["minimum_quality"] >= MIN_QUALITY and
                                                   observation["noncyclic_q20_anchors"] and
                                                   bool(observation["paths"]["q20"]))
    return observations


def corroboration_witness(links):
    """Coupled distinct original-ID/footprint matching, capped at three owners."""
    owners = {}

    def assign(read_id, visited):
        for footprint in sorted(links[read_id]):
            if footprint in visited:
                continue
            visited.add(footprint)
            if footprint not in owners or assign(owners[footprint], visited):
                owners[footprint] = read_id
                return True
        return False

    for read_id in sorted(links):
        assign(read_id, set())
        if len(owners) >= 3:
            break
    return [{"read_id": read_id, "footprint": list(footprint)} for footprint, read_id in sorted(owners.items())]


def summarize_classes(classes, ledger):
    result = []
    for item in classes:
        links, records, groups = defaultdict(set), [], defaultdict(set)
        exact, q20, ambiguous = set(), set(), set()
        for record in ledger:
            for observation in record["observations"]:
                if observation["class_id"] != item["id"]:
                    continue
                ordinal = record["ordinal"]
                exact.add(ordinal)
                if observation["minimum_quality"] >= MIN_QUALITY:
                    q20.add(ordinal)
                if not observation["unique_unfiltered_placement"]:
                    ambiguous.add(ordinal)
                if observation["qualifying_local_witness"]:
                    footprint = tuple(observation["footprint"])
                    links[record["read_id"]].add(footprint)
                    groups[record["canonical_sequence_sha256"]].add(footprint)
                    records.append({"ordinal": ordinal, "read_id": record["read_id"],
                                    "canonical_sequence_sha256": record["canonical_sequence_sha256"],
                                    "footprint": list(footprint), "class_read_strand": observation["class_read_strand"]})
        if any(len(footprints) != 1 for footprints in groups.values()):
            raise ValueError("Identical/RC sequence produced different unique class footprints")
        witness = corroboration_witness(links)
        for member in witness:
            member["records"] = [record for record in records if record["read_id"] == member["read_id"] and record["footprint"] == member["footprint"]]
        witness_groups = {record["canonical_sequence_sha256"] for member in witness for record in member["records"]}
        eligible = len(witness) >= 3 and len(groups) >= 3 and len(witness_groups) >= 3
        result.append({"id": item["id"], "region_id": item["region_id"], "kind": item["kind"],
                       "sequence": item["sequence"].decode(), "representations": item["representations"],
                       "exact_records": len(exact), "q20_exact_records": len(q20),
                       "ambiguous_exact_records": len(ambiguous), "qualifying_records": len(records),
                       "distinct_read_ids": len(links), "distinct_footprints": len(set().union(*links.values())) if links else 0,
                       "canonical_sequence_groups": [{"sha256": key, "footprints": [list(value) for value in sorted(values)]} for key, values in sorted(groups.items())],
                       "witness": witness, "experimental_local_eligibility": eligible,
                       "full_product_phase": "not_evaluated"})
    return result


def expected_fields(expected):
    required = ("records", "bases", "size_bytes", "sha256")
    # The preceding audits use bytes for the consumed-prefix byte count.
    expected = dict(expected)
    if "size_bytes" not in expected and "bytes" in expected:
        expected["size_bytes"] = expected["bytes"]
    if any(key not in expected for key in required):
        raise ValueError("Expected prefix needs records, bases, size_bytes and sha256")
    if (type(expected["records"]) is not int or not 1 <= expected["records"] <= MAX_PREFIX_RECORDS or
            any(type(expected[key]) is not int or expected[key] <= 0 for key in ("bases", "size_bytes")) or
            not isinstance(expected["sha256"], str) or len(expected["sha256"]) != 64):
        raise ValueError("Invalid expected prefix identity")
    return {key: expected[key] for key in required}


def scan_fastq(path, expected, graph, classes):
    expected = expected_fields(expected)
    path = Path(path)
    before = file_identity(path)
    hasher, size, bases = hashlib.sha256(), 0, 0
    ledger, counters, edge_visits = [], Counter(), 0
    layer_counts = {strand: {layer: Counter() for layer in ("unfiltered", "q20")} for strand in ("+", "-")}
    k = graph["k"]
    with path.open("rb") as stream:
        for ordinal in range(1, expected["records"] + 1):
            lines = [stream.readline(MAX_LINE_BYTES + 1) for _ in range(4)]
            if any(not line.endswith(b"\n") or len(line) > MAX_LINE_BYTES for line in lines):
                raise ValueError("Truncated, unterminated or overlong FASTQ record")
            header, sequence, separator, quality = [line.rstrip(b"\r\n") for line in lines]
            if (not header.startswith(b"@") or not header[1:].split() or not separator.startswith(b"+") or
                    not sequence or len(sequence) != len(quality) or len(sequence) > MAX_READ_LENGTH or
                    any(not 65 <= base <= 90 for base in sequence) or
                    any(score < 33 or score > 126 for score in quality)):
                raise ValueError("Malformed or overlong FASTQ record")
            read_id = identifier(header[1:].split()[0].decode("ascii"))
            raw = b"".join(lines)
            hasher.update(raw)
            size += len(raw)
            bases += len(sequence)
            counters["records"] += 1
            # Preserve literal class placements even if their graph is absent.
            reverse = reverse_complement(sequence)
            graph_hit = any(oriented[offset:offset + k] in graph["edges"] for oriented in (sequence, reverse) for offset in range(max(0, len(sequence) - k + 1)))
            class_hit = any(item["sequence"] in oriented for oriented in (sequence, reverse) for item in classes)
            if not graph_hit and not class_hit:
                counters["records_without_graph_edge"] += 1
                continue
            if not graph_hit:
                counters["class_only_records_without_graph_edge"] += 1
            if len(ledger) >= MAX_HISTORY_RECORDS:
                raise ValueError("Retained-history record cap exceeded")
            histories = read_histories(sequence, quality, graph)
            observations = observations_for_read(sequence, quality, classes, histories)
            for strand, layers in histories.items():
                for layer, history in layers.items():
                    counts = layer_counts[strand][layer]
                    counts["records_with_runs"] += bool(history["runs"])
                    counts["runs"] += len(history["runs"])
                    counts["break_intervals"] += len(history["breaks"])
                    for interruption in history["breaks"]:
                        for reason in interruption["reasons"]:
                            counts["broken_windows_" + reason] += interruption["end"] - interruption["start"]
                    for run in history["runs"]:
                        counts["edge_visits"] += len(run["edges"])
                        counts["branch_events"] += len(run["branch_events"])
                        counts["cyclic_node_visits"] += len(run["cyclic_node_visits"])
                        counts["repeated_edge_occurrences"] += run["repeated_edge_occurrences"]
                        counts["runs_touching_cyclic_scc"] += bool(run["cyclic_scc_ids"])
                        edge_visits += len(run["edges"])
            if edge_visits > MAX_TOTAL_EDGE_VISITS:
                raise ValueError("Total retained edge-visit cap exceeded")
            ledger.append({"ordinal": ordinal, "read_id": read_id, "raw_fastq": raw.decode("ascii"),
                           "sequence": sequence.decode(), "quality": quality.decode(),
                           "sequence_sha256": digest(sequence),
                           "canonical_sequence_sha256": digest(min(sequence, reverse)),
                           "histories": histories, "observations": observations})
    consumed = {"records": counters["records"], "bases": bases, "size_bytes": size, "sha256": hasher.hexdigest()}
    if consumed != expected:
        raise ValueError("Consumed FASTQ prefix identity differs")
    after = file_identity(path)
    if after != before:
        raise ValueError("FASTQ changed during scan")
    return {"schema_version": SCHEMA_VERSION, "status": "complete", "experimental_only": True,
            "production_admission_changed": False, "consumed_prefix": consumed,
            "staged_fastq_before": before, "staged_fastq_after": after,
            "target": graph["target"], "graph_phase": graph["phase"],
            "counts": {**counters, "retained_records": len(ledger), "retained_edge_visits": edge_visits},
            "layer_counts": layer_counts, "classes": summarize_classes(classes, ledger),
            "layer_count_scope": "retained graph-edge or literal-class matching records",
            "break_coordinate_system": "oriented read kmer starts, zero-based half-open",
            "records": ledger,
            "limits": {key: value for key, value in globals().items() if key.startswith("MAX_")},
            "limitations": ["Exact graph walks only; no gap bridging or alignment rescue",
                            "Local finite-class eligibility is not independent-molecule or copy-number proof",
                            "No concatenation between runs, records, mates, or regions; no amplicon output"]}


def publish_json(path, value):
    """Stream an fsynced temporary file, publish complete JSON without overwrite."""
    path = Path(path)
    descriptor, temporary = tempfile.mkstemp(prefix=".read-histories-", dir=path.parent)
    try:
        size = 0
        with os.fdopen(descriptor, "wb") as stream:
            for chunk in json.JSONEncoder(sort_keys=True, separators=(",", ":")).iterencode(value):
                raw = chunk.encode("ascii")
                size += len(raw)
                if size > MAX_OUTPUT_BYTES:
                    raise ValueError("Serialized output byte cap exceeded")
                stream.write(raw)
            if size + 1 > MAX_OUTPUT_BYTES:
                raise ValueError("Serialized output byte cap exceeded")
            stream.write(b"\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("graph", "manifest", "expected-prefix", "fastq", "output"):
        parser.add_argument("--" + option, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise ValueError("Output exists; no overwrite permitted")
    sources = {"graph": args.graph, "manifest": args.manifest,
               "expected_prefix": args.expected_prefix, "scanner": Path(__file__)}
    identities = {key: file_identity(path) for key, path in sources.items()}
    graph = prepare_graph(bounded_json(args.graph, MAX_GRAPH_BYTES))
    classes = prepare_manifest(bounded_json(args.manifest, MAX_MANIFEST_BYTES), graph)
    report = scan_fastq(args.fastq, bounded_json(args.expected_prefix, MAX_MANIFEST_BYTES), graph, classes)
    after = {key: file_identity(path) for key, path in sources.items()}
    if identities != after:
        raise ValueError("Graph/manifest/expected/source changed during scan")
    report["sources_before"], report["sources_after"] = identities, after
    publish_json(args.output, report)
    print(json.dumps({"status": "complete", "output": str(args.output), "counts": report["counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
