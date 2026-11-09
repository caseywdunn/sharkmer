#!/usr/bin/env python3
"""Independent synthetic checks. No biological FASTQ is opened by this script."""

import argparse
from collections import Counter, defaultdict
import copy
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import random
import sys
import tempfile
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[4]
SCANNER = ROOT / "scripts/audit_read_histories.py"
SPEC = importlib.util.spec_from_file_location("reviewed_histories", SCANNER)
PRODUCER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PRODUCER)
CHECKS = []
FIXTURES = []


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def identity(path):
    raw = path.read_bytes()
    return {"size_bytes": len(raw), "sha256": sha(raw)}


def rc(sequence):
    complements = dict(zip(b"ACGTRYSWKMBDHVN", b"TGCAYRSWMKVHDBN"))
    return bytes(complements.get(base, base) for base in reversed(sequence))


def require(condition, description):
    if not condition:
        raise AssertionError(description)
    CHECKS.append(description)


def independent_scc(document):
    """Reachability intersections; intentionally differs from producer Kosaraju."""
    nodes = {item["node_id"] for item in document["graph"]["nodes"]}
    forward = {node: set() for node in nodes}
    reverse = {node: set() for node in nodes}
    for edge in document["graph"]["edges"]:
        left, right = edge["source_node_id"], edge["target_node_id"]
        forward[left].add(right)
        reverse[right].add(left)

    def reachable(root, adjacency):
        seen, todo = {root}, [root]
        while todo:
            for neighbor in adjacency[todo.pop()]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    todo.append(neighbor)
        return seen

    groups, unassigned = [], set(nodes)
    while unassigned:
        root = min(unassigned)
        group = reachable(root, forward) & reachable(root, reverse)
        unassigned -= group
        if len(group) > 1 or root in forward[root]:
            groups.append(group)
    return groups


def graph_from_sequences(sequences, k=19):
    edges = sorted({sequence[start:start + k].decode() for sequence in sequences
                    for start in range(len(sequence) - k + 1)})
    nodes = sorted({part for edge in edges for part in (edge[:-1], edge[1:])})
    incoming = Counter(edge[1:] for edge in edges)
    outgoing = Counter(edge[:-1] for edge in edges)
    document = {"schema_version": 1, "target": {"kmer_length": k}, "graph": {
        "phase": "independent_synthetic", "node_count": len(nodes), "edge_count": len(edges),
        "nodes": [{"node_id": node, "oriented_sequence": node,
                   "incoming_degree": incoming[node], "outgoing_degree": outgoing[node],
                   "cyclic_scc_id": None} for node in nodes],
        "edges": [{"edge_id": edge, "oriented_sequence": edge,
                   "source_node_id": edge[:-1], "target_node_id": edge[1:], "count": 10,
                   "source_cyclic_scc_id": None, "target_cyclic_scc_id": None} for edge in edges]}}
    labels = {node: index for index, group in enumerate(independent_scc(document)) for node in group}
    for node in document["graph"]["nodes"]:
        node["cyclic_scc_id"] = labels.get(node["node_id"])
    for edge in document["graph"]["edges"]:
        edge["source_cyclic_scc_id"] = labels.get(edge["source_node_id"])
        edge["target_cyclic_scc_id"] = labels.get(edge["target_node_id"])
    return document


def manifest_from_sequences(sequences):
    return {"schema_version": 1, "regions": [{"id": "synthetic_region", "classes": [
        {"id": "class_" + str(index), "sequence": sequence.decode(), "sequence_sha256": sha(sequence),
         "kind": "synthetic_calibration_only"} for index, sequence in enumerate(sequences)]}]}


def independent_runs(sequence, quality, graph, q20):
    k = graph["k"]
    chunks, current = [], []
    for offset in range(len(sequence) - k + 1):
        edge = sequence[offset:offset + k]
        good = (all(base in b"ACGT" for base in edge) and edge in graph["edges"] and
                (not q20 or all(score >= 53 for score in quality[offset:offset + k])))
        if good:
            current.append((offset, edge))
        elif current:
            chunks.append(current)
            current = []
    if current:
        chunks.append(current)
    return [(chunk[0][0], chunk[-1][0] + k, [edge.decode() for _, edge in chunk]) for chunk in chunks]


def verify_histories(sequence, quality, graph, histories, label):
    for strand, oriented, scores in (("+", sequence, quality), ("-", rc(sequence), quality[::-1])):
        for layer in ("unfiltered", "q20"):
            actual = histories[strand][layer]["runs"]
            expected = independent_runs(oriented, scores, graph, layer == "q20")
            require([(run["start"], run["end"], run["edges"]) for run in actual] == expected,
                    label + ": naive window/run oracle " + strand + "/" + layer)
            for run in actual:
                spell = run["edges"][0] + "".join(edge[-1] for edge in run["edges"][1:])
                assert spell == run["sequence"] == oriented[run["start"]:run["end"]].decode()
                assert all(left[1:] == right[:-1] for left, right in zip(run["edges"], run["edges"][1:]))
                assert run["repeated_edge_occurrences"] == len(run["edges"]) - len(set(run["edges"]))
                expected_events = []
                for offset, (left, right) in enumerate(zip(run["edges"], run["edges"][1:])):
                    node = left[1:].encode()
                    if len(graph["incoming"][node]) > 1 or len(graph["outgoing"][node]) > 1:
                        expected_events.append({"edge_offset": offset, "read_kmer_start": run["start"] + offset,
                                                "incoming_edge": left, "outgoing_edge": right})
                assert run["branch_events"] == expected_events
                if layer == "q20":
                    assert min(scores[run["start"]:run["end"]]) >= 53


def record(sequence, quality, graph, classes, ordinal=1, read_id=None, label="fixture"):
    histories = PRODUCER.read_histories(sequence, quality, graph)
    verify_histories(sequence, quality, graph, histories, label)
    return {"ordinal": ordinal, "read_id": read_id or "synthetic_" + str(ordinal),
            "canonical_sequence_sha256": sha(min(sequence, rc(sequence))),
            "sequence": sequence.decode(), "quality": quality.decode(), "histories": histories,
            "observations": PRODUCER.observations_for_read(sequence, quality, classes, histories)}


def actual_class_controls(inputs):
    manifest = json.loads((inputs / "its_manifest.json").read_text())
    archived = {item["id"]: (item["sequence"], item["sequence_sha256"])
                for region in (1, 2)
                for item in json.loads((inputs / "originals" / ("region_" + str(region) + ".scanner_manifest.json")).read_text())["candidates"]}
    observed = {item["id"]: (item["sequence"], item["sequence_sha256"])
                for region in manifest["regions"] for item in region["classes"]}
    require(observed == archived, "All four class strings/checksums exactly equal original archived scanner manifests")
    pinned = {"pre": "0810b40673c4e4bdc7dcad0573bbe7f3f49e8a0b4542ff5789bd7ce2ec545bf9",
              "post": "b23e89f11d6ed653c7d8113a3f2c07506932049c8ed2397d4083a2f5fd5dfd40"}
    for phase in ("pre", "post"):
        graph_path = inputs / ("graph_" + phase + ".json")
        require(identity(graph_path)["sha256"] == pinned[phase], phase + ": graph equals preregistered archive-member identity")
        document = json.loads(graph_path.read_text())
        actual = independent_scc(document)
        declared = defaultdict(set)
        for node in document["graph"]["nodes"]:
            if node.get("cyclic_scc_id") is not None:
                declared[node["cyclic_scc_id"]].add(node["node_id"])
        require({frozenset(group) for group in actual} == {frozenset(group) for group in declared.values()},
                phase + ": independent reachability SCCs equal archived annotations")
        graph = PRODUCER.prepare_graph(document)
        classes = PRODUCER.prepare_manifest(manifest, graph)
        require(len(classes) == 4, phase + ": exactly four frozen local classes")
        for item in classes:
            sequence = item["sequence"]
            representations = []
            for strand, candidate in (("+", sequence), ("-", rc(sequence))):
                edges = [candidate[start:start + graph["k"]] for start in range(len(candidate) - graph["k"] + 1)]
                if all(edge in graph["edges"] for edge in edges):
                    representations.append(strand)
                    anchors = {candidate[:graph["k"] - 1].decode(), candidate[-(graph["k"] - 1):].decode()}
                    require(all(not anchors.intersection(group) for group in actual),
                            phase + "/" + item["id"] + ": represented anchors lie outside independent SCCs")
            require(bool(representations), phase + "/" + item["id"] + ": actual literal class is represented")
            require(representations == [rep["strand"] for rep in item["representations"] if rep["represented"]],
                    phase + "/" + item["id"] + ": orientation agrees with literal kmer oracle")
            if item["region_id"] == "region_1":
                require("-" in representations, phase + "/" + item["id"] + ": region one RC maps to archived graph")
            normalized = []
            for reverse in (False, True):
                ledger = []
                for start in (0, 1, 2):
                    synthetic = b"AT"[:start] + sequence
                    synthetic = (synthetic + b"ACGT" * 40)[:150]
                    assert len(synthetic) == 150 and synthetic[start:start + len(sequence)] == sequence
                    if reverse:
                        synthetic = rc(synthetic)
                    label = phase + "/" + item["id"] + ("/RC/" if reverse else "/F/") + str(start)
                    ledger.append(record(synthetic, b"I" * 150, graph, classes, start + 1,
                                         "synthetic_start_" + str(start), label))
                    FIXTURES.append({"id": label, "sequence": synthetic.decode(), "quality": "I" * 150,
                                     "expected_class": item["id"], "expected_footprint": [-start, 150 - start]})
                summary = next(value for value in PRODUCER.summarize_classes(classes, ledger) if value["id"] == item["id"])
                require(summary["experimental_local_eligibility"], phase + "/" + item["id"] + ": three positioned synthetic reads pass " + str(reverse))
                require(summary["distinct_read_ids"] == summary["distinct_footprints"] == len(summary["canonical_sequence_groups"]) == 3,
                        phase + "/" + item["id"] + ": three distinct coupled IDs/footprints/sequence groups " + str(reverse))
                require(sorted(value["footprint"] for value in summary["witness"]) == [[-2, 148], [-1, 149], [0, 150]],
                        phase + "/" + item["id"] + ": expected class-normalized positions " + str(reverse))
                require(summary["full_product_phase"] == "not_evaluated", phase + "/" + item["id"] + ": no full phase")
                normalized.append([(value["read_id"], value["footprint"]) for value in summary["witness"]])
                shuffled = list(reversed(ledger))
                shuffled_summary = next(value for value in PRODUCER.summarize_classes(classes, shuffled) if value["id"] == item["id"])
                require(shuffled_summary == summary, phase + "/" + item["id"] + ": evidence invariant to ledger iteration order " + str(reverse))
            require(normalized[0] == normalized[1], phase + "/" + item["id"] + ": F/RC eligibility witnesses agree")
            # Neither differing IDs nor quality clipping create independent support.
            synthetic = (sequence + b"ACGT" * 40)[:150]
            ledger = [record(synthetic if index % 2 == 0 else rc(synthetic), b"I" * 150, graph, classes,
                             index + 1, "duplicate_" + str(index), "duplicate") for index in range(4)]
            summary = next(value for value in PRODUCER.summarize_classes(classes, ledger) if value["id"] == item["id"])
            require(not summary["experimental_local_eligibility"] and summary["distinct_footprints"] == 1 and len(summary["canonical_sequence_groups"]) == 1,
                    phase + "/" + item["id"] + ": duplicate and RC records do not create diversity")
            for position in (0, len(sequence) // 2, len(sequence) - 1):
                quality = bytearray(b"I" * 150)
                quality[position] = 44  # Q11, including literal class endpoints.
                for reverse in (False, True):
                    seq, scores = (rc(synthetic), bytes(quality[::-1])) if reverse else (synthetic, bytes(quality))
                    damaged = record(seq, scores, graph, classes, label="quality_boundary")
                    matching = [value for value in damaged["observations"] if value["class_id"] == item["id"]]
                    require(bool(matching) and not any(value["qualifying_local_witness"] for value in matching),
                            phase + "/" + item["id"] + ": Q11 cannot rescue endpoint/internal base " + str((position, reverse)))
            repeated = sequence + b"N" + sequence
            quality = b"I" * (len(sequence) + 1) + b"!" * len(sequence)
            ambiguous = record(repeated, quality, graph, classes, label="ambiguous_quality")
            matching = [value for value in ambiguous["observations"] if value["class_id"] == item["id"]]
            require(len(matching) >= 2 and not any(value["qualifying_local_witness"] for value in matching),
                    phase + "/" + item["id"] + ": quality cannot select one of two unfiltered placements")


def matching_controls():
    ids, footprints = ["a", "b", "c"], [(0, 150), (-1, 149), (-2, 148)]
    all_edges = list(itertools.product(ids, footprints))
    for mask in range(512):
        links = defaultdict(set)
        for bit, (read_id, footprint) in enumerate(all_edges):
            if mask & (1 << bit):
                links[read_id].add(footprint)
        expected = 0
        for size in range(1, 4):
            if any(all(footprint in links[read_id] for read_id, footprint in zip(selected_ids, permutation))
                   for selected_ids in itertools.combinations(ids, size)
                   for permutation in itertools.permutations(footprints, size)):
                expected = size
        actual = PRODUCER.corroboration_witness(links)
        assert len(actual) == expected, (mask, actual, expected)
        assert len({row["read_id"] for row in actual}) == len(actual)
        assert len({tuple(row["footprint"]) for row in actual}) == len(actual)
    require(True, "512 exhaustive three-by-three ID/footprint graphs equal brute-force maximum matching")
    links = {"a": set(footprints), "b": {footprints[0]}, "c": {footprints[0]}}
    require(len(PRODUCER.corroboration_witness(links)) == 2, "Uncoupled three-ID/three-footprint cardinalities fail Hall control")


def synthetic_controls():
    rng = random.Random(20260926)
    dna = lambda length: bytes(rng.choice(b"ACGT") for _ in range(length))
    left, repeat, right = dna(27), dna(29), dna(27)
    sequences = [left + repeat * count + right for count in (1, 2, 3, 4)]
    graph = PRODUCER.prepare_graph(graph_from_sequences(sequences))
    classes = PRODUCER.prepare_manifest(manifest_from_sequences(sequences), graph)
    mapped = [record(sequence, b"I" * len(sequence), graph, classes, index + 1, label="repeat_" + str(index + 1))
              for index, sequence in enumerate(sequences)]
    long_runs = [value["histories"]["+"]["q20"]["runs"] for value in mapped]
    require(all(len(runs) == 1 for runs in long_runs), "All one-to-four repeat-copy fixtures remain contiguous observed walks")
    require(max(Counter(long_runs[-1][0]["edges"]).values()) >= 4, "Repeated occurrence vectors retain four visits")
    a, b = (runs[0]["edges"] for runs in long_runs[-2:])
    require(set(a) == set(b) and set(zip(a, a[1:])) == set(zip(b, b[1:])), "Three/four-copy walks have identical edge and local-pair sets")
    require(a != b and len(a) != len(b), "Ordered histories retain distinction discarded by those aggregate sets")
    def edge_counts(pool):
        return Counter(edge for run in pool for edge in set(run["edges"]))

    def link_counts(pool):
        return Counter(link for run in pool for link in {(event["incoming_edge"], event["outgoing_edge"])
                                                         for event in run["branch_events"]})

    require(edge_counts([long_runs[-2][0]] * 3) == edge_counts([long_runs[-1][0]] * 3),
            "Equal-size three/four-copy pools have identical per-read-deduplicated edge support Counter totals")
    require(link_counts([long_runs[-2][0]] * 3) == link_counts([long_runs[-1][0]] * 3),
            "Equal-size three/four-copy pools have identical per-read-deduplicated branch support Counter totals")
    sequence = sequences[0]
    midpoint = len(sequence) // 2
    for variant, label in ((sequence[:midpoint] + b"N" + sequence[midpoint + 1:], "N"),
                           (sequence[:midpoint] + b"X" + sequence[midpoint + 1:], "X"),
                           (sequence[:midpoint] + b"A" + sequence[midpoint:], "insertion"),
                           (sequence[:midpoint] + sequence[midpoint + 1:], "deletion")):
        result = record(variant, b"I" * len(variant), graph, classes, label=label)
        require(not any(value["qualifying_local_witness"] for value in result["observations"]), label + ": no filled gap or alignment rescue")
        require(len(result["histories"]["+"]["q20"]["runs"]) >= 2, label + ": observed graph histories split")
    # A second identical class constitutes an explicit class tie, not a winner.
    tied = manifest_from_sequences([sequence, sequence])
    tied_classes = PRODUCER.prepare_manifest(tied, graph)
    tied_record = record(sequence, b"I" * len(sequence), graph, tied_classes, label="class_tie")
    require(not any(value["qualifying_local_witness"] for value in tied_record["observations"]), "Identical competing class labels remain ambiguous")
    # A literal palindrome gives competing strand observations.
    half = dna(30)
    palindrome = half + rc(half)
    palindrome_graph = PRODUCER.prepare_graph(graph_from_sequences([palindrome]))
    palindrome_classes = PRODUCER.prepare_manifest(manifest_from_sequences([palindrome]), palindrome_graph)
    palindrome_record = record(palindrome, b"I" * len(palindrome), palindrome_graph, palindrome_classes, label="palindrome")
    require(len(palindrome_record["observations"]) == 2 and not any(value["qualifying_local_witness"] for value in palindrome_record["observations"]),
            "Palindromic strand tie cannot become a unique witness")
    # Missing self-loop is unrepresentable; abundant adjacent edges do not repair it.
    homopolymer = dna(27) + b"A" * 40 + dna(27)
    homo_document = graph_from_sequences([homopolymer])
    homo_document["graph"]["edges"] = [edge for edge in homo_document["graph"]["edges"] if edge["edge_id"] != "A" * 19]
    remaining_sequences = [edge["edge_id"].encode() for edge in homo_document["graph"]["edges"]]
    homo_graph = PRODUCER.prepare_graph(graph_from_sequences(remaining_sequences))
    homo_classes = PRODUCER.prepare_manifest(manifest_from_sequences([homopolymer]), homo_graph)
    require(not any(rep["represented"] for rep in homo_classes[0]["representations"]), "Omitted self-loop class explicitly unrepresented")
    homo_record = record(homopolymer, b"I" * len(homopolymer), homo_graph, homo_classes, label="omitted_selfloop")
    require(not any(value["qualifying_local_witness"] for value in homo_record["observations"]), "Omitted self-loop cannot be rescued by graph-adjacent support")
    # Repeat longer than an available read cannot become a full-window witness.
    long_class = sequences[-1]
    fragments = [record(long_class[:80], b"I" * 80, graph, classes, 1, label="long_left"),
                 record(long_class[-80:], b"I" * 80, graph, classes, 2, label="long_right")]
    summary = PRODUCER.summarize_classes(classes, fragments)[-1]
    require(not summary["experimental_local_eligibility"] and summary["qualifying_records"] == 0,
            "Separate marginal end fragments cannot establish a longer repeat window")
    # Unequal mixture: each record's sequence is preserved, not majority consensus.
    mixture = [copy.deepcopy(mapped[0]) for _ in range(7)] + [copy.deepcopy(mapped[1]) for _ in range(3)]
    for index, value in enumerate(mixture):
        value["ordinal"], value["read_id"] = index + 1, "mixture_" + str(index)
    summary = PRODUCER.summarize_classes(classes, mixture)
    require([value["qualifying_records"] for value in summary[:2]] == [7, 3], "Unequal mixture retains both observed repeat lengths")
    require(all(value["full_product_phase"] == "not_evaluated" for value in summary), "No mixture consensus is promoted to full phase")
    # Distant choices carried by a single read are absent from aggregate support.
    flank_left, middle, flank_right = dna(27), dna(62), dna(27)
    templates = [flank_left + first + middle + second + flank_right
                 for first, second in ((b"A", b"A"), (b"C", b"C"), (b"A", b"C"), (b"C", b"A"))]
    phase_graph = PRODUCER.prepare_graph(graph_from_sequences(templates))
    phase_runs = []
    for sequence in templates:
        histories = PRODUCER.read_histories(sequence, b"I" * len(sequence), phase_graph)
        verify_histories(sequence, b"I" * len(sequence), phase_graph, histories, "distant_phase")
        assert len(histories["+"]["q20"]["runs"]) == 1
        phase_runs.append(histories["+"]["q20"]["runs"][0])
    require(edge_counts(phase_runs[:2]) == edge_counts(phase_runs[2:]),
            "AA+CC versus AC+CA read pools have identical per-read-deduplicated edge support Counter totals")
    require(link_counts(phase_runs[:2]) == link_counts(phase_runs[2:]),
            "AA+CC versus AC+CA read pools have identical per-read-deduplicated branch support Counter totals")
    require({run["sequence"] for run in phase_runs[:2]} != {run["sequence"] for run in phase_runs[2:]},
            "Positioned whole-read histories preserve the phase distinction lost by those Counter totals")
    # Canonicalization includes true IUPAC reverse complements outside the class.
    ambiguous = b"R" + sequences[0] + b"K"
    pair = [record(ambiguous, b"I" * len(ambiguous), graph, classes, 1, label="IUPAC"),
            record(rc(ambiguous), b"I" * len(ambiguous), graph, classes, 2, label="IUPAC_RC")]
    require(PRODUCER.reverse_complement(ambiguous) == rc(ambiguous), "IUPAC complements agree with independent complement table")
    summary = PRODUCER.summarize_classes(classes, pair)[0]
    require(len(summary["canonical_sequence_groups"]) == 1 and summary["distinct_footprints"] == 1,
            "IUPAC/RC full-read copies share one sequence group and one normalized footprint")


def expect_failure(call, label):
    try:
        call()
    except (ValueError, FileExistsError):
        require(True, label)
    else:
        raise AssertionError("Expected failure: " + label)


def scanner_io_controls():
    rng = random.Random(314159)
    literal = lambda size: bytes(rng.choice(b"ACGT") for _ in range(size))
    represented, absent = literal(70), literal(71)
    graph = PRODUCER.prepare_graph(graph_from_sequences([represented]))
    classes = PRODUCER.prepare_manifest(manifest_from_sequences([represented, absent]), graph)
    with tempfile.TemporaryDirectory(prefix="independent-histories-") as temporary:
        directory = Path(temporary)
        fastq = directory / "synthetic.fastq"
        sequences = [represented, absent]
        raw = b"".join(b"@synthetic_" + str(index).encode() + b"\n" + sequence + b"\n+\n" + b"I" * len(sequence) + b"\n"
                       for index, sequence in enumerate(sequences))
        expected = {"records": 2, "bases": sum(map(len, sequences)), "size_bytes": len(raw), "sha256": sha(raw)}
        fastq.write_bytes(raw)
        report = PRODUCER.scan_fastq(fastq, expected, graph, classes)
        require(report["consumed_prefix"] == expected and report["counts"]["retained_records"] == 2,
                "Synthetic scanner preserves byte-prefix identity and both represented/unrepresented literal records")
        require(report["classes"][1]["exact_records"] == 1 and not report["classes"][1]["experimental_local_eligibility"],
                "Literal no-graph class survives prefilter and remains explicitly ineligible")
        require(not any(report["records"][1]["histories"][strand][layer]["runs"]
                        for strand in ("+", "-") for layer in ("unfiltered", "q20")),
                "Unrepresented literal occurrence has no fabricated graph history")
        require("".join(value["raw_fastq"] for value in report["records"]).encode() == raw,
                "Synthetic raw FASTQ ledger reproduces exact original bytes")
        for field, value in (("records", True), ("records", "2"), ("records", 1000001),
                             ("bases", False), ("size_bytes", 0), ("sha256", "0" * 64)):
            expect_failure(lambda: PRODUCER.scan_fastq(fastq, {**expected, field: value}, graph, classes),
                           "Typed prefix/hash failure for " + field + "=" + repr(value))
        for name, value in (("MAX_HISTORY_RECORDS", 1), ("MAX_TOTAL_EDGE_VISITS", 1),
                            ("MAX_READ_LENGTH", 30), ("MAX_LINE_BYTES", 20)):
            with patch.object(PRODUCER, name, value):
                expect_failure(lambda: PRODUCER.scan_fastq(fastq, expected, graph, classes), "Fail closed at " + name)
        for malformed in (raw.replace(b"\n+\n", b"\n-\n", 1), raw[:-2] + b" \n", raw[:-1],
                          raw.replace(represented, b"!" + represented[1:], 1)):
            fastq.write_bytes(malformed)
            expect_failure(lambda: PRODUCER.scan_fastq(fastq, expected, graph, classes), "Malformed FASTQ fails before success publication")
        output = directory / "output.json"
        PRODUCER.publish_json(output, {"original": True})
        original_identity = identity(output)
        expect_failure(lambda: PRODUCER.publish_json(output, {"modified": True}), "Atomic publication refuses existing output")
        require(identity(output) == original_identity, "No-clobber publication preserves original bytes")
        for limit in (1, 2):
            capped = directory / ("capped_" + str(limit) + ".json")
            with patch.object(PRODUCER, "MAX_OUTPUT_BYTES", limit):
                expect_failure(lambda: PRODUCER.publish_json(capped, {}), "Output byte cap includes final newline at " + str(limit))
            require(not capped.exists(), "Output cap cannot expose partial JSON")
        require(not list(directory.glob(".read-histories-*")), "Failed output publication cleans temporary files")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Independent review output already exists")
    before = identity(SCANNER)
    status, failure = "pass", None
    try:
        actual_class_controls(args.inputs)
        matching_controls()
        synthetic_controls()
        scanner_io_controls()
        require(identity(SCANNER) == before, "Prototype source unchanged during independent checks")
    except Exception as error:
        import traceback
        status, failure = "fail", traceback.format_exc()
        print(failure, file=sys.stderr)
    payload = {"status": status, "scope": "Synthetic calibration fixtures only; no biological FASTQ opened",
               "scanner": before, "review_source": identity(Path(__file__)), "python": sys.version,
               "inputs": {path.name: identity(path) for path in sorted(args.inputs.glob("*.json"))},
               "checks_passed": len(CHECKS), "checks": CHECKS, "synthetic_150bp_fixtures": FIXTURES,
               "failure": failure, "limitations": "Fixtures derived from archived classes test implementation sensitivity, not independent biological truth."}
    with args.output.open("x") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": status, "checks_passed": len(CHECKS), "output": str(args.output)}))
    return int(status != "pass")


if __name__ == "__main__":
    sys.exit(main())
