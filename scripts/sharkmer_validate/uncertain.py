"""Validate separately published uncertain hypotheses without promoting them."""

import hashlib
import re
import shlex
from collections import Counter


CAUSE_COUNTS = {
    "pre_prune_omitted_self_loop_node": "omitted_self_loop_marker_occurrences",
    "pre_prune_cyclic_scc_node": "cyclic_scc_marker_occurrences",
    "retained_collision_edge": "retained_collision_marker_occurrences",
}
CAUSE_IDENTITIES = {
    "pre_prune_omitted_self_loop_node": "oriented_node_sequence",
    "pre_prune_cyclic_scc_node": "oriented_node_sequence",
    "retained_collision_edge": "oriented_edge_sequence",
}


def _require(condition, message):
    if not condition:
        raise ValueError(f"Invalid uncertain products: {message}")


def _count(value):
    return type(value) is int and value >= 0


def _digest(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def uncertain_header(header):
    # Inspect every token: a duplicate later field cannot erase this label.
    return "assembly_status=uncertain" in shlex.split(header)


def validate_metadata(entry, sample_prefix):
    """Return the declared uncertain filename, or None for absent/empty output."""
    if "uncertain_products" not in entry:
        return None
    value = entry["uncertain_products"]
    _require(isinstance(value, dict), "assessment must be a mapping")
    _require(type(value.get("schema_version")) is int and value["schema_version"] == 1,
             "unsupported schema")
    for key, expected in (
        ("purpose", "uncertain_repeat_structure_graph_path_hypotheses"),
        ("coordinate_system", "zero_based_half_open_in_emitted_oriented_sequence"),
        ("read_support_evaluation", "not_evaluated"),
        ("repeat_length_evaluation", "unresolved"),
        ("search_limit_scope", "all_attempted_thresholds_for_gene"),
    ):
        _require(value.get(key) == expected, f"invalid {key}")
    for key in ("path_search_evaluated", "search_limits_reached", "retention_truncated"):
        _require(type(value.get(key)) is bool, f"invalid {key}")
    for key in ("n_products", "unique_products", "retained_path_records",
                "excluded_standard_equivalent_path_records", "exact_duplicate_path_records"):
        _require(_count(value.get(key)), f"invalid {key}")
    products = value.get("products")
    lengths = value.get("product_lengths")
    _require(isinstance(products, list) and isinstance(lengths, list), "missing product lists")
    count = value["n_products"]
    _require(count == value["unique_products"] == len(products) == len(lengths), "count mismatch")
    _require(all(_count(length) and length > 0 for length in lengths), "invalid lengths")
    _require(value["retained_path_records"] == count + value["exact_duplicate_path_records"]
             + value["excluded_standard_equivalent_path_records"], "retention accounting mismatch")
    gene = entry.get("gene_name")
    _require(isinstance(gene, str) and gene, "missing gene")
    filename = value.get("output_file")
    _require(filename == (f"{sample_prefix}_{gene}.uncertain.fasta" if count else None),
             "output path does not match sample/gene/count")

    ids, sequences, source_count = set(), set(), 0
    for index, (product, length) in enumerate(zip(products, lengths)):
        _require(isinstance(product, dict), "product must be a mapping")
        identifier = product.get("id")
        _require(identifier == f"{sample_prefix}_{gene}.uncertain_{index}" and identifier not in ids,
                 "invalid or duplicate product ID")
        ids.add(identifier)
        sequence = product.get("oriented_sequence")
        _require(isinstance(sequence, str) and re.fullmatch("[ACGT]+", sequence) is not None,
                 "invalid literal sequence")
        _require(sequence not in sequences, "duplicate oriented sequence")
        sequences.add(sequence)
        _require(_count(product.get("sequence_length"))
                 and len(sequence) == product["sequence_length"] == length, "sequence length mismatch")
        _require(product.get("sequence_sha256") == _digest(sequence), "sequence hash mismatch")
        sources = product.get("sources")
        _require(isinstance(sources, list) and sources, "missing source paths")
        source_count += len(sources)
        causes, thresholds = set(), set()
        for source in sources:
            _require(isinstance(source, dict) and _count(source.get("threshold"))
                     and source["threshold"] > 0, "invalid source threshold")
            node_limit = source.get("node_visit_limit")
            node_visits = source.get("max_observed_node_visits")
            _require(_count(node_limit) and _count(node_visits)
                     and 1 <= node_visits <= node_limit, "invalid source node visit bounds")
            thresholds.add(source["threshold"])
            markers = source.get("marker_occurrences")
            _require(isinstance(markers, list) and markers, "missing uncertainty markers")
            observed = Counter()
            spans = []
            for marker in markers:
                _require(isinstance(marker, dict), "marker must be a mapping")
                cause = marker.get("cause")
                _require(isinstance(cause, str) and cause in CAUSE_COUNTS, "unknown uncertainty cause")
                _require(marker.get("identity_kind") == CAUSE_IDENTITIES[cause],
                         "marker identity kind disagrees with cause")
                start, end = marker.get("start"), marker.get("end")
                _require(_count(start) and _count(end) and start < end <= length,
                         "marker outside emitted sequence")
                word = sequence[start:end]
                _require(marker.get("identity_sequence") == word
                         and marker.get("identity_sha256") == _digest(word), "marker identity mismatch")
                observed[cause] += 1
                causes.add(cause)
                spans.append((start, end))
            for cause, field in CAUSE_COUNTS.items():
                _require(_count(source.get(field)) and source[field] == observed[cause],
                         "marker cause count mismatch")
            expected_runs = []
            for start, end in sorted(spans):
                if expected_runs and start <= expected_runs[-1]["end"]:
                    expected_runs[-1]["end"] = max(end, expected_runs[-1]["end"])
                else:
                    expected_runs.append({"start": start, "end": end})
            runs = source.get("marker_covered_runs")
            _require(isinstance(runs, list), "missing marker-covered runs")
            for run in runs:
                _require(isinstance(run, dict) and _count(run.get("start"))
                         and _count(run.get("end")) and 0 <= run["start"] < run["end"] <= length,
                         "invalid marker-covered run")
            _require(runs == expected_runs, "marker-covered runs disagree with marker union")
        description = product.get("description")
        _require(isinstance(description, str), "missing header description")
        fields = {}
        for token in description.split():
            _require("=" in token, "malformed header field")
            key, field_value = token.split("=", 1)
            _require(key and key not in fields, "duplicate header field")
            fields[key] = field_value
        expected_fields = {
            "assembly_status": "uncertain", "uncertainty": "repeat_structure",
            "repeat_length": "unresolved", "read_support": "not_evaluated",
            "hypothesis_length": str(length), "sequence_sha256": _digest(sequence),
            "causes": ",".join(sorted(causes)),
            "thresholds": ",".join(str(threshold) for threshold in sorted(thresholds)),
            "source_path_records": str(len(sources)), "enumeration": "bounded",
            "limit_scope": "gene", "search_limits_reached": str(value["search_limits_reached"]).lower(),
            "retention_truncated": str(value["retention_truncated"]).lower(),
        }
        _require(all(fields.get(key) == expected for key, expected in expected_fields.items()),
                 "header disagrees with assessment")
    _require(source_count == value["retained_path_records"]
             - value["excluded_standard_equivalent_path_records"], "source path count mismatch")
    return filename


def validate_fasta(entry, output_dir):
    """Bind every uncertain FASTA record to its separate structured metadata."""
    value = entry.get("uncertain_products")
    if not value or not value.get("output_file"):
        return
    expected = [(product["id"] + " " + product["description"], product["oriented_sequence"])
                for product in value["products"]]
    observed, header, sequence = [], None, []
    with (output_dir / value["output_file"]).open() as stream:
        for raw in stream:
            line = raw.strip()
            if line.startswith(">"):
                if header is not None:
                    observed.append((header, "".join(sequence)))
                header, sequence = line[1:], []
            elif line:
                _require(header is not None, "sequence before FASTA header")
                sequence.append(line)
    if header is not None:
        observed.append((header, "".join(sequence)))
    _require(observed == expected, "FASTA records disagree with structured candidates")
