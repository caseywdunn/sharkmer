import hashlib
import json
import random
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml


WORKSPACE = Path("/tmp/sharkmer-release-comparison")
BUILDS_PATH = WORKSPACE / "builds.json"
OUTPUT_ROOT = WORKSPACE / f"supplemental-synthetic-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sequence_sha256(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def write_fastq(path, sequences):
    with path.open("w") as output:
        for read_index, sequence in enumerate(sequences):
            output.write(f"@read{read_index}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def parse_fasta(path):
    records = []
    header = None
    sequence_lines = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                sequence = "".join(sequence_lines)
                records.append(
                    {
                        "header": header,
                        "length": len(sequence),
                        "sha256": sequence_sha256(sequence),
                        "sequence": sequence,
                    }
                )
            header = line[1:]
            sequence_lines = []
        else:
            sequence_lines.append(line)
    if header is not None:
        sequence = "".join(sequence_lines)
        records.append(
            {
                "header": header,
                "length": len(sequence),
                "sha256": sequence_sha256(sequence),
                "sequence": sequence,
            }
        )
    return records


def make_threshold_input(inputs_directory):
    random_generator = random.Random(13120260909)
    prefix = "".join(random_generator.choice("ACGT") for _ in range(25))
    suffix = "".join(random_generator.choice("ACGT") for _ in range(25))
    abundant = prefix + "".join(random_generator.choice("ACGT") for _ in range(350)) + suffix
    rare = prefix + "".join(random_generator.choice("ACGT") for _ in range(130)) + suffix
    reverse_primer = suffix[-15:].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    input_path = inputs_directory / "threshold.fastq"
    write_fastq(input_path, [abundant] * 100 + [rare] * 4)
    return {
        "case": "threshold",
        "input_path": input_path,
        "sample": "threshold",
        "gene": "target",
        "primer": (
            f"name=target,forward={prefix[:15]},reverse={reverse_primer},trim=15,"
            "mismatches=0,min-length=170,max-length=210,dedup-edit-threshold=0"
        ),
        "expected_sequences": [rare],
        "description": "Known threshold regression: rare valid product after abundant invalid path.",
    }


def make_repeat_inputs(inputs_directory):
    random_generator = random.Random(13220260909)
    prefix = "".join(random_generator.choice("ACGT") for _ in range(59)) + "C"
    suffix = "G" + "".join(random_generator.choice("ACGT") for _ in range(59))
    reverse_primer = suffix[-15:].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    cases = []
    for label, repeat_sequence in (("a18", "A" * 18), ("a19", "A" * 19), ("a40", "A" * 40), ("ac40", "AC" * 20)):
        sequence = prefix + repeat_sequence + suffix
        input_path = inputs_directory / f"repeat-{label}.fastq"
        write_fastq(input_path, [sequence] * 4)
        cases.append(
            {
                "case": f"repeat-{label}",
                "input_path": input_path,
                "sample": "repeat",
                "gene": "target",
                "primer": (
                    f"name=target,forward={prefix[:15]},reverse={reverse_primer},trim=15,"
                    "mismatches=0,min-length=100,max-length=220,dedup-edit-threshold=0"
                ),
                "expected_sequences": [sequence],
                "description": "Known synthetic repeat regression; candidate behavior may conservatively withhold uncertain repeats.",
            }
        )
    return cases


def run_case(version_name, binary, binary_sha256, case, output_directory):
    case_directory = output_directory / version_name / case["case"]
    case_directory.mkdir(parents=True, exist_ok=False)
    command = [
        str(binary),
        "-s",
        case["sample"],
        "-o",
        str(case_directory),
        "-k",
        "19",
        "-t",
        "2",
        "--pcr-primers",
        case["primer"],
        str(case["input_path"]),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=90)
    (case_directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    (case_directory / "stdout.log").write_text(completed.stdout)
    (case_directory / "stderr.log").write_text(completed.stderr)
    stats_path = case_directory / f"{case['sample']}.stats.yaml"
    fasta_path = case_directory / f"{case['sample']}_{case['gene']}.fasta"
    stats = yaml.safe_load(stats_path.read_text()) if stats_path.is_file() else None
    products = parse_fasta(fasta_path) if fasta_path.is_file() else []
    observed_sequences = [product["sequence"] for product in products]
    expected_sequences = case["expected_sequences"]
    return {
        "version": version_name,
        "binary": str(binary),
        "binary_sha256": binary_sha256,
        "case": case["case"],
        "description": case["description"],
        "command": command,
        "returncode": completed.returncode,
        "input": {
            "path": str(case["input_path"]),
            "sha256": sha256(case["input_path"]),
            "size_bytes": case["input_path"].stat().st_size,
        },
        "expected": [
            {"length": len(sequence), "sha256": sequence_sha256(sequence)}
            for sequence in expected_sequences
        ],
        "expected_output_path": str(fasta_path),
        "stats_path": str(stats_path),
        "stats_sha256": sha256(stats_path) if stats_path.is_file() else None,
        "stats_pcr_results": stats.get("pcr_results") if isinstance(stats, dict) else None,
        "products": products,
        "all_products_exact": observed_sequences == expected_sequences,
        "has_incorrect_product": any(sequence not in expected_sequences for sequence in observed_sequences),
        "raw_output_directory": str(case_directory),
    }


def main():
    builds = json.loads(BUILDS_PATH.read_text())
    if OUTPUT_ROOT.exists():
        raise RuntimeError(f"Supplemental output already exists: {OUTPUT_ROOT}")
    inputs_directory = OUTPUT_ROOT / "inputs"
    inputs_directory.mkdir(parents=True)
    cases = [make_threshold_input(inputs_directory), *make_repeat_inputs(inputs_directory)]
    results = []
    for version_name, build in builds["versions"].items():
        binary = Path(build["binary_path"])
        if sha256(binary) != build["binary_sha256"]:
            raise RuntimeError(f"Binary checksum changed: {version_name}")
        for case in cases:
            results.append(run_case(version_name, binary, build["binary_sha256"], case, OUTPUT_ROOT))
    summary = {
        "schema_version": 1,
        "label": "supplemental known synthetic regressions; not biological calibration or timed benchmarking",
        "build_manifest": str(BUILDS_PATH),
        "output_root": str(OUTPUT_ROOT),
        "results": results,
    }
    (OUTPUT_ROOT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    compact = [
        {
            "version": result["version"],
            "case": result["case"],
            "returncode": result["returncode"],
            "all_products_exact": result["all_products_exact"],
            "has_incorrect_product": result["has_incorrect_product"],
            "product_lengths": [product["length"] for product in result["products"]],
            "product_hashes": [product["sha256"] for product in result["products"]],
            "input_sha256": result["input"]["sha256"],
            "binary_sha256": result["binary_sha256"],
        }
        for result in results
    ]
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
