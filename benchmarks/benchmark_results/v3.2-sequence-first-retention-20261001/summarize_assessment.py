#!/usr/bin/env python3
"""Produce compact tables from the completed paired comparison."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path


def table(path, fields, rows):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.summary.read_text())
    if source["status"] != "pass" or len(source["cells"]) != 10:
        raise ValueError("Need completed ten-sample comparison")
    paired = {row["accession"]: row for row in source["paired"]["cells"]}
    if source["totals"]["standard"] != 64 or source["totals"]["old_uncertain"] != 62:
        raise ValueError("Old inclusive output did not reproduce 64 standard/62 uncertain records")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    counts, resources, gains = [], [], []
    for cell in source["cells"]:
        accession = cell["accession"]
        exact = paired[accession]
        by_gene = Counter(row["gene"] for row in cell["uncertain_gains"])
        counts.append({"accession": accession, "standard": cell["standard"],
                       "old_uncertain": cell["old_uncertain"], "new_uncertain": cell["new_uncertain"],
                       "new_uncertain_identities": len(cell["uncertain_gains"]),
                       "gained_targets": ",".join(f"{gene}:{count}" for gene, count in sorted(by_gene.items()))})
        wall = exact["resource_changes"]["wall_seconds"]
        rss = exact["resource_changes"]["peak_rss_bytes"]
        resources.append({"accession": accession,
                          "old_wall_s": exact["roles"]["baseline"]["metrics"]["wall_seconds"]["median"],
                          "new_wall_s": exact["roles"]["candidate"]["metrics"]["wall_seconds"]["median"],
                          "wall_percent": wall["percent"], "wall_priority": wall["priority_review"],
                          "old_rss_mib": exact["roles"]["baseline"]["metrics"]["peak_rss_bytes"]["median"] / 2**20,
                          "new_rss_mib": exact["roles"]["candidate"]["metrics"]["peak_rss_bytes"]["median"] / 2**20,
                          "rss_percent": rss["percent"], "rss_priority": rss["priority_review"]})
        gains.extend({"accession": accession, **row} for row in cell["uncertain_gains"])
    table(output / "sample_counts.tsv",
          ["accession", "standard", "old_uncertain", "new_uncertain", "new_uncertain_identities", "gained_targets"],
          counts)
    table(output / "resources.tsv",
          ["accession", "old_wall_s", "new_wall_s", "wall_percent", "wall_priority",
           "old_rss_mib", "new_rss_mib", "rss_percent", "rss_priority"], resources)
    table(output / "uncertain_gains.tsv",
          ["accession", "gene", "length", "sha256"], gains)
    total = source["totals"]["standard"] + source["totals"]["new_uncertain"]
    result = {"status": "pass", "old_total": 126, "new_total": total,
              "standard": source["totals"]["standard"],
              "old_uncertain": source["totals"]["old_uncertain"],
              "new_uncertain": source["totals"]["new_uncertain"],
              "new_uncertain_fraction": source["totals"]["new_uncertain"] / total,
              "gained_uncertain_identities": source["totals"]["gained_uncertain"],
              "sum_median_wall_old_s": source["paired"]["median_wall_sum_seconds"]["baseline"],
              "sum_median_wall_new_s": source["paired"]["median_wall_sum_seconds"]["candidate"],
              "wall_priority_accessions": [row["accession"] for row in resources if row["wall_priority"]],
              "rss_priority_accessions": [row["accession"] for row in resources if row["rss_priority"]]}
    if sum(row["new_uncertain_identities"] for row in counts) != result["gained_uncertain_identities"]:
        raise ValueError("Gained identity totals disagree")
    (output / "report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
