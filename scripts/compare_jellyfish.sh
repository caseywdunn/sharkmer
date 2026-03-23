#!/usr/bin/env bash
#
# Compare sharkmer and jellyfish kmer histograms to verify they produce
# identical results. Requires both sharkmer and jellyfish on PATH.
#
# Usage: ./scripts/compare_jellyfish.sh <fastq_file> [k] [histo_max]

set -euo pipefail

FASTQ="${1:?Usage: $0 <fastq_file> [k] [histo_max]}"
K="${2:-21}"
HISTO_MAX="${3:-10000}"

if [[ ! -f "$FASTQ" ]]; then
    echo "Error: file not found: $FASTQ" >&2
    exit 1
fi

for cmd in sharkmer jellyfish; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "Error: $cmd not found on PATH" >&2
        exit 1
    fi
done

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "Comparing sharkmer and jellyfish histograms"
echo "  FASTQ:    $FASTQ"
echo "  k:        $K"
echo "  histo_max: $HISTO_MAX"
echo "  tmpdir:   $TMPDIR"
echo

# --- sharkmer ---
echo "Running sharkmer..."
sharkmer -s compare -k "$K" --chunks 1 --histo-max "$HISTO_MAX" -o "$TMPDIR/sharkmer" "$FASTQ" 2>/dev/null

# Extract data lines (skip comment, header, and zero-frequency rows)
grep -v '^#' "$TMPDIR/sharkmer/compare.final.histo" \
    | grep -v '^count' \
    | awk -F'\t' '$2 > 0 {print $1, $2}' \
    > "$TMPDIR/sharkmer.hist"

# --- jellyfish ---
echo "Running jellyfish..."

# Decompress if gzipped (jellyfish cannot read from process substitution)
if [[ "$FASTQ" == *.gz ]]; then
    gunzip -c "$FASTQ" > "$TMPDIR/input.fastq"
    JELLY_INPUT="$TMPDIR/input.fastq"
else
    JELLY_INPUT="$FASTQ"
fi

jellyfish count -m "$K" -s 100M -t 1 -C "$JELLY_INPUT" -o "$TMPDIR/jf.jf"
jellyfish histo -h "$HISTO_MAX" "$TMPDIR/jf.jf" > "$TMPDIR/jellyfish.hist"

# --- Compare ---
echo
echo "Comparing histograms..."

# Both should now be space-separated: count frequency
# jellyfish output is already in that format
# Diff the two
if diff -q "$TMPDIR/sharkmer.hist" "$TMPDIR/jellyfish.hist" > /dev/null 2>&1; then
    echo "PASS: Histograms are identical."
    exit 0
else
    echo "FAIL: Histograms differ."
    echo
    echo "First differences (< sharkmer, > jellyfish):"
    diff "$TMPDIR/sharkmer.hist" "$TMPDIR/jellyfish.hist" | head -20
    exit 1
fi
