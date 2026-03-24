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

# --- Timing helper ---
# On macOS, /usr/bin/time -l reports peak RSS in bytes.
# On Linux, /usr/bin/time -v reports peak RSS in KB.
if [[ "$(uname)" == "Darwin" ]]; then
    TIME_CMD=(/usr/bin/time -l)
    parse_peak_rss() {
        # macOS: "maximum resident set size" line reports bytes
        local bytes
        bytes=$(grep 'maximum resident set size' "$1" | awk '{print $1}')
        if [[ -n "$bytes" ]]; then
            echo "$(( bytes / 1024 / 1024 )) MB"
        else
            echo "N/A"
        fi
    }
    parse_wall_time() {
        # macOS /usr/bin/time outputs "N.NN real" on the first line of stderr
        grep 'real' "$1" | awk '{print $1 " s"}'
    }
else
    TIME_CMD=(/usr/bin/time -v)
    parse_peak_rss() {
        # Linux: "Maximum resident set size (kbytes)" line
        local kb
        kb=$(grep 'Maximum resident set size' "$1" | awk '{print $NF}')
        if [[ -n "$kb" ]]; then
            echo "$(( kb / 1024 )) MB"
        else
            echo "N/A"
        fi
    }
    parse_wall_time() {
        grep 'Elapsed (wall clock)' "$1" | sed 's/.*: //'
    }
fi

# --- sharkmer ---
echo "Running sharkmer..."
"${TIME_CMD[@]}" sharkmer -s compare -k "$K" --chunks 1 --histo-max "$HISTO_MAX" -o "$TMPDIR/sharkmer" "$FASTQ" 2>"$TMPDIR/sharkmer_time.txt"
SHARKMER_TIME=$(parse_wall_time "$TMPDIR/sharkmer_time.txt")
SHARKMER_RAM=$(parse_peak_rss "$TMPDIR/sharkmer_time.txt")

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

"${TIME_CMD[@]}" jellyfish count -m "$K" -s 100M -t 1 -C "$JELLY_INPUT" -o "$TMPDIR/jf.jf" 2>"$TMPDIR/jf_count_time.txt"
JF_COUNT_TIME=$(parse_wall_time "$TMPDIR/jf_count_time.txt")
JF_COUNT_RAM=$(parse_peak_rss "$TMPDIR/jf_count_time.txt")

"${TIME_CMD[@]}" jellyfish histo -h "$HISTO_MAX" "$TMPDIR/jf.jf" 2>"$TMPDIR/jf_histo_time.txt" > "$TMPDIR/jellyfish.hist"
JF_HISTO_TIME=$(parse_wall_time "$TMPDIR/jf_histo_time.txt")
JF_HISTO_RAM=$(parse_peak_rss "$TMPDIR/jf_histo_time.txt")

# --- Compare ---
echo
echo "=== Performance ==="
echo "  sharkmer:         time=$SHARKMER_TIME  peak_ram=$SHARKMER_RAM"
echo "  jellyfish count:  time=$JF_COUNT_TIME  peak_ram=$JF_COUNT_RAM"
echo "  jellyfish histo:  time=$JF_HISTO_TIME  peak_ram=$JF_HISTO_RAM"
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
