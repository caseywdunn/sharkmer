#!/usr/bin/env zsh
set -euo pipefail

# Download the first N reads from ENA FASTQ files for a run accession.
#
# Usage:
#   ./sra_download.sh ERR571460 1000          # R1 only (default)
#   ./sra_download.sh ERR571460 1000 both     # R1 and R2
#
# Output (R1 only, default):
#   ERR571460_1_1000.fastq
#
# Output (both mates):
#   ERR571460_1_1000.fastq
#   ERR571460_2_1000.fastq
#
# Notes:
# - This streams ENA-hosted FASTQ files and stops as soon as N reads have
#   been written, so it does not intentionally download the full run.
# - If N is larger than the number of reads available, all reads are written.
# - FASTQ records are 4 lines per read.
# - Benchmarks only use R1, so defaulting to R1 saves bandwidth and time.

usage() {
  cat >&2 <<'EOF'
Usage:
  sra_download.sh ACCESSION NREADS [both]

Arguments:
  ACCESSION   ENA/SRA run accession, e.g. ERR571460
  NREADS      Number of reads to keep from the start of each FASTQ file
  both        Optional: download both R1 and R2 (default: R1 only)

Examples:
  sra_download.sh ERR571460 1000          # R1 only
  sra_download.sh ERR571460 1000 both     # both mates
EOF
  exit 1
}

ena_first_reads() {
  local acc="$1"
  local n="$2"
  local mode="${3:-r1}"
  local report
  local fastq_field

  [[ -n "$acc" && -n "$n" ]] || usage
  [[ "$n" =~ '^[0-9]+$' ]] || {
    echo "Error: NREADS must be a positive integer." >&2
    exit 1
  }
  (( n > 0 )) || {
    echo "Error: NREADS must be greater than zero." >&2
    exit 1
  }

  report="$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=run_accession,fastq_ftp")" || {
    echo "Error: could not retrieve metadata from ENA for accession '${acc}'." >&2
    exit 1
  }

  fastq_field="$(printf '%s\n' "$report" | awk -F'\t' 'NR==2 {print $2}')"

  if [[ -z "${fastq_field:-}" ]]; then
    echo "Error: accession '${acc}' was not found, or ENA does not report FASTQ files for it." >&2
    echo "Try:" >&2
    echo "  curl -s 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=run_accession,fastq_ftp,submitted_ftp'" >&2
    exit 1
  fi

  printf '%s\n' "$fastq_field" \
    | tr ';' '\n' \
    | awk 'NF' \
    | nl -v 1 -w 1 -s $'\t' \
    | while IFS=$'\t' read -r mate url; do
        # Skip R2 unless 'both' was requested
        if [[ "$mode" != "both" && "$mate" -gt 1 ]]; then
          echo "Skipping mate ${mate} (pass 'both' to download all mates)" >&2
          continue
        fi

        local outfile="${acc}_${mate}_${n}.fastq"

        echo "Downloading mate ${mate} from https://${url}" >&2
        echo "Writing first ${n} reads to ${outfile}" >&2

        # Expected behavior: once awk has written N reads, it exits, which
        # closes the pipe upstream. curl/gzip may then exit nonzero because
        # the downstream pipe closed early.
        {
          curl -fsSL "https://${url}" \
            | gzip -dc \
            | awk -v n="$n" 'NR <= 4*n { print } NR == 4*n { exit }'
        } > "$outfile" || true

        local lines reads
        lines=$(wc -l < "$outfile")
        reads=$(( lines / 4 ))

        if (( reads < n )); then
          echo "Note: ${outfile} contains ${reads} reads because the source file ended before ${n} reads were available." >&2
        else
          echo "Done: ${outfile} contains ${reads} reads." >&2
        fi
      done
}

ena_first_reads "${1:-}" "${2:-}" "${3:-r1}"
