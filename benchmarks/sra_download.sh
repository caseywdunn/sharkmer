#!/usr/bin/env zsh
set -euo pipefail

# Download the first N reads from each ENA FASTQ file for a run accession.
#
# Usage:
#   ./sra_download.sh ERR571460 1000
#
# Output for paired-end runs:
#   ERR571460_1_1000.fastq
#   ERR571460_2_1000.fastq
#
# Output for single-end runs:
#   ERR571460_1_1000.fastq
#
# Notes:
# - This streams ENA-hosted FASTQ files and stops as soon as N reads have
#   been written, so it does not intentionally download the full run.
# - If N is larger than the number of reads available, all reads are written.
# - FASTQ records are 4 lines per read.

usage() {
  cat >&2 <<'EOF'
Usage:
  sra_download.sh ACCESSION NREADS

Arguments:
  ACCESSION   ENA/SRA run accession, e.g. ERR571460
  NREADS      Number of reads to keep from the start of each FASTQ file

Examples:
  sra_download.sh ERR571460 1000
  sra_download.sh SRR12345678 50000
EOF
  exit 1
}

ena_first_reads() {
  local acc="$1"
  local n="$2"
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

ena_first_reads "${1:-}" "${2:-}"
