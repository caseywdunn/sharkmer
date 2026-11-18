#!/bin/bash
#SBATCH --job-name=v4-direct-paired
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=tmp/v4-direct-ingestion-20261110/logs/paired-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
module reset
export PATH="/nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin:$PATH"
export LC_ALL=C
date --iso-8601=seconds
hostname
lscpu
python --version
sha256sum benchmarks/compare_count_backends.py tmp/v4-direct-ingestion-20261110/candidate-v1-source.json tmp/v4-direct-ingestion-20261110/bin/sharkmer-candidate-v1 tmp/v4-performance-baseline-20261109/bin/sharkmer-v3.2.1
python - <<'PY'
from pathlib import Path
import hashlib,json
for path, expected in json.loads(Path('tmp/v4-direct-ingestion-20261110/candidate-v1-source.json').read_text()).items():
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected, path
PY
python benchmarks/compare_count_backends.py \
  --protocol benchmarks/benchmark_results/v4-performance-baseline-20261109/protocol.json \
  --baseline tmp/v4-performance-baseline-20261109/bin/sharkmer-v3.2.1 \
  --candidate tmp/v4-direct-ingestion-20261110/bin/sharkmer-candidate-v1 \
  --output tmp/v4-direct-ingestion-20261110/paired-v1
date --iso-8601=seconds
sha256sum benchmarks/compare_count_backends.py tmp/v4-direct-ingestion-20261110/bin/sharkmer-candidate-v1
