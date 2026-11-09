#!/bin/bash
#SBATCH --job-name=sharkmer-conda-setup
#SBATCH --partition=devel
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=tmp/release-3.2.0-20261109/packaging/logs/setup-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
module reset
source tmp/release-3.2.0-20261109/packaging/environment.sh
date --iso-8601=seconds
/apps/software/system/software/miniconda/26.9.4/bin/conda create \
    --yes --prefix "$PACKAGING_ROOT/envs/build" --override-channels \
    --channel conda-forge --no-default-packages python=3.12 conda-build
"$PACKAGING_ROOT/envs/build/bin/conda" list --prefix "$PACKAGING_ROOT/envs/build" \
    --explicit > "$PACKAGING_ROOT/build-environment.explicit.txt"
"$PACKAGING_ROOT/envs/build/bin/conda-build" --version
date --iso-8601=seconds
