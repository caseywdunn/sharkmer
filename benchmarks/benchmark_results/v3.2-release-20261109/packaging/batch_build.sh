#!/bin/bash
#SBATCH --job-name=sharkmer-conda-build
#SBATCH --partition=devel
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tmp/release-3.2.0-20261109/packaging/logs/build-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
module reset
source tmp/release-3.2.0-20261109/packaging/environment.sh
export PATH="$PACKAGING_ROOT/envs/build/bin:$PATH"
python - <<'PY'
import os
from pathlib import Path
import re
recipe = Path(os.environ['PACKAGING_ROOT']) / 'recipes/sharkmer/meta.yaml'
text = recipe.read_text()
assert '{% set version = "3.2.0" %}' in text, 'Expected final 3.2.0 recipe'
assert re.search(r'^  sha256: [0-9a-f]{64}$', text, re.M), 'Final source SHA required'
PY
date --iso-8601=seconds
sha256sum "$PACKAGING_ROOT/recipes/sharkmer/meta.yaml" "$PACKAGING_ROOT/recipes/sharkmer/build.sh"
conda build "$PACKAGING_ROOT/recipes/sharkmer" \
    --override-channels --channel conda-forge --channel bioconda \
    --no-anaconda-upload --croot "$CONDA_BLD_PATH" \
    --output-folder "$PACKAGING_ROOT/packages"
date --iso-8601=seconds
