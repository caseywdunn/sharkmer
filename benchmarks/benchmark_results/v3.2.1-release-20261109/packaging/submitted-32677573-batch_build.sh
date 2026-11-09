#!/bin/bash
#SBATCH --job-name=sharkmer-conda-3.2.1
#SBATCH --partition=devel
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tmp/release-3.2.1-20261109/packaging/logs/build-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
module reset
source tmp/release-3.2.1-20261109/packaging/environment.sh
export PATH="$PACKAGING_PREVIOUS/envs/build/bin:$PATH"
python - <<'PY'
import os
from pathlib import Path
import re
recipe = Path(os.environ['PACKAGING_ROOT']) / 'recipes/sharkmer/meta.yaml'
text = recipe.read_text()
assert '{% set version = "3.2.1" %}' in text, 'Expected final 3.2.1 recipe'
assert re.search(r'^  sha256: [0-9a-f]{64}$', text, re.M), 'Final source SHA required'
PY
date --iso-8601=seconds
conda-build --version
sha256sum "$PACKAGING_ROOT/recipes/sharkmer/meta.yaml" "$PACKAGING_ROOT/recipes/sharkmer/build.sh"
conda build "$PACKAGING_ROOT/recipes/sharkmer" \
    --override-channels --channel conda-forge --channel bioconda \
    --no-anaconda-upload --croot "$CONDA_BLD_PATH" \
    --output-folder "$PACKAGING_ROOT/packages" 2>&1 | tee "$PACKAGING_ROOT/build-output.log"
python - <<'PY'
import os
from pathlib import Path
import re
root = Path(os.environ['PACKAGING_ROOT'])
log = (root / 'build-output.log').read_text()
test = re.search(r'TEST START:(.*?)TEST END:', log, re.S)
assert test is not None, 'Completed installed-package smoke-test block missing'
assert '+ sharkmer --version' in test.group(1), 'Package version command missing'
assert re.search(r'^sharkmer 3\.2\.1$', test.group(1), re.M), 'Packaged version is not3.2.1'
(root / 'package-version.txt').write_text('sharkmer 3.2.1\n')
print('Installed package version verified: sharkmer 3.2.1')
PY
date --iso-8601=seconds
