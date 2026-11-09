#!/bin/bash
# Reuse installed tooling and cached packages; keep new writable state separate.
PACKAGING_ROOT=/nfs/roberts/project/pi_cwd7/cwd7/sharkmer/tmp/release-3.2.1-20261109/packaging
PACKAGING_PREVIOUS=/nfs/roberts/project/pi_cwd7/cwd7/sharkmer/tmp/release-3.2.0-20261109/packaging
export PACKAGING_ROOT PACKAGING_PREVIOUS
export CONDARC="$PACKAGING_ROOT/condarc"
export CONDA_PKGS_DIRS="$PACKAGING_ROOT/conda-pkgs:$PACKAGING_PREVIOUS/conda-pkgs"
export CONDA_ENVS_PATH="$PACKAGING_ROOT/envs"
export CONDA_REGISTER_ENVS=false
export CONDA_NUMBER_CHANNEL_NOTICES=0
export XDG_CACHE_HOME="$PACKAGING_ROOT/cache"
export TMPDIR="$PACKAGING_ROOT/tmp"
export PYTHONNOUSERSITE=1
export CONDA_BLD_PATH="$PACKAGING_ROOT/conda-bld"
export CARGO_HOME="$PACKAGING_ROOT/cargo"
export CARGO_BUILD_JOBS=2 CPU_COUNT=2
mkdir -p "$PACKAGING_ROOT/conda-pkgs" "$CONDA_ENVS_PATH" "$XDG_CACHE_HOME" "$TMPDIR" "$CARGO_HOME"
