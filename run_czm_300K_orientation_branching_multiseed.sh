#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-/opt/homebrew/Caskroom/miniconda/base/envs/${CONDA_ENV}/bin/python}
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN=$(command -v python)
fi

CLASS=${CLASS:-ceramic}
T_K=${T_K:-300}
THETAS=${THETAS:-"0 15 30 45"}
BRANCH_THETA=${BRANCH_THETA:-30}
SEEDS=${SEEDS:-"1201 1202 1203 1204 1205"}
OUTROOT=${OUTROOT:-runs/czm_Rcurve_300K_orientation_branching_${CLASS}_multiseed_v3}

"$PYTHON_BIN" run_czm_300K_orientation_branching_multiseed.py \
  --project "$(pwd)" \
  --python-bin "$PYTHON_BIN" \
  --class "$CLASS" \
  --temperature "$T_K" \
  --thetas "$THETAS" \
  --branch-theta "$BRANCH_THETA" \
  --seeds "$SEEDS" \
  --outroot "$OUTROOT" \
  "$@"
