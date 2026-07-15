#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
OUTROOT=${OUTROOT:-runs/mpz_v9_14_event_driven_remesh_700K_v1}
SEEDS=${SEEDS:-"1"}
CLASSES=${CLASSES:-"ceramic weakT DBTT"}
T_K=${T_K:-700}
TARGET_EXT_UM=${TARGET_EXT_UM:-50}
STEPS=${STEPS:-12000}
NX=${NX:-36}
NY=${NY:-72}
TIP_H_FINE=${TIP_H_FINE:-1e-6}
TIP_RATIO=${TIP_RATIO:-1.20}
DU=${DU:-2e-7}
DT=${DT:-8.4}
MPZ_LENGTH_UM=${MPZ_LENGTH_UM:-100}
MPZ_N_BINS=${MPZ_N_BINS:-200}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-5}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-10}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

"$PYTHON_BIN" -m pytest -q tests/test_event_driven_remesh_v914.py
mkdir -p "$OUTROOT"

EVENT_STATISTICS=deterministic \
STOCHASTIC_EMISSION=0 \
PROPAGATION_CONTROL=raw \
RNG_COUPLING=common \
"$PYTHON_BIN" run_mpz_v9_14_event_driven_remesh.py \
  --outroot "$OUTROOT" \
  --seeds "$SEEDS" \
  --classes "$CLASSES" \
  --T-K "$T_K" \
  --target-extension-um "$TARGET_EXT_UM" \
  --steps "$STEPS" \
  --nx "$NX" --ny "$NY" \
  --tip-h-fine "$TIP_H_FINE" --tip-ratio "$TIP_RATIO" \
  --dU "$DU" --dt "$DT" \
  --mpz-length-um "$MPZ_LENGTH_UM" --mpz-n-bins "$MPZ_N_BINS" \
  --save-snapshots "$SAVE_SNAPSHOTS" --snapshot-cols "$SNAPSHOT_COLS" \
  --snapshot-by-extension-um "$SNAPSHOT_BY_EXT_UM" \
  --event-statistics deterministic --no-stochastic-emission \
  --propagation-control raw --rng-coupling common \
  2>&1 | tee "$OUTROOT/driver.log"
