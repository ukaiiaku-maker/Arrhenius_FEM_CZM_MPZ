#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
MATERIAL=${MATERIAL:-weakT}
T_K=${T_K:-700}
OUTROOT=${OUTROOT:-runs/v10_0_2_progressive_${MATERIAL}_${T_K}K_5um_v1}
STEPS=${STEPS:-50000}
NX=${NX:-36}
NY=${NY:-72}
TIP_H_FINE=${TIP_H_FINE:-1e-6}
TIP_RATIO=${TIP_RATIO:-1.20}
DU=${DU:-2e-7}
DT=${DT:-8.4}
MPZ_LENGTH_UM=${MPZ_LENGTH_UM:-100}
MPZ_N_BINS=${MPZ_N_BINS:-200}
MAX_TRIAL_DAMAGE_CHANGE=${MAX_TRIAL_DAMAGE_CHANGE:-0.02}
MIN_TRIAL_RETRY_DT_S=${MIN_TRIAL_RETRY_DT_S:-1e-18}
MAX_TRIAL_RETRIES=${MAX_TRIAL_RETRIES:-64}
MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL=${MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL:-10000}
MIN_TRIANGLE_QUALITY=${MIN_TRIANGLE_QUALITY:-0.035}
MIN_CHILD_AREA_RATIO=${MIN_CHILD_AREA_RATIO:-0.08}
MAX_TIP_H_OVER_DA=${MAX_TIP_H_OVER_DA:-0.75}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

if [[ -e "$OUTROOT" ]]; then
  echo "ERROR: output path already exists: $OUTROOT"
  echo "Choose a new versioned OUTROOT; this runner does not overwrite prior results."
  exit 1
fi

"$PYTHON_BIN" -m pip install -e . --no-deps
"$PYTHON_BIN" -m compileall -q arrhenius_fracture
"$PYTHON_BIN" -m py_compile \
  audit_v10_0_no_artificial_controls.py \
  audit_v10_0_2_progressive_smoke.py
"$PYTHON_BIN" -m pytest -q \
  tests/test_pf_equivalent_manifest_v10.py \
  tests/test_kinetic_campaign_czm_v10.py \
  tests/test_cohesive_trial_state_v10.py \
  tests/test_kinetic_cohesive_stepper_v10.py \
  tests/test_kinetic_campaign_reset_v1001.py \
  tests/test_kinetic_event_lifecycle_v1002.py \
  tests/test_progressive_run_2d_transform_v10.py \
  tests/test_progressive_event_lifecycle_transform_v1002.py

ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM=5 \
ARRHENIUS_PREFINED_MODE_I_CORRIDOR=1 \
ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY="$MIN_TRIANGLE_QUALITY" \
ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO="$MIN_CHILD_AREA_RATIO" \
ARRHENIUS_MAX_TIP_H_OVER_DA="$MAX_TIP_H_OVER_DA" \
ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE="$MAX_TRIAL_DAMAGE_CHANGE" \
ARRHENIUS_MIN_TRIAL_RETRY_DT_S="$MIN_TRIAL_RETRY_DT_S" \
ARRHENIUS_MAX_TRIAL_RETRIES="$MAX_TRIAL_RETRIES" \
ARRHENIUS_MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL="$MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL" \
"$PYTHON_BIN" -m arrhenius_fracture.mode_i_first_passage_v10_0_2_progressive \
  --v10-material-class "$MATERIAL" \
  --czm-opening-coupling clock_linear \
  --mode 2d \
  --temperatures "$T_K" \
  --steps "$STEPS" \
  --nx "$NX" --ny "$NY" \
  --tip-h-fine "$TIP_H_FINE" --tip-ratio "$TIP_RATIO" \
  --dU "$DU" --dt "$DT" \
  --da-phys-um 5 \
  --target-crack-extension-um 5 \
  --crystal-aniso --crystal-compete --no-crystal-branch --max-fronts 1 \
  --crack-backend adaptive_czm \
  --mpz-length-um "$MPZ_LENGTH_UM" --mpz-n-bins "$MPZ_N_BINS" \
  --save-snapshots 5 --snapshot-cols 5 \
  --snapshot-by-crack-extension-um 5 \
  --out "$OUTROOT"

"$PYTHON_BIN" audit_v10_0_no_artificial_controls.py \
  "$OUTROOT" --require-progressive
"$PYTHON_BIN" audit_v10_0_2_progressive_smoke.py \
  "$OUTROOT" --target-um 5

cat <<EOF
V10.0.2 PROGRESSIVE ONE-SEGMENT GATE PASSED
material=$MATERIAL
T_K=$T_K
out=$OUTROOT
The next gate is cohesive-penalty convergence; do not launch a longer or multi-temperature run yet.
EOF
