#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
MATERIAL=${MATERIAL:-DBTT}
T_K=${T_K:-700}
THETA_DEG=${THETA_DEG:-45}
TARGET_UM=${TARGET_UM:-100}
DA_UM=${DA_UM:-5}
OUTROOT=${OUTROOT:-runs/v10_0_5_2_DBTT_700K_theta45_100um_200bins_v1}
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
TENSOR_DRIVE_PROBE_RADIUS_M=${TENSOR_DRIVE_PROBE_RADIUS_M:-1e-5}
TENSOR_DRIVE_SECTOR_HALF_ANGLE_DEG=${TENSOR_DRIVE_SECTOR_HALF_ANGLE_DEG:-25}
TENSOR_DRIVE_MIN_ELEMENTS=${TENSOR_DRIVE_MIN_ELEMENTS:-3}
SNAPSHOT_INTERVAL_UM=${SNAPSHOT_INTERVAL_UM:-10}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-11}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

if [[ "$TARGET_UM" != "100" ]]; then
  echo "ERROR: this guarded runner is certified only for TARGET_UM=100"
  exit 1
fi
if [[ "$DA_UM" != "5" ]]; then
  echo "ERROR: this guarded runner requires DA_UM=5"
  exit 1
fi
if [[ "$MPZ_N_BINS" != "200" ]]; then
  echo "ERROR: this guarded runner requires MPZ_N_BINS=200"
  exit 1
fi
if [[ -e "$OUTROOT" ]]; then
  echo "ERROR: output path already exists: $OUTROOT"
  echo "Choose a new versioned OUTROOT; this runner never reuses a prior directory."
  exit 1
fi

CONDA_ENV="$CONDA_ENV" PYTHON_BIN="$PYTHON_BIN" \
  bash run_v10_0_3_integration_tests.sh

ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM="$TARGET_UM" \
ARRHENIUS_PREFINED_MODE_I_CORRIDOR=1 \
ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY="$MIN_TRIANGLE_QUALITY" \
ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO="$MIN_CHILD_AREA_RATIO" \
ARRHENIUS_MAX_TIP_H_OVER_DA="$MAX_TIP_H_OVER_DA" \
ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE="$MAX_TRIAL_DAMAGE_CHANGE" \
ARRHENIUS_MIN_TRIAL_RETRY_DT_S="$MIN_TRIAL_RETRY_DT_S" \
ARRHENIUS_MAX_TRIAL_RETRIES="$MAX_TRIAL_RETRIES" \
ARRHENIUS_MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL="$MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL" \
ARRHENIUS_TENSOR_DRIVE_PROBE_RADIUS_M="$TENSOR_DRIVE_PROBE_RADIUS_M" \
ARRHENIUS_TENSOR_DRIVE_SECTOR_HALF_ANGLE_DEG="$TENSOR_DRIVE_SECTOR_HALF_ANGLE_DEG" \
ARRHENIUS_TENSOR_DRIVE_MIN_ELEMENTS="$TENSOR_DRIVE_MIN_ELEMENTS" \
"$PYTHON_BIN" -m arrhenius_fracture.mode_i_first_passage_v10_0_5_2_parallel \
  --v10-material-class "$MATERIAL" \
  --czm-opening-coupling clock_linear \
  --mode 2d \
  --temperatures "$T_K" \
  --steps "$STEPS" \
  --nx "$NX" --ny "$NY" \
  --tip-h-fine "$TIP_H_FINE" --tip-ratio "$TIP_RATIO" \
  --dU "$DU" --dt "$DT" \
  --da-phys "${DA_UM}e-6" \
  --target-crack-extension-um "$TARGET_UM" \
  --crystal-aniso --crystal-compete \
  --crystal-theta-deg "$THETA_DEG" \
  --max-fronts 1 \
  --crack-backend adaptive_czm \
  --mpz-length-um "$MPZ_LENGTH_UM" --mpz-n-bins "$MPZ_N_BINS" \
  --save-snapshots "$SAVE_SNAPSHOTS" --snapshot-cols 6 \
  --snapshot-by-crack-extension-um "$SNAPSHOT_INTERVAL_UM" \
  --out "$OUTROOT"

"$PYTHON_BIN" audit_v10_0_3_progressive_integration.py \
  "$OUTROOT" --target-um "$TARGET_UM"

"$PYTHON_BIN" normalize_v10_0_3_1_reporting.py "$OUTROOT"
"$PYTHON_BIN" normalize_v10_0_5_1_slip_trace_reporting.py "$OUTROOT"

"$PYTHON_BIN" audit_v10_0_5_2_long_growth.py \
  "$OUTROOT" \
  --target-um "$TARGET_UM" \
  --expected-mpz-bins "$MPZ_N_BINS" \
  --da-um "$DA_UM"

cat <<EOF
V10.0.5.2 DBTT 700 K 100 UM MULTICOMMIT GATE PASSED
out=$OUTROOT
material=$MATERIAL
T_K=$T_K
theta_deg=$THETA_DEG
target_um=$TARGET_UM
physical_increment_um=$DA_UM
expected_commits=20
mpz_length_um=$MPZ_LENGTH_UM
mpz_n_bins=$MPZ_N_BINS
parallel_coupling=v10.0.5
slip_trace_reporting=v10.0.5.1
complete_channel_diagnostics=v10.0.5.2
No material-response classification or reparameterization criterion was applied.
The 500 um production run remains unauthorized until this gate is reviewed.
EOF
