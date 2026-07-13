#!/usr/bin/env bash
set -euo pipefail

# Three-rate production campaign built on the no-branch 500-um R-curve sweep.
# Rate is changed by holding nominal dU fixed and scaling dt by 1/rate_factor.
# The sharp-front adaptive controller then scales dU and dt together by the
# same trial fraction, preserving the requested rate under step refinement.

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PARAMETERS=${PARAMETERS:-four_class_exp_floor_exact_model_inputs.csv}
OUTROOT=${OUTROOT:-runs/four_class_exp_floor_CZM_rates_no_branch_500um_theta45}
CLASSES=${CLASSES:-"ceramic peak weakT DBTT"}
TEMPS=${TEMPS:-"300 400 500 600 700 800 900 1000 1100 1200"}
THETA=${THETA:-45}
TARGET_EXT_UM=${TARGET_EXT_UM:-500}
RATE_FACTORS=${RATE_FACTORS:-"1 10 100"}
BASE_DU=${BASE_DU:-2e-7}
BASE_DT=${BASE_DT:-8.4}
LONG_STEPS=${LONG_STEPS:-20000}
MAX_JOBS=${MAX_JOBS:-1}
FORCE=${FORCE:-0}

# Snapshot defaults: initial state plus extension-triggered snapshots through
# the 500-um propagation history.  The solver writes field_snapshots_<T>K.png.
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-12}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-6}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-50}

# A cheap three-case rate preflight before the 120-case production campaign.
PREFLIGHT_RATE_SMOKE=${PREFLIGHT_RATE_SMOKE:-1}
PREFLIGHT_CLASS=${PREFLIGHT_CLASS:-weakT}
PREFLIGHT_TEMP=${PREFLIGHT_TEMP:-900}
PREFLIGHT_EXT_UM=${PREFLIGHT_EXT_UM:-25}

if [[ ! -f run_four_class_exp_floor_czm_500um_sweep.py ]]; then
  echo "ERROR: run_four_class_exp_floor_czm_500um_sweep.py not found." >&2
  echo "Apply the four-class 500-um R-curve sweep patch/package first." >&2
  exit 2
fi
if [[ ! -f "$PARAMETERS" ]]; then
  echo "ERROR: parameter table not found: $PARAMETERS" >&2
  exit 2
fi

if [[ -n "${PYTHON_BIN:-}" ]]; then
  PY="$PYTHON_BIN"
else
  if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found; set PYTHON_BIN explicitly" >&2
    exit 2
  fi
  PY="$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' 2>&1 | tr -d '\r' | awk 'NF {last=$0} END {print last}')"
fi
if [[ -z "$PY" || ! -x "$PY" ]]; then
  echo "ERROR: could not resolve Python interpreter" >&2
  exit 2
fi

mkdir -p "$OUTROOT"

# Write an auditable campaign map including the nominal engineering strain-rate
# ratios for the default 4-mm specimen height.
OUTROOT="$OUTROOT" RATE_FACTORS="$RATE_FACTORS" BASE_DU="$BASE_DU" BASE_DT="$BASE_DT" \
SAVE_SNAPSHOTS="$SAVE_SNAPSHOTS" SNAPSHOT_COLS="$SNAPSHOT_COLS" SNAPSHOT_BY_EXT_UM="$SNAPSHOT_BY_EXT_UM" \
"$PY" - <<'PY'
import json, os
from pathlib import Path
root=Path(os.environ['OUTROOT'])
factors=[float(x) for x in os.environ['RATE_FACTORS'].replace(',',' ').split()]
dU=float(os.environ['BASE_DU']); dt0=float(os.environ['BASE_DT']); Ly=4e-3
rows=[]
for f in factors:
    dt=dt0/f
    rows.append({
        'rate_factor': f,
        'dU_nominal_m': dU,
        'dt_nominal_s': dt,
        'opening_rate_m_per_s': dU/dt,
        'engineering_strain_rate_per_s': (dU/dt)/Ly,
    })
(root/'rate_campaign_config.json').write_text(json.dumps({
    'rate_definition': 'fixed nominal dU; nominal dt=BASE_DT/rate_factor',
    'adaptive_rule': 'accepted dU and dt are multiplied by the same adaptive_frac',
    'specimen_height_m': Ly,
    'rates': rows,
    'snapshots': {
        'save_snapshots': int(os.environ.get('SAVE_SNAPSHOTS','12')),
        'snapshot_cols': int(os.environ.get('SNAPSHOT_COLS','6')),
        'snapshot_by_extension_um': float(os.environ.get('SNAPSHOT_BY_EXT_UM','50')),
    },
}, indent=2))
PY

run_one_rate() {
  local factor="$1"
  local dt
  dt="$($PY - <<PY
print(${BASE_DT}/float(${factor}))
PY
)"
  local label="rate_${factor}x"
  local root="$OUTROOT/$label"
  echo
  echo "========================================================================"
  echo " RATE ${factor}x | dU=${BASE_DU} m | dt=${dt} s | out=${root}"
  echo "========================================================================"

  local args=(
    --parameters "$PARAMETERS"
    --outroot "$root"
    --classes "$CLASSES"
    --temps "$TEMPS"
    --theta "$THETA"
    --target-ext-um "$TARGET_EXT_UM"
    --long-steps "$LONG_STEPS"
    --max-jobs "$MAX_JOBS"
    --conda-env "$CONDA_ENV"
    --python-bin "$PY"
    --dU "$BASE_DU"
    --dt "$dt"
    --save-snapshots "$SAVE_SNAPSHOTS"
    --snapshot-cols "$SNAPSHOT_COLS"
    --snapshot-by-ext-um "$SNAPSHOT_BY_EXT_UM"
  )
  [[ "$FORCE" == "1" ]] && args+=(--force)
  "$PY" run_four_class_exp_floor_czm_500um_sweep.py "${args[@]}"
}

run_preflight() {
  local smoke_root="$OUTROOT/_rate_preflight"
  rm -rf "$smoke_root"
  echo "=== Three-rate adaptive timestep preflight: ${PREFLIGHT_CLASS}, ${PREFLIGHT_TEMP} K, ${PREFLIGHT_EXT_UM} um ==="
  for factor in $RATE_FACTORS; do
    local dt
    dt="$($PY - <<PY
print(${BASE_DT}/float(${factor}))
PY
)"
    "$PY" run_four_class_exp_floor_czm_500um_sweep.py \
      --parameters "$PARAMETERS" \
      --outroot "$smoke_root/rate_${factor}x" \
      --classes "$PREFLIGHT_CLASS" \
      --temps "$PREFLIGHT_TEMP" \
      --theta "$THETA" \
      --target-ext-um "$PREFLIGHT_EXT_UM" \
      --long-steps "$LONG_STEPS" \
      --max-jobs 1 \
      --conda-env "$CONDA_ENV" \
      --python-bin "$PY" \
      --dU "$BASE_DU" \
      --dt "$dt" \
      --save-snapshots 3 \
      --snapshot-cols 3 \
      --snapshot-by-ext-um 10
  done
  "$PY" audit_four_class_rate_sweep.py \
    --root "$smoke_root" \
    --rate-factors "$RATE_FACTORS" \
    --base-dU "$BASE_DU" \
    --base-dt "$BASE_DT"
  echo "Rate preflight passed."
}

if [[ "$PREFLIGHT_RATE_SMOKE" == "1" ]]; then
  run_preflight
fi

FAILED=0
for factor in $RATE_FACTORS; do
  if ! run_one_rate "$factor"; then
    echo "WARNING: ${factor}x campaign had failed or incomplete cases; continuing." >&2
    FAILED=1
  fi
done

# Audit all completed cases even if one campaign returned nonzero.
if ! "$PY" audit_four_class_rate_sweep.py \
  --root "$OUTROOT" \
  --rate-factors "$RATE_FACTORS" \
  --base-dU "$BASE_DU" \
  --base-dt "$BASE_DT"; then
  echo "WARNING: adaptive timestep audit reported one or more problems." >&2
  FAILED=1
fi

# Combine rate summaries and make direct 1x/10x/100x comparison figures.
if ! "$PY" summarize_four_class_rate_sweep.py \
  --root "$OUTROOT" \
  --rate-factors "$RATE_FACTORS"; then
  echo "WARNING: rate-comparison summary/plotting failed." >&2
  FAILED=1
fi

if [[ "$FAILED" == "1" ]]; then
  exit 1
fi

echo "=== All requested rate campaigns completed and passed timestep audit ==="
