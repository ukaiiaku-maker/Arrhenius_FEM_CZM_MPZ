#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
MATERIAL=${MATERIAL:-weakT}
T_K=${T_K:-700}
THETA_DEG=${THETA_DEG:-45}
OUTROOT=${OUTROOT:-runs/v10_0_5_1_parallel_${MATERIAL}_${T_K}K_theta${THETA_DEG}_5um_v1}
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

CONDA_ENV="$CONDA_ENV" PYTHON_BIN="$PYTHON_BIN" \
  bash run_v10_0_3_integration_tests.sh

ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM=5 \
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
"$PYTHON_BIN" -m arrhenius_fracture.mode_i_first_passage_v10_0_5_parallel \
  --v10-material-class "$MATERIAL" \
  --czm-opening-coupling clock_linear \
  --mode 2d \
  --temperatures "$T_K" \
  --steps "$STEPS" \
  --nx "$NX" --ny "$NY" \
  --tip-h-fine "$TIP_H_FINE" --tip-ratio "$TIP_RATIO" \
  --dU "$DU" --dt "$DT" \
  --da-phys 5e-6 \
  --target-crack-extension-um 5 \
  --crystal-aniso --crystal-compete \
  --crystal-theta-deg "$THETA_DEG" \
  --max-fronts 1 \
  --crack-backend adaptive_czm \
  --mpz-length-um "$MPZ_LENGTH_UM" --mpz-n-bins "$MPZ_N_BINS" \
  --save-snapshots 5 --snapshot-cols 5 \
  --snapshot-by-crack-extension-um 5 \
  --out "$OUTROOT"

"$PYTHON_BIN" audit_v10_0_3_progressive_integration.py \
  "$OUTROOT" --target-um 5

"$PYTHON_BIN" normalize_v10_0_3_1_reporting.py "$OUTROOT"

"$PYTHON_BIN" - "$OUTROOT" <<'PY'
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
audit = json.loads(
    (root / "parallel_opening_emission_v10_0_5_audit.json").read_text()
)
assert audit["implementation_certified"] is True
assert audit["capture_count"] > 0
assert audit["nonzero_emission_drive_capture_count"] > 0
assert audit["drive_factor_normalization_or_clipping_active"] is False
assert audit["directional_multiplier_applied_after_hazard"] is False
assert audit["fit_derived_shielding_cap_active"] is False
assert (root / "mode_i_v10_0_5_results.json").is_file()
print("V10.0.5 PARALLEL COUPLING AUDIT VERIFIED")
print(json.dumps({
    "capture_count": audit["capture_count"],
    "nonzero_drive_captures": audit["nonzero_emission_drive_capture_count"],
    "min_drive_factor": audit["min_emission_drive_factor"],
    "max_drive_factor": audit["max_emission_drive_factor"],
    "emission_observation_required": False,
}, indent=2))
PY

"$PYTHON_BIN" normalize_v10_0_5_1_slip_trace_reporting.py "$OUTROOT"

cat <<EOF
V10.0.5.1 AUDITED REDUCED 2-D SLIP-TRACE COUPLING CHECKPOINT PASSED
material=$MATERIAL
T_K=$T_K
theta_deg=$THETA_DEG
out=$OUTROOT
integration_kernel=v10.0.3
reporting_normalization=v10.0.3.1
parallel_coupling=v10.0.5
slip_trace_reporting=v10.0.5.1
The plastic channels are reduced 2-D slip-trace channels, not a complete 3-D BCC slip-system model.
Nonzero emission is not required for implementation certification.
No material-response classification or reparameterization criterion was applied.
EOF
