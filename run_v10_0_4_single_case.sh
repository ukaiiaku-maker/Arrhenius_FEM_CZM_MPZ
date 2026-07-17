#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
MATERIAL=${MATERIAL:-weakT}
T_K=${T_K:-700}
THETA_DEG=${THETA_DEG:-45}
OUTROOT=${OUTROOT:?Set OUTROOT to a new versioned output directory}
CZM_PENALTY_NORMAL=${CZM_PENALTY_NORMAL:-1e18}
CZM_PENALTY_TANGENT=${CZM_PENALTY_TANGENT:-1e18}
RUN_TESTS=${RUN_TESTS:-1}
NO_PLOTS=${NO_PLOTS:-1}
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

if [[ "$THETA_DEG" != "45" && "$THETA_DEG" != "45.0" ]]; then
  echo "ERROR: v10.0.4 completion gates are certified only for THETA_DEG=45"
  exit 1
fi

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

if [[ -e "$OUTROOT" ]]; then
  echo "ERROR: output path already exists: $OUTROOT"
  exit 1
fi

if [[ "$RUN_TESTS" == "1" ]]; then
  CONDA_ENV="$CONDA_ENV" PYTHON_BIN="$PYTHON_BIN" \
    bash run_v10_0_3_integration_tests.sh
fi

EXTRA_ARGS=()
if [[ "$NO_PLOTS" == "1" ]]; then
  EXTRA_ARGS+=(--no-plots --save-snapshots 0)
else
  EXTRA_ARGS+=(--save-snapshots 5 --snapshot-cols 5 --snapshot-by-crack-extension-um 5)
fi

ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM=5 \
ARRHENIUS_PREFINED_MODE_I_CORRIDOR=1 \
ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY="$MIN_TRIANGLE_QUALITY" \
ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO="$MIN_CHILD_AREA_RATIO" \
ARRHENIUS_MAX_TIP_H_OVER_DA="$MAX_TIP_H_OVER_DA" \
ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE="$MAX_TRIAL_DAMAGE_CHANGE" \
ARRHENIUS_MIN_TRIAL_RETRY_DT_S="$MIN_TRIAL_RETRY_DT_S" \
ARRHENIUS_MAX_TRIAL_RETRIES="$MAX_TRIAL_RETRIES" \
ARRHENIUS_MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL="$MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL" \
"$PYTHON_BIN" -m arrhenius_fracture.mode_i_first_passage_v10_0_3_progressive \
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
  --czm-penalty-normal "$CZM_PENALTY_NORMAL" \
  --czm-penalty-tangent "$CZM_PENALTY_TANGENT" \
  --mpz-length-um "$MPZ_LENGTH_UM" --mpz-n-bins "$MPZ_N_BINS" \
  "${EXTRA_ARGS[@]}" \
  --out "$OUTROOT"

"$PYTHON_BIN" audit_v10_0_3_progressive_integration.py \
  "$OUTROOT" --target-um 5

"$PYTHON_BIN" normalize_v10_0_3_1_reporting.py "$OUTROOT"

OUTROOT="$OUTROOT" MATERIAL="$MATERIAL" T_K="$T_K" THETA_DEG="$THETA_DEG" \
CZM_PENALTY_NORMAL="$CZM_PENALTY_NORMAL" \
CZM_PENALTY_TANGENT="$CZM_PENALTY_TANGENT" \
"$PYTHON_BIN" - <<'PY'
import json
import os
import subprocess
from pathlib import Path

root = Path(os.environ["OUTROOT"])
cert = json.loads((root / "v10_0_3_progressive_integration_certification.json").read_text())
summary = json.loads((root / "anisotropic_calibrated_tip_first_passage_summary.json").read_text())
results = json.loads((root / "mode_i_v10_0_3_1_results.json").read_text())
try:
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
except Exception:
    commit = None
payload = {
    "schema": "v10_0_4_single_case_metadata",
    "orchestration_point_release": "10.0.4",
    "integration_kernel": "10.0.3",
    "reporting_layer": "10.0.3.1",
    "git_commit": commit,
    "material": os.environ["MATERIAL"],
    "temperature_K": float(os.environ["T_K"]),
    "crystal_theta_deg": float(os.environ["THETA_DEG"]),
    "czm_penalty_normal_Pa_per_m": float(os.environ["CZM_PENALTY_NORMAL"]),
    "czm_penalty_tangent_Pa_per_m": float(os.environ["CZM_PENALTY_TANGENT"]),
    "target_extension_um": 5.0,
    "certified": bool(cert.get("certified", False)),
    "certification_file": "v10_0_3_progressive_integration_certification.json",
    "parameter_fingerprint_sha256": summary.get("parameter_fingerprint_sha256"),
    "front_state_model": summary.get("front_state_model"),
    "B_final": summary.get("B_final"),
    "source_population_bound": summary.get("source_population_bound"),
    "result": results[0] if results else None,
}
(root / "v10_0_4_case_metadata.json").write_text(json.dumps(payload, indent=2, default=str))
print("V10.0.4 SINGLE CASE CERTIFIED")
print(json.dumps({k: payload[k] for k in (
    "material", "temperature_K", "czm_penalty_normal_Pa_per_m",
    "czm_penalty_tangent_Pa_per_m", "certified")}, indent=2))
PY
