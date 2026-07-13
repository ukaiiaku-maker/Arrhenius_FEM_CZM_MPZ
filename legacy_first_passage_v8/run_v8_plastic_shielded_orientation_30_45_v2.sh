#!/usr/bin/env bash
set -euo pipefail

# Focused crystallographic fatigue-path study using the current v8 production
# workflow and the plastic-shielded parameterization.
#
# Runs:
#   plastic_shielded_case64_M1
#   Kmax = 7.0 MPa sqrt(m)
#   theta = 30 deg
#   theta = 45 deg
#
# Crystal anisotropy and natural branching are already enabled in the v8 runner.
# This wrapper changes only crystal orientation.

PROD_DRIVER="${PROD_DRIVER:-run_v8_material_response_production_2d_oriented.sh}"

if [[ ! -f "${PROD_DRIVER}" ]]; then
  echo "ERROR: ${PROD_DRIVER} not found." >&2
  echo "Place all three files from the orientation package in the fatigue-PF project root." >&2
  exit 2
fi

if [[ ! -f "run_v8_compare_1d_2d_K_sweep_oriented.py" ]]; then
  echo "ERROR: run_v8_compare_1d_2d_K_sweep_oriented.py not found." >&2
  exit 2
fi

run_theta() {
  local theta="$1"
  local outroot="runs/v8_orientation_plastic_shielded_K7/theta${theta}"

  echo
  echo "======================================================================"
  echo "plastic_shielded_case64_M1 | Kmax=7.0 | theta=${theta} deg"
  echo "output: ${outroot}"
  echo "======================================================================"

  CASE_FILTER="plastic_shielded_case64_M1" \
  KLIST_OVERRIDE="7.0" \
  OUTROOT="${outroot}" \
  PRODUCTION_LEVEL="full" \
  TARGET_EXT_UM="${TARGET_EXT_UM:-1470}" \
  CYCLES_MAX="${CYCLES_MAX:-2e14}" \
  BLOCKS="${BLOCKS:-25000}" \
  SNAPSHOTS="${SNAPSHOTS:-24}" \
  SNAPSHOT_BY_EXT_UM="${SNAPSHOT_BY_EXT_UM:-25}" \
  MAKE_2D_PLOTS="${MAKE_2D_PLOTS:-1}" \
  CRYSTAL_THETA_DEG="${theta}" \
  bash "${PROD_DRIVER}"
}

run_theta 30
run_theta 45

echo
echo "Both orientation runs completed."
