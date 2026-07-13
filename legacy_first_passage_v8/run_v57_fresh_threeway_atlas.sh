#!/usr/bin/env bash
set -euo pipefail

# Fresh end-to-end run.  This is intentionally separate from the incremental
# extension wrapper so an existing V5.6 atlas is never overwritten by accident.
BASE_OUT="${BASE_OUT:-runs/sn_v1_barrier_phenomena_map_v5_7_fresh_base}"
EXT_OUT="${EXT_OUT:-runs/sn_v1_barrier_phenomena_map_v5_7_fresh_extension}"
N_SURFACES="${N_SURFACES:-3840}"
CANDIDATE_POOL="${CANDIDATE_POOL:-2048}"
DESIGN_BATCH_SIZE="${DESIGN_BATCH_SIZE:-256}"
CYCLES_MAX="${CYCLES_MAX:-1e12}"
FATIGUE_TEMPS="${FATIGUE_TEMPS:-100 200 300 400 500 600 700 800 900}"
FRACTURE_TEMPS="${FRACTURE_TEMPS:-100 200 300 400 500 600 700 800 900}"
FORCE_FRESH="${FORCE_FRESH:-0}"

export PYTHONNOUSERSITE=1

if [[ -e "${BASE_OUT}/independent_exp_floor_design_v5_6.csv" || -e "${EXT_OUT}/fracture_monotonic_points_v5_7.csv" ]]; then
  if [[ "${FORCE_FRESH}" == "1" ]]; then
    rm -rf "${BASE_OUT}" "${EXT_OUT}"
  else
    echo "ERROR: fresh-run output already exists. Set FORCE_FRESH=1 to delete it explicitly." >&2
    exit 2
  fi
fi

N_SURFACES="${N_SURFACES}" \
CANDIDATE_POOL="${CANDIDATE_POOL}" \
DESIGN_BATCH_SIZE="${DESIGN_BATCH_SIZE}" \
CYCLES_MAX="${CYCLES_MAX}" \
FATIGUE_TEMPS="${FATIGUE_TEMPS}" \
OUT="${BASE_OUT}" \
bash run_sn_v1_representative_exp_floor_map.sh

ATLAS_DIR="${BASE_OUT}" \
OUT="${EXT_OUT}" \
FRACTURE_TEMPS="${FRACTURE_TEMPS}" \
bash run_v57_extend_existing_atlas.sh
