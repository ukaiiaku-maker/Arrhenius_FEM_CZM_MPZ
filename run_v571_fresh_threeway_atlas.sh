#!/usr/bin/env bash
set -euo pipefail

FRESH_ATLAS_OUT="${FRESH_ATLAS_OUT:-runs/sn_v1_barrier_phenomena_map_v5_6_fresh}"
EXTENSION_OUT="${EXTENSION_OUT:-runs/sn_v1_barrier_phenomena_extension_v5_7_1_fresh}"
FORCE_FRESH="${FORCE_FRESH:-0}"

if [[ -e "${FRESH_ATLAS_OUT}" && "${FORCE_FRESH}" != "1" ]]; then
  echo "ERROR: ${FRESH_ATLAS_OUT} already exists." >&2
  echo "Use the incremental extension workflow, or set FORCE_FRESH=1 deliberately." >&2
  exit 2
fi

N_SURFACES="${N_SURFACES:-3840}" \
CANDIDATE_POOL="${CANDIDATE_POOL:-2048}" \
DESIGN_BATCH_SIZE="${DESIGN_BATCH_SIZE:-256}" \
CYCLES_MAX="${CYCLES_MAX:-1e12}" \
FATIGUE_TEMPS="${FATIGUE_TEMPS:-100 200 300 400 500 600 700}" \
OUT="${FRESH_ATLAS_OUT}" \
bash run_sn_v1_representative_exp_floor_map.sh

ATLAS_DIR="${FRESH_ATLAS_OUT}" \
OUT="${EXTENSION_OUT}" \
bash run_v571_extend_existing_atlas.sh
