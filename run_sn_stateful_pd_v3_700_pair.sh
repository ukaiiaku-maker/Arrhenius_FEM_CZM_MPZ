#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
OUT="${OUT:-runs/sn_stateful_pd_v3p1_700_pair}" \
STRESSES="700" \
CASES="no_shield shielded" \
MESH_SEED="${MESH_SEED:-1}" \
PD_SEED="${PD_SEED:-1}" \
SNAPSHOT_EVERY="${SNAPSHOT_EVERY:-25}" \
PRINT_EVERY="${PRINT_EVERY:-5}" \
bash run_sn_stateful_pd_pilot.sh
