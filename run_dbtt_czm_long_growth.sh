#!/usr/bin/env bash
set -euo pipefail

# Convenience wrapper for propagation/morphology runs of the anisotropic
# Arrhenius FEM/CZM DBTT case.  The underlying temperature driver remains the
# single source of physics parameters; this wrapper only selects long-growth
# stopping and snapshot policies.

export LONG_GROWTH=1
export TARGET_EXT_UM="${TARGET_EXT_UM:-1400}"
export LONG_STEPS="${LONG_STEPS:-20000}"
export SAVE_SNAPSHOTS="${SAVE_SNAPSHOTS:-12}"
export SNAPSHOT_COLS="${SNAPSHOT_COLS:-6}"
export SNAPSHOT_BY_EXT_UM="${SNAPSHOT_BY_EXT_UM:-100}"
export MAX_JOBS="${MAX_JOBS:-1}"
export THETA="${THETA:-30}"
export TEMPS="${TEMPS:-300 400 500 600 700 800 900 1000}"
export OUTROOT="${OUTROOT:-runs/dbtt_czm_theta${THETA}_long_growth}"

exec bash run_dbtt_czm_orientation_temperature_test.sh
