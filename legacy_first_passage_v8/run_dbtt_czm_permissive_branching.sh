#!/usr/bin/env bash
set -euo pipefail

# Long-growth DBTT/CZM run with a deliberately permissive branch-birth preset.
# The underlying driver remains the single source of physics parameters.
export LONG_GROWTH=1
export PERMISSIVE_BRANCHING=1
export ENABLE_BRANCHING=1
export BRANCH_FP_MIN_RATIO="${BRANCH_FP_MIN_RATIO:-0.82}"
export BRANCH_CLOCK_TARGET="${BRANCH_CLOCK_TARGET:-0.75}"
export BRANCH_SECONDARY_MIN_K_RATIO="${BRANCH_SECONDARY_MIN_K_RATIO:-0.80}"
export BRANCH_SPACING="${BRANCH_SPACING:-8.0}"
export TARGET_EXT_UM="${TARGET_EXT_UM:-750}"
export LONG_STEPS="${LONG_STEPS:-20000}"
export SAVE_SNAPSHOTS="${SAVE_SNAPSHOTS:-12}"
export SNAPSHOT_COLS="${SNAPSHOT_COLS:-6}"
export SNAPSHOT_BY_EXT_UM="${SNAPSHOT_BY_EXT_UM:-75}"
export MAX_JOBS="${MAX_JOBS:-1}"
export MAX_FRONTS="${MAX_FRONTS:-8}"
export RETIRE_STAGNANT_BRANCHES="${RETIRE_STAGNANT_BRANCHES:-1}"
export THETA="${THETA:-45}"
export TEMPS="${TEMPS:-900}"
export OUTROOT="${OUTROOT:-runs/dbtt_czm_theta${THETA}_permissive_branching}"

exec bash run_dbtt_czm_orientation_temperature_test.sh
