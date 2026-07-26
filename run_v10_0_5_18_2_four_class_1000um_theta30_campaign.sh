#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}

export OPTIONS=${OPTIONS:-"v913_paper_peak01_0242980_persistent_sites v913_paper_dbtt01_0202500_persistent_sites v913_paper_weakT01_0129902_persistent_sites v913_paper_ceramic01_0077080_persistent_sites"}
export TEMPERATURES=${TEMPERATURES:-"300 400 500 600 700 800 850 900 950 1000 1050 1100 1150 1200"}
export THETA=${THETA:-30}
export TARGET_EXT_UM=${TARGET_EXT_UM:-1000}
export STEPS=${STEPS:-2000000}
export MAX_JOBS=${MAX_JOBS:-2}
export SKIP_FINISHED=${SKIP_FINISHED:-1}
export PRINT_EVERY=${PRINT_EVERY:-100}
export GENERATE_SOLVER_PLOTS=${GENERATE_SOLVER_PLOTS:-0}
export SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-10}
export SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-100}
export SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
export PYTHONUNBUFFERED=${PYTHONUNBUFFERED:-1}
export CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_18_2_four_class_14T_1000um_theta30_v1}

exec bash "$ROOT/run_v10_0_5_18_2_four_class_focused_100um.sh"
