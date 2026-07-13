#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
OUTROOT="${OUTROOT:-runs/mixed_mode_fem_czm_v2_event_controlled_500K}"
CLASSES="${CLASSES:-ceramic DBTT}"
TARGET_PSI="${TARGET_PSI:--60 -45 -30 -15 0 15 30 45 60}"
INTERVAL="${MONITOR_INTERVAL:-15}"
STALE_MINUTES="${STALE_MINUTES:-15}"
SHOW_ACTIVE_LOG_LINES="${SHOW_ACTIVE_LOG_LINES:-0}"

exec conda run --no-capture-output -n "$CONDA_ENV" \
  python -u monitor_mixed_mode_fem_czm_v2_1.py \
  --root "$OUTROOT" \
  --classes "$CLASSES" \
  --target-psi-deg="$TARGET_PSI" \
  --interval "$INTERVAL" \
  --stale-minutes "$STALE_MINUTES" \
  --show-active-log-lines "$SHOW_ACTIVE_LOG_LINES" \
  "$@"
