#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
OUT="${OUT:-runs/v1_two_barrier_dbtt_fatigue_map_corrected}"
SCOPE="${SCOPE:-core}"
TEMPS="${TEMPS:-300 400 500 600 700 900}"
KLIST="${KLIST:-0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.5 5.5 6.5 7.5 8.5 10 12 14}"
CYCLES_MAX="${CYCLES_MAX:-2e14}"
MAX_BLOCKS="${MAX_BLOCKS:-10000}"
MONO_KMAX="${MONO_KMAX:-40}"
CASE_FILTER="${CASE_FILTER:-}"
RESUME="${RESUME:-1}"

ARGS=(
  run_v1_two_barrier_dbtt_fatigue_map_corrected.py
  --case-table selected_v1_temperature_cases_corrected.csv
  --out "$OUT"
  --scope "$SCOPE"
  --temperatures $TEMPS
  --Kmax-MPa-sqrt-m $KLIST
  --cycles-max "$CYCLES_MAX"
  --max-blocks "$MAX_BLOCKS"
  --monotonic-Kmax-MPa "$MONO_KMAX"
)

if [[ "$RESUME" == "1" || "$RESUME" == "true" ]]; then
  ARGS+=(--resume)
fi

if [[ -n "$CASE_FILTER" ]]; then
  # shellcheck disable=SC2206
  FILTER_WORDS=( ${CASE_FILTER//,/ } )
  ARGS+=(--case-filter "${FILTER_WORDS[@]}")
fi

"$PYTHON_BIN" "${ARGS[@]}"
