#!/usr/bin/env bash
set -euo pipefail

OUT="${OUT:-runs/v1_two_barrier_dbtt_fatigue_map}"
SCOPE="${SCOPE:-core}"
TEMPS="${TEMPS:-300 400 500 600 700 900}"
KLIST="${KLIST:-3.5 4.5 5.5 6.5 7.5 8.5 10 12 14}"
CYCLES_MAX="${CYCLES_MAX:-2e14}"
MAX_BLOCKS="${MAX_BLOCKS:-5000}"
RESUME="${RESUME:-1}"

python -m compileall -q arrhenius_fracture

RESUME_ARG=""
if [[ "$RESUME" == "1" || "$RESUME" == "true" ]]; then RESUME_ARG="--resume"; fi

python run_v1_two_barrier_dbtt_fatigue_map.py \
  --case-table selected_v1_temperature_cases.csv \
  --out "${OUT}" \
  --scope "${SCOPE}" \
  --temperatures ${TEMPS} \
  --Kmax-MPa-sqrt-m ${KLIST} \
  --cycles-max "${CYCLES_MAX}" \
  --max-blocks "${MAX_BLOCKS}" \
  ${RESUME_ARG}
