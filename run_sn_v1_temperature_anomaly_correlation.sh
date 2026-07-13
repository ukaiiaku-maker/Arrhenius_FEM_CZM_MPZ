#!/usr/bin/env bash
set -euo pipefail
python -m compileall -q arrhenius_fracture

python run_sn_v1_temperature_anomaly_correlation.py \
  --out "${OUT:-runs/sn_v1_temperature_anomaly_correlation}" \
  --systems ${SYSTEMS:-W[100] Ta[111] Cu Al0.7CoCrFeNi-BCC Al0.7CoCrFeNi-FCC} \
  --entropy-multipliers ${ENTROPY_MULTIPLIERS:-0.0 0.5 1.0 1.5} \
  --temperatures ${TEMPS:-300 350 400 450 500 600 700 800 900 1000} \
  --sn-stresses ${SN_STRESSES:-200 250 300 350 400 450 500 550 600 650 700 750 800 900} \
  --strength-epsdot "${STRENGTH_EPSDOT:-1e-4}" \
  --cycles-max "${CYCLES_MAX:-1e10}" \
  --max-blocks "${MAX_BLOCKS:-5000}" \
  --sn-n-phase "${SN_N_PHASE:-48}" \
  --sn-block-cycles "${SN_BLOCK_CYCLES:-1e8}" \
  --sn-target-dep "${SN_TARGET_DEP:-5e-4}" \
  --sn-target-rho-rel "${SN_TARGET_RHO_REL:-0.10}" \
  --sn-target-dB "${SN_TARGET_DB:-0.10}" \
  --case "${CASE:-shielded}"
