#!/usr/bin/env bash
set -euo pipefail
python -m compileall -q arrhenius_fracture

OUT="${OUT:-runs/sn_v1_arrhenius_two_case}"
python -m arrhenius_fracture.sn_v1_arrhenius \
  --out "${OUT}" \
  --T "${T:-300}" --R "${R:-0.1}" --frequency-Hz "${FREQ:-1000}" \
  --sigma-a-MPa ${STRESSES:-250 300 350 400 450 500 550 600 650 700 750 800 900} \
  --cycles-max "${CYCLES_MAX:-1e10}" \
  --block-cycles "${BLOCK_CYCLES:-1e8}" \
  --max-blocks "${MAX_BLOCKS:-5000}" \
  --n-phase "${N_PHASE:-64}" \
  --target-dep-eq-block "${TARGET_DEP_BLOCK:-2e-4}" \
  --target-rho-rel-block "${TARGET_RHO_REL:-0.05}" \
  --target-dB-nuc "${TARGET_DB_NUC:-0.05}" \
  --exp-system "${EXP_SYSTEM:-W[100]}" \
  --emit-energy-scale "${EMIT_E_SCALE:-0.75}" \
  --emit-entropy-scale "${EMIT_S_SCALE:-0.75}" \
  --peierls-energy-scale "${PEIERLS_E_SCALE:-0.00375}" \
  --peierls-entropy-scale "${PEIERLS_S_SCALE:-0.00375}" \
  --taylor-energy-scale "${TAYLOR_E_SCALE:-0.015}" \
  --taylor-entropy-scale "${TAYLOR_S_SCALE:-0.015}"

python postprocess_sn_v1_arrhenius_two_case.py --root "${OUT}"
