#!/usr/bin/env bash
set -euo pipefail

python -m compileall -q arrhenius_fracture

OUT="${OUT:-runs/sn_pf2d_arrhenius_fullplastic_production}"
T="${T:-300}"
R="${R:-0.1}"
FREQ="${FREQ:-1000}"
STRESSES="${STRESSES:-500 550 600 700}"

python -m arrhenius_fracture.sn_pf2d_fullplastic \
  --out "${OUT}" \
  --T "${T}" --R "${R}" --frequency-Hz "${FREQ}" \
  --sigma-a-MPa ${STRESSES} \
  --cycles-max "${CYCLES_MAX:-1e10}" \
  --block-cycles "${BLOCK_CYCLES:-1e8}" \
  --max-blocks "${MAX_BLOCKS:-3000}" \
  --target-dep-eq-block "${TARGET_DEP_BLOCK:-2e-4}" \
  --target-rho-rel-block "${TARGET_RHO_REL:-0.05}" \
  --target-dB-nuc "${TARGET_DB_NUC:-0.05}" \
  --nx "${NX:-36}" --ny "${NY:-72}" \
  --root-h-fine "${ROOT_H_FINE:-30e-6}" \
  --plastic-n-phase "${PLASTIC_PHASES:-12}" \
  --hazard-n-phase "${HAZARD_PHASES:-16}" \
  --exp-system "${EXP_SYSTEM:-W[100]}" \
  --emit-energy-scale "${EMIT_E_SCALE:-0.75}" \
  --emit-entropy-scale "${EMIT_S_SCALE:-0.75}" \
  --peierls-energy-scale "${PEIERLS_E_SCALE:-0.00375}" \
  --peierls-entropy-scale "${PEIERLS_S_SCALE:-0.00375}" \
  --taylor-energy-scale "${TAYLOR_E_SCALE:-0.015}" \
  --taylor-entropy-scale "${TAYLOR_S_SCALE:-0.015}" \
  --plastic-event-strain "${PLASTIC_EVENT_STRAIN:-1e-5}" \
  --phi-taylor-max "${PHI_TAYLOR_MAX:-20}" \
  --morph-band-length-m "${MORPH_BAND_LENGTH:-100e-6}" \
  --morph-decay-length-m "${MORPH_DECAY_LENGTH:-150e-6}" \
  --snapshot-every "${SNAPSHOT_EVERY:-50}" \
  --print-every "${PRINT_EVERY:-10}"

python postprocess_sn_fullplastic.py --root "${OUT}"
