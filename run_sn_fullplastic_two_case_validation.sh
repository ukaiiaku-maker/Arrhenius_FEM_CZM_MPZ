#!/usr/bin/env bash
set -euo pipefail

python -m compileall -q arrhenius_fracture

OUT="${OUT:-runs/sn_pf2d_fullplastic_two_case_validation}"
T="${T:-300}"
R="${R:-0.1}"
FREQ="${FREQ:-1000}"
STRESSES="${STRESSES:-600 700 800}"

python -m arrhenius_fracture.sn_pf2d_fullplastic \
  --out "${OUT}" \
  --T "${T}" --R "${R}" --frequency-Hz "${FREQ}" \
  --sigma-a-MPa ${STRESSES} \
  --cycles-max "${CYCLES_MAX:-1e9}" \
  --block-cycles "${BLOCK_CYCLES:-1e8}" \
  --max-blocks "${MAX_BLOCKS:-1500}" \
  --target-dep-eq-block "${TARGET_DEP_BLOCK:-2e-4}" \
  --target-rho-rel-block "${TARGET_RHO_REL:-0.05}" \
  --target-dB-nuc "${TARGET_DB_NUC:-0.05}" \
  --nx "${NX:-36}" --ny "${NY:-72}" \
  --root-h-fine "${ROOT_H_FINE:-30e-6}" \
  --plastic-n-phase "${PLASTIC_PHASES:-12}" \
  --hazard-n-phase "${HAZARD_PHASES:-16}" \
  --morph-band-length-m "${MORPH_BAND_LENGTH:-100e-6}" \
  --morph-decay-length-m "${MORPH_DECAY_LENGTH:-150e-6}" \
  --snapshot-every "${SNAPSHOT_EVERY:-50}" \
  --print-every "${PRINT_EVERY:-10}"

python postprocess_sn_fullplastic.py --root "${OUT}"
