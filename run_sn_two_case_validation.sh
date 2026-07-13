#!/usr/bin/env bash
set -euo pipefail

# Validation run designed to expose shielding separation without launching a
# full-resolution production atlas.  Unlike the 10^4-cycle smoke test, this
# run allows the plastic state to evolve enough for shielding feedback to act.

python -m compileall -q arrhenius_fracture

V1_OUT="${V1_OUT:-runs/sn_v1_two_case_validation}"
PF2D_OUT="${PF2D_OUT:-runs/sn_pf2d_two_case_validation}"
T="${T:-300}"
R="${R:-0.1}"
FREQ="${FREQ:-1000}"

python -m arrhenius_fracture.sn_v1 \
  --out "${V1_OUT}" \
  --T "${T}" --R "${R}" --frequency-Hz "${FREQ}" \
  --sigma-a-MPa ${V1_STRESSES:-250 300 350 400 450 500 550 600 650 700 750 800 900} \
  --cycles-max "${V1_CYCLES_MAX:-1e10}" \
  --max-blocks "${V1_MAX_BLOCKS:-5000}"

# Small-mesh field validation at the amplitude where V1 predicts the strongest
# case separation.  Default horizon is long enough to build P~O(0.1-0.4).
python -m arrhenius_fracture.sn_pf2d \
  --out "${PF2D_OUT}" \
  --T "${T}" --R "${R}" --frequency-Hz "${FREQ}" \
  --sigma-a-MPa ${PF2D_STRESSES:-600} \
  --cycles-max "${PF2D_CYCLES_MAX:-1e8}" \
  --block-cycles "${PF2D_BLOCK_CYCLES:-1e6}" \
  --max-blocks "${PF2D_MAX_BLOCKS:-120}" \
  --nx "${NX:-18}" --ny "${NY:-36}" \
  --root-h-fine "${ROOT_H_FINE:-30e-6}" \
  --snapshot-every "${SNAPSHOT_EVERY:-1000}" \
  --print-every "${PRINT_EVERY:-20}"

python postprocess_sn_two_case.py --v1-root "${V1_OUT}" --pf2d-root "${PF2D_OUT}"
