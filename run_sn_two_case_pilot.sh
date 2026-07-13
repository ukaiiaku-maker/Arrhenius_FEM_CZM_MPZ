#!/usr/bin/env bash
set -euo pipefail

# Two-case S-N pilot: reduced V1 sweep first, then selected full-field 2-D points.
# Both cases use the same barrier family. Only the plastic-state -> crack-opening
# shielding feedback differs.

python -m compileall -q arrhenius_fracture

V1_OUT="${V1_OUT:-runs/sn_v1_two_case}"
PF2D_OUT="${PF2D_OUT:-runs/sn_pf2d_two_case}"
T="${T:-300}"
R="${R:-0.1}"
FREQ="${FREQ:-1000}"

python -m arrhenius_fracture.sn_v1 \
  --out "${V1_OUT}" \
  --T "${T}" --R "${R}" --frequency-Hz "${FREQ}" \
  --sigma-a-MPa ${V1_STRESSES:-250 300 350 400 450 500 550 600 650 700 750 800 900} \
  --cycles-max "${V1_CYCLES_MAX:-1e10}" \
  --max-blocks "${V1_MAX_BLOCKS:-5000}"

# The 2-D pilot is intentionally smaller. Start with three amplitudes spanning
# the V1 transition and expand after the field evolution is inspected.
python -m arrhenius_fracture.sn_pf2d \
  --out "${PF2D_OUT}" \
  --T "${T}" --R "${R}" --frequency-Hz "${FREQ}" \
  --sigma-a-MPa ${PF2D_STRESSES:-500 600 700} \
  --cycles-max "${PF2D_CYCLES_MAX:-1e8}" \
  --max-blocks "${PF2D_MAX_BLOCKS:-1500}" \
  --nx "${NX:-36}" --ny "${NY:-72}" \
  --root-h-fine "${ROOT_H_FINE:-30e-6}" \
  --snapshot-every "${SNAPSHOT_EVERY:-100}"

python postprocess_sn_two_case.py --v1-root "${V1_OUT}" --pf2d-root "${PF2D_OUT}"
