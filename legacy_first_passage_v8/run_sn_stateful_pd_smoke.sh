#!/usr/bin/env bash
set -euo pipefail

python -m compileall -q arrhenius_fracture

OUT="${OUT:-runs/sn_stateful_pd_smoke}"
python -m arrhenius_fracture.sn_pd2d_stateful \
  --out "$OUT" \
  --cases "${CASE:-no_shield}" \
  --T "${T:-300}" --R "${R:-0.1}" --frequency-Hz "${FREQ:-1000}" \
  --sigma-a-MPa ${STRESSES:-600} \
  --cycles-max "${CYCLES_MAX:-3e4}" \
  --block-cycles "${BLOCK_CYCLES:-1e4}" \
  --max-blocks "${MAX_BLOCKS:-8}" \
  --nx "${NX:-14}" --ny "${NY:-28}" \
  --root-h-fine "${ROOT_H_FINE:-45e-6}" \
  --pd-horizon-m "${PD_HORIZON:-150e-6}" \
  --pd-patch-radius-m "${PD_PATCH_RADIUS:-0.40e-3}" \
  --pd-boundary-shell-m "${PD_BOUNDARY_SHELL:-130e-6}" \
  --nu-grow-s "${NU_GROW_S:-2}" \
  --nu-link-s "${NU_LINK_S:-1}" \
  --grow-stress-GPa "${GROW_STRESS_GPA:-0.8}" \
  --link-stress-GPa "${LINK_STRESS_GPA:-0.8}" \
  --plastic-n-phase "${PLASTIC_PHASES:-6}" \
  --hazard-n-phase "${HAZARD_PHASES:-8}" \
  --snapshot-every "${SNAPSHOT_EVERY:-5}" \
  --print-every 1
