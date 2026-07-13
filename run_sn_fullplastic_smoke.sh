#!/usr/bin/env bash
set -euo pipefail
OUT="${OUT:-runs/sn_pf2d_fullplastic_smoke}"
python -m arrhenius_fracture.sn_pf2d_fullplastic \
  --out "${OUT}" --sigma-a-MPa 700 \
  --cycles-max 1e7 --block-cycles 1e7 --max-blocks 3 \
  --nx 12 --ny 24 --root-h-fine 80e-6 \
  --plastic-n-phase 6 --hazard-n-phase 4 \
  --snapshot-every 0 --print-every 1
python postprocess_sn_fullplastic.py --root "${OUT}"
