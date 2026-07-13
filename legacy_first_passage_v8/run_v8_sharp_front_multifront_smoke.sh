#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-runs/v8_sharp_front_multifront_smoke}
python -m arrhenius_fracture.sharp_front \
  --mode 2d \
  --temperatures 300 \
  --steps 5 \
  --nx 24 --ny 48 \
  --tip-h-fine 1.0e-6 \
  --tip-ratio 1.25 \
  --n-stagger 1 \
  --crystal-aniso \
  --crystal-branch \
  --j-decomposition cluster \
  --adaptive-events \
  --sigma-cap-GPa 0 \
  --dN-cap inf \
  --emb-sat-frac 1 \
  --out "$OUT"
