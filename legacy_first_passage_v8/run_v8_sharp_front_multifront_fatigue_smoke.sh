#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-runs/v8_sharp_front_multifront_fatigue_smoke}
python -m arrhenius_fracture.sharp_front \
  --mode 2d \
  --fatigue-cycles \
  --temperatures 300 \
  --steps 5 \
  --nx 24 --ny 48 \
  --tip-h-fine 1.0e-6 \
  --tip-ratio 1.25 \
  --n-stagger 1 \
  --crystal-aniso \
  --crystal-branch \
  --j-decomposition cluster \
  --sigma-cap-GPa 0 \
  --dN-cap inf \
  --R 0.1 \
  --frequency-Hz 1000 \
  --block-cycles 10 \
  --min-block-cycles 1e-6 \
  --target-dB 0.2 \
  --target-dN-store 0.25 \
  --storage-model escape_limited \
  --out "$OUT"
