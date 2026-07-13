#!/usr/bin/env bash
set -euo pipefail

python -m arrhenius_fracture.fatigue_sharp_front \
  --temperatures 300 500 700 \
  --Kmax-MPa-sqrt-m 20 \
  --R 0.1 \
  --frequency-Hz 1000 \
  --cycles-max 1e6 \
  --block-cycles 1e4 \
  --max-block-cycles 1e5 \
  --target-dB 0.10 \
  --target-dN-store 0.10 \
  --n-phase 96 \
  --print-every 25 \
  --out runs/fatigue_v1_smoke
