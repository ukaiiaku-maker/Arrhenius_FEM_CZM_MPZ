#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-runs/v8_fatigue_v1_tuning_smoke}
python -m arrhenius_fracture.fatigue_sharp_front \
  --temperatures 300 500 700 \
  --Kmax-MPa-sqrt-m 20 \
  --R 0.1 \
  --frequency-Hz 1000 \
  --cycles-max 1e4 \
  --block-cycles 10 \
  --max-blocks 2000 \
  --min-block-cycles 1e-6 \
  --storage-model escape_limited \
  --dN-cap inf \
  --no-plots \
  --out "$OUT"
