#!/usr/bin/env bash
set -euo pipefail

# Small probe for entropy/plasticity-controlled process-zone evolution.
# Increase cycles-max and Kmax after the smoke test behaves sensibly.
python -m arrhenius_fracture.fatigue_sharp_front \
  --temperatures 300 400 500 600 700 800 900 1000 \
  --Kmax-MPa-sqrt-m 18 \
  --R 0.1 \
  --frequency-Hz 1000 \
  --cycles-max 1e8 \
  --block-cycles 1e5 \
  --max-block-cycles 1e6 \
  --target-dB 0.10 \
  --target-dN-store 0.10 \
  --n-phase 128 \
  --storage-model escape_limited \
  --out runs/fatigue_v1_temperature_probe
