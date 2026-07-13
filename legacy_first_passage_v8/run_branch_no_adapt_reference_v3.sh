#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-branch_no_adapt_reference_v3}
mkdir -p "$OUT"
# Reference run: no adaptive stepping and no caps. This may sever the ligament in
# one step if the tip is genuinely unstable at the accepted FEM/J state.
python3 -m arrhenius_fracture.sharp_front --mode 2d \
  --temperatures 900 \
  --steps 20 --nx 50 --ny 100 \
  --tip-h-fine 0.6e-6 --tip-ratio 1.25 \
  --dU 2e-7 --dt 8.4 --n-stagger 2 \
  --save-snapshots 4 --snapshot-cols 4 --print-every 1 \
  --crystal-aniso --crystal-compete --crystal-branch \
  --crystal-material branchy --crystal-theta-deg 45 \
  --cleave-gamma-aniso 2.0 --branch-overdrive-ratio 0.70 \
  --branch-share-mode hazard --branch-hazard-sharpness 2.0 \
  --branch-energy-share hazard-budget \
  --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
  --multihit-m 3 --multihit-tau 1e-6 \
  --cleave-H0-eV 2.6 --cleave-shield-chi 0.2 \
  --emb-sat-frac 1 --n-sat 2000 \
  --out "$OUT/ref_T900_th45"
