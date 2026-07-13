#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-clusterJ_v8_smoke}
mkdir -p "$OUT"
python3 -m arrhenius_fracture.sharp_front \
  --mode 2d --nx 18 --ny 36 --tip-h-fine 0 --tip-ratio 1.25 \
  --dU 2e-7 --dt 8.4 --steps 5 --n-stagger 1 \
  --save-snapshots 2 --snapshot-cols 2 --print-every 1 \
  --crystal-aniso --crystal-compete --crystal-branch --crystal-material branchy \
  --cleave-gamma-aniso 2.0 --branch-overdrive-ratio 0.95 \
  --branch-fp-min-ratio 0.98 --branch-clock-target 1.0 \
  --branch-secondary-min-K-ratio 0.90 --branch-secondary-min-K-MPa 10.0 \
  --branch-starved-suppression-radius 120e-6 --branch-spacing 50 \
  --retire-stagnant-branches --branch-stagnant-lag 120e-6 \
  --branch-stagnant-K-ratio 0.25 --branch-stagnant-lambda 1e-25 \
  --branch-stagnant-steps 25 --branch-stagnant-no-fire-steps 80 \
  --branch-share-mode hazard --branch-hazard-sharpness 2.0 --branch-energy-share hazard-budget \
  --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
  --multihit-m 3 --multihit-tau 1e-6 --cleave-H0-eV 2.6 --cleave-shield-chi 0.2 \
  --emb-sat-frac 1 --n-sat 2000 --adaptive-events --adaptive-event-target 0.35 \
  --adaptive-min-frac 1e-8 --adaptive-grow 4.0 --max-fronts 32 \
  --da-phys 2e-6 --j-decomposition cluster --rJ-cluster 10e-6 --rJ-outer 12e-6 \
  --branch-resolve-length 20e-6 --local-J-clearance-factor 1.0 --min-J-active-elems 12 \
  --local-J-handoff-min-K-ratio 0.60 --temperatures 300 --crystal-theta-deg 30 --out "$OUT"
