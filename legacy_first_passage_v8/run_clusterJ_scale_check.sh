#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-clusterJ_scale_check}
mkdir -p "$OUT"
BASE="--mode 2d --nx 50 --ny 100 \
  --tip-h-fine 0.45e-6 --tip-ratio 1.25 \
  --dU 2e-7 --dt 8.4 --steps 1200 --n-stagger 2 \
  --save-snapshots 4 --snapshot-cols 4 --print-every 200 \
  --crystal-aniso --crystal-compete --crystal-branch --crystal-material branchy \
  --cleave-gamma-aniso 2.0 --branch-overdrive-ratio 0.70 \
  --branch-share-mode hazard --branch-hazard-sharpness 2.0 \
  --branch-energy-share hazard-budget \
  --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
  --multihit-m 3 --multihit-tau 1e-6 \
  --cleave-H0-eV 2.6 --cleave-shield-chi 0.2 \
  --emb-sat-frac 1 --n-sat 2000 \
  --adaptive-events --adaptive-event-target 0.35 --adaptive-min-frac 1e-8 \
  --adaptive-grow 4.0 --max-fronts 16 --branch-spacing 10 \
  --da-phys 2e-6 --j-decomposition cluster \
  --local-J-clearance-factor 1.0 --min-J-active-elems 12 \
  --local-J-handoff-min-K-ratio 0.25"
for RJCL in 5e-6 10e-6 15e-6; do
  for RJLOC in 8e-6 12e-6; do
    RES=$(python3 - <<PY
r=float('$RJLOC')
print(f'{max(20e-6, 1.7*r):.6g}')
PY
)
    TAG="T300_th45_rJcl${RJCL}_rJloc${RJLOC}"
    TAG=${TAG//./p}
    echo "=== $OUT/$TAG ==="
    python3 -m arrhenius_fracture.sharp_front $BASE \
      --temperatures 300 --crystal-theta-deg 45 \
      --rJ-cluster "$RJCL" --rJ-outer "$RJLOC" --branch-resolve-length "$RES" \
      --out "$OUT/$TAG"
  done
done
python3 analyze_clusterJ_growth.py "$OUT" > "$OUT/clusterJ_summary.csv" || true
