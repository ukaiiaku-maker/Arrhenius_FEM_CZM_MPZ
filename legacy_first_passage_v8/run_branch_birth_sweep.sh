#!/usr/bin/env bash
set -euo pipefail

# Branch-birth / competing-hazard sweep for the sharp-tip model.
# Purpose: map whether a second near-critical cleavage hazard exists at initiation.
# This is the quantitative first-passage/Kc sweep. It may sever quickly after branching;
# use run_branch_morphology_sweep.sh for visual morphology-resolved traces.

ROOT=${1:-branch_birth_sweep_$(date +%Y%m%d_%H%M)}
mkdir -p "$ROOT"

# Shared resolved 2-D settings. Keep embrittlement cap OFF; use physical ledger saturation.
COMMON=(
  --mode 2d
  --nx 50 --ny 100
  --tip-h-fine 0.6e-6 --tip-ratio 1.25
  --n-stagger 2
  --save-snapshots 4 --snapshot-cols 4
  --print-every 20
  --crystal-aniso --crystal-compete --crystal-branch
  --crystal-material branchy
  --branch-share-mode hazard
  --branch-hazard-sharpness 2.0
  --branch-energy-share hazard-budget
  --emit-S-T-c0-kB=-20
  --emit-S-T-c1=0.02
  --emit-S-sigma-max-kB=8
  --multihit-m 3 --multihit-tau 1e-6
  --emb-sat-frac 1
  --n-sat 2000
  --cleave-shield-chi 0.2
  --cleave-H0-eV 2.6
  --v-rayleigh 2600
)

# Core sweep axes. Start here; expand only after confirming cost.
TEMPS=(500 700 850 900 950 1100)
THETAS=(0 15 30 37.5 45 52.5 60)
GAMMAS=(1.0 2.0 3.0)
RATIOS=(0.55 0.70 0.85)

# Rate-matched loading family. dU/dt ~ constant so the effective K-rate is comparable.
# dU=1e-6 gives better K resolution than the 4e-6 smoke test; increase steps for high-K cases.
DU=1e-6
DT=42
STEPS=140

for gamma in "${GAMMAS[@]}"; do
  for ratio in "${RATIOS[@]}"; do
    for theta in "${THETAS[@]}"; do
      tag="g${gamma}_r${ratio}_th${theta}"
      out="$ROOT/$tag"
      echo "=== running $tag ==="
      python3 -m arrhenius_fracture.sharp_front "${COMMON[@]}" \
        --temperatures "${TEMPS[@]}" \
        --crystal-theta-deg "$theta" \
        --cleave-gamma-aniso "$gamma" \
        --branch-overdrive-ratio "$ratio" \
        --dU "$DU" --dt "$DT" --steps "$STEPS" \
        --out "$out" \
        > "$out.log" 2>&1
    done
  done
 done

python3 analyze_branch_sweep.py "$ROOT" --out "$ROOT/branch_sweep_summary.csv"
echo "Done. Summary: $ROOT/branch_sweep_summary.csv"
