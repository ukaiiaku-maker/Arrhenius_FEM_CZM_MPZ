# Suggested sharp-front branch/deflection sweeps

All commands use the sharp-interface driver only.  They are intended to separate three questions:

1. W-like material: does the crack deflect but usually remain single-front?
2. Branch-prone model material: does a second material class produce co-critical hazards and branching?
3. Regime coupling: does branching change with T, shielding/recovery, and loading rate?

## Shared resolved 2-D flags

```bash
COMMON="--mode 2d \
  --nx 50 --ny 100 \
  --tip-h-fine 0.6e-6 --tip-ratio 1.25 \
  --n-stagger 2 --save-snapshots 4 --snapshot-cols 4 \
  --crystal-aniso --crystal-compete --crystal-branch \
  --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
  --multihit-m 3 --multihit-tau 1e-6 \
  --emb-sat-frac 1 --n-sat 2000 --v-rayleigh 2600"
```

## 1. W-like baseline: deflection-dominated

```bash
for THETA in 0 15 30 45; do
  python3 -m arrhenius_fracture.sharp_front $COMMON \
    --crystal-material w \
    --crystal-theta-deg $THETA \
    --cleave-H0-eV 3.0 --cleave-shield-chi 0.6 \
    --temperatures 600 800 900 1000 1200 \
    --dU 2e-6 --dt 84 --steps 130 \
    --out run2d_W_theta${THETA}
done
```

Expected: deflection and path meander, but few or no true branches unless the local overdrive field becomes multi-lobed.

## 2. Branch-prone model material: anisotropy/material class sweep

```bash
for THETA in 30 45; do
  for RATIO in 0.70 0.80 0.90; do
    python3 -m arrhenius_fracture.sharp_front $COMMON \
      --crystal-material branchy \
      --crystal-theta-deg $THETA \
      --branch-overdrive-ratio $RATIO \
      --branch-hazard-sharpness 2 \
      --branch-energy-share hazard-budget \
      --cleave-H0-eV 3.0 --cleave-shield-chi 0.6 \
      --temperatures 800 900 1000 \
      --dU 2e-6 --dt 84 --steps 130 \
      --out run2d_branchy_theta${THETA}_ratio${RATIO}
  done
done
```

Expected: branch frequency should increase as `branch-overdrive-ratio` decreases and as the branch-prone material preset is used.  Use `branch_diagnostics_*K.csv` to verify whether branching was caused by true near-co-critical lobes rather than by the topology bookkeeping.

## 3. Branching across fracture regimes

Ceramic-like, no shielding/recovery:

```bash
python3 -m arrhenius_fracture.sharp_front $COMMON \
  --crystal-material branchy --crystal-theta-deg 45 \
  --cleave-H0-eV 2.6 --cleave-shield-chi 0.0 --n-sat inf \
  --temperatures 500 700 900 1100 1300 \
  --dU 5e-7 --dt 21 --steps 120 \
  --out run2d_branchy_ceramic
```

Toughness peak, mild shielding, unbounded ledger:

```bash
python3 -m arrhenius_fracture.sharp_front $COMMON \
  --crystal-material branchy --crystal-theta-deg 45 \
  --cleave-H0-eV 3.6 --cleave-shield-chi 0.10 --n-sat inf \
  --temperatures 600 800 900 1000 1200 \
  --dU 1e-6 --dt 42 --steps 120 \
  --out run2d_branchy_peak
```

Weakly temperature-dependent:

```bash
python3 -m arrhenius_fracture.sharp_front $COMMON \
  --crystal-material branchy --crystal-theta-deg 45 \
  --cleave-H0-eV 4.0 --cleave-shield-chi 0.20 --n-sat 1500 \
  --temperatures 500 700 900 1100 1300 \
  --dU 1e-6 --dt 42 --steps 90 \
  --out run2d_branchy_weakT
```

DBTT:

```bash
python3 -m arrhenius_fracture.sharp_front $COMMON \
  --crystal-material branchy --crystal-theta-deg 45 \
  --cleave-H0-eV 3.0 --cleave-shield-chi 0.60 --n-sat 2000 \
  --temperatures 600 800 900 1000 1200 \
  --dU 2e-6 --dt 84 --steps 130 \
  --out run2d_branchy_dbtt
```

## Branch diagnostics to inspect

For every temperature, inspect:

- `branch_diagnostics_<T>K.csv`
- `crack_path_<T>K.csv`
- `crack_path_branch_<T>K.csv`, when present
- `summary.json`

The diagnostic columns to watch first are:

- `n_candidates`
- `metric2_over_metric1`
- `branch_spawned`
- `share_w1`, `share_w2`
- `KJ1_Pa_sqrtm`, `KJ2_Pa_sqrtm`
- `lambda_c1`, `lambda_c2`
- `advance1_m`, `advance2_m`

A legitimate branch should have `n_candidates >= 2`, `metric2_over_metric1` near the chosen branch threshold, `branch_spawned = 1`, and nonzero daughter-front advance/path length.
