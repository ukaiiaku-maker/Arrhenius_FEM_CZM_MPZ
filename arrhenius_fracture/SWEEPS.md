# Sweep & convergence command set

All runs from the package parent dir. Fixed emission kinetics for comparability:
`--emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 --multihit-m 3 --multihit-tau 1e-6`.
Always run regime studies with the artificial barrier cap OFF (`--emb-sat-frac 1`)
and bound the ledger physically with `--n-sat` (or `--recover-k`).

Shorthand:
```
KIN="--emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 --multihit-m 3 --multihit-tau 1e-6"
RAMP="--Kdot 0.005 --Kmax 40 --dt 1.0"
TS="--temperatures 500 600 700 800 900 1000 1100 1200 1300"
```

## 1. One command per named regime (physical, no caps)
```
# ceramic  (no recovery: embrittlement wins high-T)
python3 -m arrhenius_fracture.sharp_front --mode 1d $TS $RAMP $KIN \
    --cleave-H0-eV 2.6 --cleave-shield-chi 0.0 --emb-sat-frac 1 --n-sat inf --out out/ceramic

# toughness peak  (no recovery + mild shielding)
python3 -m arrhenius_fracture.sharp_front --mode 1d $TS $RAMP $KIN \
    --cleave-H0-eV 3.6 --cleave-shield-chi 0.10 --emb-sat-frac 1 --n-sat inf --out out/peak

# weakly T-dependent  (recovery + light shielding tuned to cancel intrinsic softening)
python3 -m arrhenius_fracture.sharp_front --mode 1d $TS $RAMP $KIN \
    --cleave-H0-eV 4.0 --cleave-shield-chi 0.20 --emb-sat-frac 1 --n-sat 1500 --out out/weakT

# DBTT  (recovery bounds dG_emb -> tough shelf; shielding sets onset)
python3 -m arrhenius_fracture.sharp_front --mode 1d $TS $RAMP $KIN \
    --cleave-H0-eV 6.0 --cleave-shield-chi 0.60 --emb-sat-frac 1 --n-sat 2000 --out out/dbtt
```

## 2. Axis sweeps (bash loops)
```
# (a) shielding axis: onset / low-T toughening
for CHI in 0.0 0.2 0.4 0.6 0.9 1.3; do
  python3 -m arrhenius_fracture.sharp_front --mode 1d $TS $RAMP $KIN \
    --cleave-H0-eV 4.0 --emb-sat-frac 1 --n-sat 2000 --cleave-shield-chi $CHI \
    --out out/chi_$CHI ; done

# (b) recovery axis: shelf survival (inf = no recovery -> ceramic)
for NS in inf 8000 3000 1500 800; do
  python3 -m arrhenius_fracture.sharp_front --mode 1d $TS $RAMP $KIN \
    --cleave-H0-eV 4.0 --emb-sat-frac 1 --cleave-shield-chi 0.4 --n-sat $NS \
    --out out/nsat_$NS ; done

# (b') recovery via annihilation rate instead of a hard ceiling
for K in 0.0 0.01 0.03 0.05 0.1; do
  python3 -m arrhenius_fracture.sharp_front --mode 1d $TS $RAMP $KIN \
    --cleave-H0-eV 6.0 --emb-sat-frac 1 --cleave-shield-chi 0.6 --recover-k $K \
    --out out/reck_$K ; done

# (c) cleavage barrier: shelf height/position
for H0 in 2.6 3.4 4.2 5.0 6.0; do
  python3 -m arrhenius_fracture.sharp_front --mode 1d $TS $RAMP $KIN \
    --emb-sat-frac 1 --n-sat 2000 --cleave-shield-chi 0.6 --cleave-H0-eV $H0 \
    --out out/H0_$H0 ; done

# (d) emission onset temperature (entropy crossover)
for TC0 in -10 -20 -30; do for TC1 in 0.01 0.02 0.04; do
  python3 -m arrhenius_fracture.sharp_front --mode 1d $TS $RAMP --multihit-m 3 --multihit-tau 1e-6 \
    --emit-S-T-c0-kB=$TC0 --emit-S-T-c1=$TC1 --emit-S-sigma-max-kB=8 \
    --cleave-H0-eV 6.0 --emb-sat-frac 1 --n-sat 2000 --cleave-shield-chi 0.6 \
    --out out/onset_${TC0}_${TC1} ; done; done

# (e) loading rate (kinetic competition; shifts transition T)
for KD in 0.001 0.005 0.02 0.1; do
  python3 -m arrhenius_fracture.sharp_front --mode 1d $TS --Kdot $KD --Kmax 40 --dt 1.0 $KIN \
    --cleave-H0-eV 6.0 --emb-sat-frac 1 --n-sat 2000 --cleave-shield-chi 0.6 \
    --out out/kdot_$KD ; done
```

## 3. 2D phase map + representative curves (one figure)
```
REGIME_OUT=./regime_map python3 -m arrhenius_fracture.regime_map
# edits inside regime_map.py select the chi x n_sat grid and H0_B.
```

## 4. Convergence / robustness (must pass before trusting a regime point)
```
python3 -m arrhenius_fracture.convergence
```
Checks, at a physical DBTT point (chi=0.6, H0=6, n_sat=2e3, cap off):
 - dt in {4,2,1,0.5,0.25}  -> Kc must be flat (observed drift 0.2% over 16x).
 - dt with dN_cap=1e9       -> isolate clock from emission throttle (identical).
 - dN_cap in {25..1e9}      -> Kc flat (with recovery, throttle is irrelevant).
 - Kmax in {30,40,60,80}    -> arrests are genuine, not window censoring (identical).
 - tau_c in {1e-7,1e-6,1e-5}-> MODEL sensitivity ~5%/decade (calibrate, don't "converge").

Additional convergence to run when changing the tip model:
 - r0 / L_pz (process-zone radius and pile-up length): rescale together and confirm
   the dimensionless regime boundaries (chi vs n_sat) are invariant.
 - da (advance increment, 1D) and the 2D mesh hbar: Kc at first advance should be
   insensitive (first passage is set by B>=1, not by da).
