# 2-D field run: the regime physics in the spatial FEM model

The regime map (regime_map.py) is the 1-D single-point engine ramped in K. The
2-D mode (`--mode 2d`, run_2d) is the full plane-strain FEM: it assembles
mechanics, runs the Arrhenius-Taylor plasticity update with density transport,
gets K from a J-integral around the moving tip, and advances the crack (kills
elements) on the engine's cleavage events, depositing the emitted wake density
into the rho field. It builds the engine from the same args, so chi_shield,
emb_sat_frac, n_sat, recover_k, v_rayleigh all propagate.

## Run (process-zone-resolved, rate-matched to the 1-D regime)
```
python3 -m arrhenius_fracture.sharp_front --mode 2d \
  --temperatures 800 950 --nx 50 --ny 100 --tip-h-fine 0.6e-6 --tip-ratio 1.25 \
  --steps 130 --dU 2e-6 --dt 84 --n-stagger 2 --save-snapshots 6 \
  --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
  --multihit-m 3 --multihit-tau 1e-6 \
  --cleave-H0-eV 3.0 --cleave-shield-chi 0.6 --emb-sat-frac 1 --n-sat 2000 \
  --v-rayleigh 2600 --out run2d
```

## Two requirements for a CREDIBLE 2-D run (both easy to get wrong)
1. RESOLVE THE PROCESS ZONE. Use a tip-graded mesh with tip_h_fine < L_pz
   (here hbar_tip=0.86 um < L_pz=1 um). On a uniform coarse mesh the back stress
   and J-contour are mesh-floored and the shielding physics is lost.
2. MATCH THE LOADING RATE. K_J is driven by dU/dt; dK_J/dt must equal the 1-D
   Kdot or the rate-dependent transition shifts. Here dt=84 s gives
   dK_J/dt ~ 0.005 MPa*sqrt(m)/s, matching the regime map. (At dt=1 the run is
   ~70x too fast and BOTH temperatures come out brittle.)

## Result (reproduces the 1-D DBTT from the fields)
| T    | Kc (J-integral) | N_em at fracture | peak eq. plastic strain | reading |
|------|-----------------|------------------|-------------------------|---------|
| 800K | 14.1 MPa*sqrt(m)| ~340             | ~0.017                  | brittle lower shelf |
| 950K | 37.8 MPa*sqrt(m)| ~1800            | ~0.05                   | tough upper shelf   |

1-D map at the same point predicted 13 and 41 -> the 2-D fields reproduce the
transition (the ~10% offset is real field/J-integral vs analytic-tip difference).
At 950 K the crack is pinned at a0 while the emission-driven dislocation field and
plastic zone grow and shield the tip, releasing only at K_J~38.

## Fields output per temperature (field_snapshots_<T>K.png, 4 columns = evolution)
 - damage d (the crack), log10 rho (dislocation density), sigma1 (max principal
   stress), eq. plastic strain. Plus diagnostics: energetics, dislocations,
   hazard_clocks, load_displacement, tip_state, toughness; and steps_<T>K.csv
   (full per-step history). summary.json + toughness_vs_temperature.png.

## Caveats / next refinements for the 2-D fields
 - Legacy BULK multiplication (default) makes rho domain-filling at high T (net-
   section yield of the small specimen). For a tip-LOCALIZED dislocation field use
   sources-only plasticity: --bulk-mult-frac 0 --tip-source-rho-per-emit <s>
   --rho-transport-c <c> (optionally --exhaustion). This ties the bulk field to the
   emission rate, consistent with the shielding/embrittlement story.
 - 2-D-specific clamps to audit before quantitative 2-D claims (per AUDIT_*.md):
   rho_cap (transport clip), mobile_rho_floor, peierls_floor, and the da_phys/mesh
   advance increment. None affect the 1-D regime conclusions.
 - Boundary stress/strain concentrations at the loaded edges are mesh/BC artifacts;
   a larger domain or SSY boundary layer would clean them up.
 - Post-initiation the crack severs in one load step (fast fracture); v_rayleigh
   bounds the velocity but at dt=84 s the per-step cap is non-binding. For a graded
   post-initiation R-curve, reduce dt after initiation (adaptive stepping).
