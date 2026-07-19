# Sharp-front branching patch notes

This patch keeps the model sharp-interface only.  No smeared variational-fracture criterion is introduced.  The goal is to make branching enter through competing sharp-tip hazards rather than through smeared damage topology.

## What changed

1. **Continuous branch advance**
   - Removed the effective `n_fire // 2` branch split.
   - Branches now share a continuous advance length.  A normal first-passage event with `n_fire = 1` can therefore move both the parent and daughter fronts.
   - `--branch-energy-share hazard-budget` keeps a finite shared length budget when both fronts fire in the same step.  `--branch-energy-share none` lets each front advance independently after birth.

2. **Conservative daughter-front ledger inheritance**
   - Added `FrontEngine.clone_split(daughter_fraction)`.
   - At branch birth, the parent tip ledger is split rather than duplicated or reset:
     - `N_em`
     - cleavage renewal clock `B`
     - emission work ledger `W_emit`
   - The daughter gets its own `FrontEngine` and then evolves independently under its local `K_J`, local direction, and local hazards.

3. **Branch birth from competing hazard statistics**
   - Branch birth share is computed from the two competing directional overdrives using
     `share_i ∝ overdrive_i^branch_hazard_sharpness`.
   - `--branch-share-mode hazard` is the default.
   - `--branch-share-mode equal` is available as a diagnostic.

4. **Branch direction is preserved on the birth step**
   - The daughter front now advances along the secondary lobe that caused the branch.
   - This prevents the identical stress state at the common origin from immediately collapsing both fronts back onto the primary winner.

5. **Wake and source deposition are front-local**
   - Primary wake density is deposited at the primary tip.
   - Daughter wake density is deposited at the daughter tip.
   - If `--tip-source-rho-per-emit` is used, source density is also deposited at the local front that generated it.

6. **Branch diagnostics**
   - Each 2-D temperature run writes `branch_diagnostics_<T>K.csv` with:
     - number of competing candidates
     - primary and secondary angles
     - primary and secondary overdrive/metric
     - `metric2/metric1`
     - branch active/spawned flags
     - branch shares
     - primary/daughter `n_fire`
     - actual front advances
     - primary/daughter coordinates
     - primary/daughter cleavage hazard rates
   - Each run also writes `run_args.json` so no-branch/branch outcomes can be audited against the exact command line.

7. **Second material class for branch-prone sweeps**
   - Added `--crystal-material branchy`.
   - This is a model material class, not tungsten.  It exists to stress-test branching with stronger anisotropy and more branch-prone co-critical hazards.
   - Branchy preset defaults:
     - `C44 = 320 GPa` if not explicitly provided
     - `--cleave-gamma-aniso 2.0`
     - `--branch-overdrive-ratio 0.80`
     - `--branch-ratio 0.85`
     - `--crystal-include-110`
     - `--gamma-110-rel 1.15`
   - The W/default preset remains near-isotropic and conservative:
     - `--cleave-gamma-aniso 0.3`
     - `--branch-overdrive-ratio 0.9`
     - `--branch-ratio 0.92`
     - no secondary `{110}` cleavage traces unless requested.

## Important interpretation

This patch makes branching possible in the sharp-front code, but it does not make branching automatic.  Branching still requires the local directional hazard/overdrive landscape to become multi-lobed or near-co-critical.  The default W-like material is expected to deflect more often than it branches.  The `branchy` material is deliberately more anisotropic and should be used as a controlled second material class for branch-statistics sweeps.

## Smoke-test command

```bash
python3 -m arrhenius_fracture.sharp_front --mode 2d \
  --temperatures 900 \
  --steps 6 --nx 20 --ny 40 \
  --tip-h-fine 1e-6 --tip-ratio 1.4 \
  --dU 4e-6 --dt 84 --n-stagger 1 \
  --print-every 1 --save-snapshots 1 \
  --crystal-aniso --crystal-compete --crystal-branch \
  --crystal-material branchy --crystal-theta-deg 45 \
  --branch-overdrive-ratio 0.4 \
  --branch-energy-share hazard-budget \
  --cleave-H0-eV 2.6 --cleave-shield-chi 0.2 \
  --n-sat 2000 --emb-sat-frac 1 \
  --out test_branch2d
```

Expected qualitative result: a `BRANCH at ...` message and two crack-path files with opposite-sign y-deflections:

- `crack_path_900K.csv`
- `crack_path_branch_900K.csv`

This is only a smoke test.  It intentionally uses a coarse mesh and aggressive forcing.  Use the resolved sweep commands below for physics comparisons.
