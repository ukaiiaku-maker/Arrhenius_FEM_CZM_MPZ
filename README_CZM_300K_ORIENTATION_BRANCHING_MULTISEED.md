# 300 K FEM/CZM multi-seed orientation + weak-branching R-curve sweep v3

This v3 package is a correction to v2. The v2 runner tried to avoid topology problems by over-refining
the local crack-tip mesh. In practice that pushed the 300 K ceramic orientation case into repeated
`local_hrefine_error:degenerate elements` loops and the monitor did not catch them because it only
recognized two-element veto signatures. v3 changes both parts:

- Uses the previously robust coarse/stable FEM/CZM settings by default:
  - `nx=12`, `ny=24`
  - `tip_h_fine=5e-6`
  - `tip_ratio=1.3`
  - `da_phys=5e-6`
  - `adaptive_event_target=0.35`
  - `adaptive_min_frac=1e-8`
  - `adaptive_grow=4.0`
  - `czm_max_angle_error_deg=60`
- Turns off the slow/finer retry by default.
- Turns off the wall-clock timeout by default.
- Detects repeated geometry vetoes with either one-element or two-element signatures:
  - `[759]`
  - `[1916, 1964]`
- Kills a bad topology trap after 8 repeated veto signatures rather than waiting hours.

## Recommended orientation run

```bash
CONDA_ENV=arrhenius-fem-czm \
CLASS=ceramic \
T_K=300 \
THETAS="0 15 30 45" \
RUN_BRANCH=0 \
SEEDS="1201 1202 1203 1204 1205" \
TARGET_EXT_UM=1000 \
OUTROOT=runs/czm_Rcurve_300K_orientation_ceramic_5seed_stable_v3 \
bash run_czm_300K_orientation_branching_multiseed.sh
```

## Recommended weak-branching run

```bash
CONDA_ENV=arrhenius-fem-czm \
CLASS=ceramic \
T_K=300 \
THETAS="" \
RUN_BRANCH=1 \
BRANCH_THETA=30 \
SEEDS="1201 1202 1203 1204 1205" \
TARGET_EXT_UM=1000 \
OUTROOT=runs/czm_Rcurve_300K_branching_ceramic_theta30_5seed_stable_v3 \
bash run_czm_300K_orientation_branching_multiseed.sh
```

## Plot

```bash
PYTHON_BIN=/opt/homebrew/Caskroom/miniconda/base/envs/arrhenius-fem-czm/bin/python \
ROOT=runs/czm_Rcurve_300K_orientation_ceramic_5seed_stable_v3 \
bash run_plot_czm_300K_orientation_branching_multiseed_Rcurves.sh
```

This is still a stochastic topology calculation. v3 is designed to avoid letting a single meshing trap run
for hours; it will mark the case incomplete and move to the next seed/condition.
