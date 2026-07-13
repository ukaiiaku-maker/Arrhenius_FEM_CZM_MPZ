# Mixed-mode FEM/CZM v2: event-controlled deterministic first passage

## Why v2

The v1 campaign calibrated boundary loading angle only in a small elastic solve.  The achieved phase angle then drifted substantially by first passage or censoring.  v1 also used three exponential hazard seeds, which spent runs on threshold scatter rather than controlling geometry.

v2 makes two structural changes:

1. **Deterministic first passage:** cumulative cleavage threshold is fixed at `H=1`.  One run is used per control iteration; there is no default seed sweep.
2. **Closed-loop event-state control:** for every requested phase angle, the solver reruns with a corrected boundary angle until the projected phase at the actual first-passage state (or right-censored endpoint) is within tolerance.

## Mode decomposition

v2 replaces the forward-wedge two-component fit with a robust full-field fit of all three crack-local stress components to the leading Mode-I and Mode-II Williams fields.  It includes constant nonsingular stresses, Huber IRLS, and overlapping radial annuli.  The output reports annulus spread, normalized residual, conditioning, and a `projection_reliable` flag.

This is exact for isotropic LEFM.  For the anisotropic crystal calculation it remains an engineering projection, so reliability metrics must be retained in all plots.  A fully anisotropic interaction integral would require deeper changes to the FE kernel and auxiliary fields; it is not silently approximated here.

## Install

Copy the package contents into the `Arrhenius_FEM_CZM` project root.  New filenames are used throughout; v1 is not overwritten.

## Verify

```bash
chmod +x verify_mixed_mode_fem_czm_v2.sh run_mixed_mode_fem_czm_v2_campaign.sh
CONDA_ENV=arrhenius-fem-czm bash verify_mixed_mode_fem_czm_v2.sh
```

## Preflight

```bash
CONDA_ENV=arrhenius-fem-czm \
CLASSES="ceramic" \
TARGET_PSI="-30 0 30" \
MAX_CONTROL_ITERS=4 \
OUTROOT=runs/mixed_mode_fem_czm_v2_preflight \
bash run_mixed_mode_fem_czm_v2_campaign.sh
```

## Production screen

```bash
CONDA_ENV=arrhenius-fem-czm \
PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
CLASSES="ceramic DBTT" \
TARGET_PSI="-60 -45 -30 -15 0 15 30 45 60" \
T_K=500 \
PSI_TOL_DEG=2 \
MAX_CONTROL_ITERS=5 \
MAX_JOBS=1 \
OUTROOT=runs/mixed_mode_fem_czm_v2_event_controlled_500K \
bash run_mixed_mode_fem_czm_v2_campaign.sh
```

Each phase-angle condition stores every control trial, a control-history CSV, and one selected final summary.  DBTT cases that do not fracture remain explicitly right-censored; their phase angle is controlled at the terminal state rather than mislabeled as first passage.
