# Mixed-mode FEM/CZM v6: exact production-backend control

## Purpose

V6 corrects the remaining mode-control defect identified in the v5 anisotropic
preflight. V5 calibrated loading on a separately assembled static damaged notch,
while physical first passage used the adaptive-CZM backend. At the first physical
step, the requested `-30 deg` case therefore measured only about `+1.65 deg`, and
its phase evolved to about `+14.4 deg` before fracture.

V6 performs every elastic calibration probe through the same production path used
by the fracture calculation:

- the same `sharp_front.run_2d` driver;
- `crack_backend=adaptive_czm`;
- the same mesh and refinement;
- the same CZM initialization and penalties;
- cubic anisotropic elasticity and crystal orientation;
- the same process-zone traction probe and J settings;
- the same mixed boundary-condition implementation.

No standalone `d[bnd.notch_nodes]=1` calibration model is used.

## Loading parameterization

V5 controlled the boundary condition through an angle capped at `89.9 deg`. The
required opening displacement can be much smaller than the sliding displacement,
so that cap can prevent a reachable mode mixity from being represented.

V6 passes direct normalized coefficients:

```text
U_open  = U_total * q_open
U_shear = U_total * q_shear
q_open^2 + q_shear^2 = 1
```

For event-state control it uses the unconstrained coordinate

```text
z = asinh(U_shear/U_open)
```

Large shear/opening ratios can therefore be represented without an artificial
near-90-degree saturation.

## Calibration sequence

For each campaign calibration, V6 runs:

1. one exact-backend opening basis probe;
2. one exact-backend sliding basis probe;
3. one exact-backend Mode-I reference probe;
4. one exact-backend verification probe for each requested target phase.

Every probe is a one-step adaptive-CZM production run with mechanically passive
high Arrhenius barriers. The generated CSV must report

```text
first_production_step_verified = True
phase_converged = True
```

before a long physical calculation is allowed to start.

For a three-angle preflight this means six one-step probe runs before the six
physical class/angle cases. This is expected and can make the start appear slower
than earlier versions.

## Event-state phase controller

The local phase can evolve with load because of plasticity, cohesive compliance,
and state evolution. V6 therefore measures phase at first passage (or at the
right-censored endpoint) and updates `z` empirically.

The controller uses, in order:

1. a sign-changing event-state bracket and safeguarded secant update;
2. an empirical secant from prior completed iterations;
3. the exact-backend one-step response matrix only for the first local update.

The calibration matrix no longer dictates all later corrections. If the target is
not reachable in the geometry, the result remains explicitly labeled
`event_phase_mismatch` or `right_censored_phase_mismatch`.

## Absolute hazard scale

V6 preserves the corrected v5 kinetics:

- domain-integral `KJ` supplies the calibrated sharp-tip magnitude;
- anisotropic FEM stress supplies dimensionless cleavage/slip directional factors;
- raw finite-radius traction does not replace the calibrated sharp-tip stress;
- barrier fingerprints and cumulative clocks remain audited by class.

## Installation

Copy the contents of this folder into the root of `Arrhenius_FEM_CZM`. The
versioned files do not overwrite v5.

## Verification

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM
chmod +x verify_mixed_mode_fem_czm_v6.sh
chmod +x run_mixed_mode_fem_czm_v6_campaign.sh

CONDA_ENV=arrhenius-fem-czm \
bash verify_mixed_mode_fem_czm_v6.sh
```

Expected:

```text
Ran 12 tests
OK
MIXED_MODE_V6 verification OK
```

## Required three-angle preflight

```bash
rm -rf runs/mixed_mode_fem_czm_v6_preflight_cal
rm -rf runs/mixed_mode_fem_czm_v6_preflight

CONDA_ENV=arrhenius-fem-czm \
PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
CLASSES="ceramic DBTT" \
TARGET_PSI="-30 0 30" \
T_K=500 \
CRYSTAL_THETA_DEG=45 \
TRACTION_PROBE_RADIUS_M=1e-5 \
EVENT_PSI_TOL_DEG=2 \
MAX_CONTROL_ITERS=6 \
MAX_JOBS=1 \
RECALIBRATE=1 \
CALROOT=runs/mixed_mode_fem_czm_v6_preflight_cal \
OUTROOT=runs/mixed_mode_fem_czm_v6_preflight \
bash run_mixed_mode_fem_czm_v6_campaign.sh
```

Inspect:

```text
runs/mixed_mode_fem_czm_v6_preflight_cal/mixed_mode_loading_calibration_v6.csv
runs/mixed_mode_fem_czm_v6_preflight/<class>/<phase>/mixed_mode_control_history_v6.csv
runs/mixed_mode_fem_czm_v6_preflight/campaign_status_v6.csv
runs/mixed_mode_fem_czm_v6_preflight/plots_v6/
```

Acceptance for all three targets:

```text
calibration_first_production_step_verified = True
event_phase_control_converged = True
status = event  OR  right_censored
```

A status containing `phase_mismatch` is diagnostic and should not enter the
controlled mixed-mode envelope.

## Full first-passage campaign

Run only after the preflight passes:

```bash
rm -rf runs/mixed_mode_fem_czm_v6_full_cal
rm -rf runs/mixed_mode_fem_czm_v6_anisotropic_500K

CONDA_ENV=arrhenius-fem-czm \
PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
CLASSES="ceramic DBTT" \
TARGET_PSI="-60 -45 -30 -15 0 15 30 45 60" \
T_K=500 \
CRYSTAL_THETA_DEG=45 \
TRACTION_PROBE_RADIUS_M=1e-5 \
EVENT_PSI_TOL_DEG=2 \
MAX_CONTROL_ITERS=6 \
MAX_JOBS=1 \
RECALIBRATE=1 \
CALROOT=runs/mixed_mode_fem_czm_v6_full_cal \
OUTROOT=runs/mixed_mode_fem_czm_v6_anisotropic_500K \
bash run_mixed_mode_fem_czm_v6_campaign.sh
```

The longer crack-extension and median-threshold campaign should be built only
after this controlled first-passage map is validated.
