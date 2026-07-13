# Mixed-mode FEM/CZM v4 — anisotropic traction-controlled first passage

## Purpose

V4 activates the existing cubic anisotropic FEM and crystallographic cleavage competition.  It does **not** obtain kinetics by partitioning an isotropic equivalent K into `KI=KJ cos(psi)` and `KII=KJ sin(psi)`.

The controlled mode coordinate is the finite-radius process-zone traction phase

`psi_sigma = atan2(tau_tn, sigma_nn)`

on the original crack plane.  The phase is calibrated from two anisotropic elastic basis solves.  During the fracture run, the actual anisotropic FEM stress tensor is sampled at a fixed physical radius.  Cleavage uses resolved opening on the selected crystal plane divided by `sqrt(gamma_rel)`; emission uses the maximum resolved shear on the BCC slip traces.  The domain J integral remains an energy-release audit.

This is the appropriate production-ready intermediate step for first passage.  It is more direct than an isotropic Williams projection and does not pretend that isotropic KI/KII are valid for a rotated cubic crystal.  A full Stroh/anisotropic interaction integral remains desirable for geometry-independent anisotropic SIF reporting, especially once the crack is extended far enough that its direction and geometry evolve.

## Default material

Tungsten single-crystal constants: C11=523 GPa, C12=203 GPa, C44=160 GPa (Zener ratio approximately 1).  Consequently elastic anisotropy alone is expected to be weak for W; the stronger directional effect comes from discrete cleavage planes and the finite cleavage-energy anisotropy.  The default crystal angle is 0 degrees and `cleave_gamma_aniso=0.3`.

## Install

Copy this directory's contents into the project root.  All files use v4 names and do not replace v3.5.

## Verify

```bash
CONDA_ENV=arrhenius-fem-czm bash verify_mixed_mode_fem_czm_v4.sh
```

## Three-angle preflight

```bash
rm -rf runs/mixed_mode_fem_czm_v4_anisotropic_preflight_cal
rm -rf runs/mixed_mode_fem_czm_v4_anisotropic_preflight

CONDA_ENV=arrhenius-fem-czm \
CLASSES="ceramic" \
TARGET_PSI="-30 0 30" \
CRYSTAL_THETA_DEG=0 \
RECALIBRATE=1 \
CALROOT=runs/mixed_mode_fem_czm_v4_anisotropic_preflight_cal \
OUTROOT=runs/mixed_mode_fem_czm_v4_anisotropic_preflight \
bash run_mixed_mode_fem_czm_v4_campaign.sh
```

Inspect `traction_phase_error_first_deg`, `traction_probe_reliable`, and the +/- symmetry before the full campaign.

## Full campaign

```bash
CONDA_ENV=arrhenius-fem-czm \
PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
CLASSES="ceramic DBTT" \
TARGET_PSI="-60 -45 -30 -15 0 15 30 45 60" \
T_K=500 \
CRYSTAL_THETA_DEG=0 \
MAX_JOBS=1 \
RECALIBRATE=1 \
CALROOT=runs/mixed_mode_fem_czm_v4_anisotropic_calibration_500K \
OUTROOT=runs/mixed_mode_fem_czm_v4_anisotropic_500K \
bash run_mixed_mode_fem_czm_v4_campaign.sh
```

No seed sweep is used; H=1 is deterministic.
