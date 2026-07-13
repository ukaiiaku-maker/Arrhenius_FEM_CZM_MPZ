# Mixed-mode FEM/CZM v3: J-consistent isotropic first passage

## Why v2 results are not production data

The v2 campaign used an isotropic Williams stress fit while running rotated cubic anisotropic elasticity. In the submitted campaign, most fits had normalized residuals near 0.9 and large annulus-to-annulus phase spread. The same unreliable fit was used both for the controller and for the cleavage/emission drives. This caused three separate problems:

1. a numerically poor fit changed the physical hazard;
2. an earlier phase-accurate iteration could be rejected in favor of a less accurate row with a passing Boolean flag;
3. event, right-censored, phase-error, and projection-error states were collapsed into `not_converged`.

The v2 anisotropic mixed-mode results should therefore be treated as diagnostics only.

## v3 formulation

v3 is a reliable isotropic geometry screen. It uses the domain-integral amplitude `KJ` and the elastically calibrated target phase angle:

```text
KI  = KJ cos(psi_target)
KII = KJ sin(psi_target)
```

Therefore `KI^2 + KII^2 = KJ^2` exactly. These values drive maximum-hoop cleavage and shear-assisted emission. The multi-annulus Williams fit is retained only as an independent diagnostic and cannot change the hazard.

The phase angle is calibrated once in an isotropic linear-elastic solve on the same mesh. Before the first crack event the geometry and elastic operator are unchanged, so the phase does not require event-state feedback iterations.

## Scope

This version deliberately refuses anisotropic quantitative control. A rotated cubic crystal requires a proper anisotropic interaction integral or Stroh auxiliary fields. It is better to run a correct isotropic mode-mixity screen than to assign false precision to an isotropic fit of anisotropic stresses.

## Installation

Copy the package contents into the `Arrhenius_FEM_CZM` project root. All filenames are version-specific.

## Verify

```bash
chmod +x verify_mixed_mode_fem_czm_v3.sh run_mixed_mode_fem_czm_v3_campaign.sh
CONDA_ENV=arrhenius-fem-czm bash verify_mixed_mode_fem_czm_v3.sh
```

Expected ending:

```text
Ran 5 tests
OK
MIXED_MODE_V3 verification OK
```

## Preflight

```bash
CONDA_ENV=arrhenius-fem-czm \
CLASSES="ceramic" \
TARGET_PSI="-30 0 30" \
OUTROOT=runs/mixed_mode_fem_czm_v3_preflight \
bash run_mixed_mode_fem_czm_v3_campaign.sh
```

## Full short screen

```bash
CONDA_ENV=arrhenius-fem-czm \
PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
CLASSES="ceramic DBTT" \
TARGET_PSI="-60 -45 -30 -15 0 15 30 45 60" \
T_K=500 \
MAX_JOBS=1 \
OUTROOT=runs/mixed_mode_fem_czm_v3_J_consistent_500K \
bash run_mixed_mode_fem_czm_v3_campaign.sh
```

There is one deterministic run per class and phase angle. Status is `event`, `right_censored`, or `failed`; censoring is not called nonconvergence.
