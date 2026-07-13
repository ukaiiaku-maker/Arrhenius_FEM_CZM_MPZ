# Mixed-mode FEM/CZM first-passage campaign v1

## Purpose

This version-specific package adds mixed-mode loading to the active FEM/CZM sharp-front solver without replacing `sharp_front.py`.

It provides:

1. Combined opening/sliding displacement loading on the fixed cracked-plate geometry.
2. Elastic calibration from boundary loading angle `alpha` to the achieved FEM stress-intensity phase angle
   `psi = atan2(KII, KI)`.
3. Near-tip extraction of `KI` and `KII` from the local FEM stress field.
4. An opening-sensitive cleavage drive based on the maximum-hoop component.
5. A shear-assisted emission drive.
6. Seeded exponential cumulative-hazard thresholds for stochastic first passage.
7. Short first-passage campaign and plotting scripts.

The existing EXP-floor barriers, plasticity, adaptive CZM insertion, crystal competition, and local/cluster J machinery remain in the active solver.

## Version-specific files

Copy the files without renaming them:

```text
Arrhenius_FEM_CZM/
├── arrhenius_fracture/
│   └── mixed_mode_first_passage_v1.py
├── tests/
│   └── test_mixed_mode_first_passage_v1.py
├── calibrate_mixed_mode_loading_v1.py
├── run_mixed_mode_fem_czm_v1_campaign.py
├── plot_mixed_mode_fem_czm_v1_results.py
├── run_mixed_mode_fem_czm_v1_campaign.sh
└── verify_mixed_mode_fem_czm_v1.py
```

No existing Python module or runner is overwritten.

## Model definition

Boundary loading is prescribed through total relative displacement amplitude `A` and loading angle `alpha`:

```text
Un = A cos(alpha)
Us = A sin(alpha)
```

The boundary solver first obtains the opening equilibrium with free lateral contraction, then superposes the sliding displacement while retaining the opening solution's lateral shape.

At each J-integral evaluation, local stresses in a forward crack-tip wedge are fitted to

```text
sigma_22 = KI / sqrt(2 pi r) + T22
sigma_12 = KII / sqrt(2 pi r) + T12
```

The cleavage drive is the maximum positive hoop-stress coefficient

```text
Kopen(theta) = KI cos^3(theta/2) - 3 KII sin(theta/2) cos^2(theta/2).
```

The emission drive is

```text
Kemit = sqrt(Kopen^2 + w_shear KII^2),
```

where `w_shear` defaults to one and is explicitly recorded.

For solver seed `s`, the first cumulative-hazard threshold is

```text
H* = -ln(U_s),  U_s ~ Uniform(0,1).
```

The barrier parameters are unchanged. The threshold is a realization of the same Poisson first-passage model, not an adjusted material parameter.

## Verify installation

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM
conda activate arrhenius-fem-czm

rm -rf arrhenius_fracture/__pycache__ tests/__pycache__
python verify_mixed_mode_fem_czm_v1.py
```

Expected ending:

```text
Ran 4 tests
OK
MIXED_MODE_V1 verification OK
```

## Elastic calibration only

Start with the isotropic calibration:

```bash
python calibrate_mixed_mode_loading_v1.py \
  --out runs/mixed_mode_fem_czm_v1_elastic_calibration \
  --target-psi-deg "-60 -45 -30 -15 0 15 30 45 60"
```

Review:

```text
mixed_mode_loading_calibration.csv
mixed_mode_loading_calibration.png
```

The achieved `psi` should be close to each target. Positive boundary sliding may correspond to negative `KII` under the solver sign convention; the calibration absorbs this sign convention.

## Initial short campaign

The default campaign runs `ceramic` and `DBTT` at 500 K, nine phase angles, and three solver seeds:

```bash
CONDA_ENV=arrhenius-fem-czm \
PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
CLASSES="ceramic DBTT" \
SEEDS="1101 1102 1103" \
MAX_JOBS=1 \
OUTROOT=runs/mixed_mode_fem_czm_v1_first_passage_500K \
bash run_mixed_mode_fem_czm_v1_campaign.sh
```

For a smaller preflight:

```bash
python run_mixed_mode_fem_czm_v1_campaign.py \
  --parameter-table four_class_exp_floor_exact_model_inputs.csv \
  --calibration-csv runs/mixed_mode_fem_czm_v1_elastic_calibration/mixed_mode_loading_calibration.csv \
  --classes "ceramic" \
  --target-psi-deg "-30 0 30" \
  --seeds "1101" \
  --T-K 500 \
  --outroot runs/mixed_mode_fem_czm_v1_preflight \
  --max-jobs 1
```

Each case stops after its first accepted crack event and writes:

```text
mixed_mode_first_passage_summary.json
mixed_mode_first_passage_summary.csv
mixed_mode_projection_calls.csv
steps_0500K.csv
crack_path_500K.csv
run.log
command.txt
```

## Main plots

The plotting script generates:

```text
KI_KII_first_passage_envelope.png
Kopen_vs_mode_phase.png
kink_angle_vs_mode_phase.png
mixed_mode_first_passage_all_cases.csv
mixed_mode_first_passage_grouped.csv
```

## Required validation before production interpretation

1. Isotropic `+psi/-psi` symmetry.
2. 16/32/64 or equivalent near-tip angular/radial projection convergence.
3. Root mesh refinement of `KI`, `KII`, and first-passage load.
4. Loading-angle calibration repeated for the anisotropic crystal orientation.
5. Comparison of maximum-hoop kink direction with the inserted first CZM segment.
6. Pure or near-pure Mode II should not be used until crack-face contact is added.

## Current limitations

- `KI/KII` extraction uses the isotropic leading Williams field. It is exact for the isotropic verification and an approximate diagnostic for cubic anisotropy.
- v1 uses a scalar opening drive and a shear-assisted emission drive rather than a full tensorial activation-volume law.
- The fixed top/bottom displacement geometry is a controlled computational mixed-mode test, not yet an experimental specimen geometry.
- Branching is intentionally disabled; the campaign is first passage only.
- The adaptive CZM backend must already exist in the active project.
