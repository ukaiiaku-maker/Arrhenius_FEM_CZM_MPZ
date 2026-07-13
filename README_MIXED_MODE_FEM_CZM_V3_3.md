# Mixed-mode FEM/CZM v3.3 — circular phase-spread calibration

This release retains the v3.2 signed two-basis calibration and fixes the
calibration failure observed for the nominal Mode-I target.

## Corrected defect

For the raw opening basis, the fitted annular phase angles lie near the
`-180/+180 degree` branch cut because the raw Mode-I sign is opposite to the
physical opening convention.  V3.2 computed the annular spread with a linear
median absolute deviation, so nearby angles such as `+179` and `-179` were
incorrectly treated as nearly 360 degrees apart.  This produced the false
`phase_spread = 266.7 deg` rejection even though the achieved target was
`psi = 0.004 deg`.

V3.3 uses a wrapped circular median absolute deviation:

- angular differences are reduced to `[-180, 180)`;
- the robust center is selected from the observed annular phases;
- a true mixture of opposite modes, such as 0 and 180 degrees, remains rejected.

The output records
`mode_projection_psi_spread_method = wrapped_circular_MAD`.

## Verification

```bash
CONDA_ENV=arrhenius-fem-czm bash verify_mixed_mode_fem_czm_v3_3.sh
```

Expected ending:

```text
Ran 13 tests
OK
MIXED_MODE_V3_3 verification OK
```

## Calibration preflight

```bash
rm -rf runs/mixed_mode_fem_czm_v3_3_circular_phase_preflight

conda run -n arrhenius-fem-czm \
  python calibrate_mixed_mode_loading_v3_3.py \
  --out runs/mixed_mode_fem_czm_v3_3_circular_phase_preflight \
  --target-psi-deg="-30 0 30"
```

The zero-degree row should no longer fail merely because annular phase values
straddle the branch cut.  `amplitude_fit_ok=False` remains acceptable because
the authoritative amplitude comes from the domain J integral.

## Three-condition fracture preflight

```bash
CONDA_ENV=arrhenius-fem-czm \
CLASSES="ceramic" \
TARGET_PSI="-30 0 30" \
CALROOT=runs/mixed_mode_fem_czm_v3_3_circular_phase_preflight \
OUTROOT=runs/mixed_mode_fem_czm_v3_3_preflight \
bash run_mixed_mode_fem_czm_v3_3_campaign.sh
```

## Full campaign

```bash
CONDA_ENV=arrhenius-fem-czm \
PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
CLASSES="ceramic DBTT" \
TARGET_PSI="-60 -45 -30 -15 0 15 30 45 60" \
T_K=500 \
MAX_JOBS=1 \
CALROOT=runs/mixed_mode_fem_czm_v3_3_circular_phase_calibration \
OUTROOT=runs/mixed_mode_fem_czm_v3_3_circular_phase_500K \
bash run_mixed_mode_fem_czm_v3_3_campaign.sh
```
