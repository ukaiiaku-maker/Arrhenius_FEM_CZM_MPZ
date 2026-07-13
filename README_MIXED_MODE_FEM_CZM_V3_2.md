# Mixed-mode FEM/CZM v3.2 — signed-basis phase calibration

This version fixes the v3.1 elastic calibration failure seen when the Williams
projection used the opposite Mode-I sign convention from the imposed opening
displacement.  It also raises the admissible loading-angle limit from 85 to
89.9 degrees because this geometry has an opening/sliding basis-amplitude ratio
of approximately 31:1.

## Physical convention

The pure-opening elastic basis defines positive Mode I.  The pure positive-sliding
basis defines positive Mode II.  The raw Williams outputs are automatically
multiplied by two stored sign factors before phase ratios are calculated.

Both the raw and normalized basis matrices are written to JSON and CSV.  This is
a sign-convention normalization only; it does not change the FEM solution.

## Verification

```bash
CONDA_ENV=arrhenius-fem-czm bash verify_mixed_mode_fem_czm_v3_2.sh
```

## Calibration preflight

```bash
rm -rf runs/mixed_mode_fem_czm_v3_2_signed_basis_preflight
conda run -n arrhenius-fem-czm python calibrate_mixed_mode_loading_v3_2.py \
  --out runs/mixed_mode_fem_czm_v3_2_signed_basis_preflight \
  --target-psi-deg="-30 0 30"
```

For the basis matrix reported from v3.1, approximate initial boundary angles are
-87.21, -0.20, and +86.37 degrees.  Large displacement angles are expected because
unit sliding displacement produces much less fitted mode amplitude than unit opening.

## Campaign

```bash
CONDA_ENV=arrhenius-fem-czm \
CLASSES="ceramic DBTT" \
TARGET_PSI="-60 -45 -30 -15 0 15 30 45 60" \
T_K=500 \
MAX_JOBS=1 \
CALROOT=runs/mixed_mode_fem_czm_v3_2_signed_basis_calibration \
OUTROOT=runs/mixed_mode_fem_czm_v3_2_signed_basis_500K \
bash run_mixed_mode_fem_czm_v3_2_campaign.sh
```
