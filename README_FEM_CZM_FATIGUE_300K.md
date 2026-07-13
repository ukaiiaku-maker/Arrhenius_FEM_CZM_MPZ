# FEM/CZM six-case fatigue campaign v4

This package retains the v8 cyclic hazard controller and uses the adaptive-CZM crack backend.

## Important smoke-test interpretation

The smoke test is only a wiring and numerical-preflight test. With `BLOCKS=80`, `Kmax="4 6 8"`, and a 50 um target, it may produce no accepted crack advance. Such a run cannot establish a fatigue threshold or Paris slope. No-growth cases are plotted as one-increment upper bounds.

## Production-style plots

The analysis now reproduces the prior atlas conventions:

- `analysis/atlas_2d_da_dN_vs_DeltaK_with_upper_bounds.png`
- `analysis/atlas_2d_local_da_dN_vs_local_DeltaK.png`
- `analysis/atlas_2d_paris_points.csv`
- `analysis/atlas_2d_local_paris_points.csv`

The integrated rate is total accepted extension divided by total cycles. If no advance occurs, the plotted upper bound is one physical advance increment divided by total cycles. Local points use `da_block_m/fatigue_cycles` versus `(1-R)KJ` for accepted growth blocks.

## Replot existing outputs only

```bash
PYTHON_BIN=/opt/homebrew/Caskroom/miniconda/base/envs/arrhenius-fem-czm/bin/python \
$PYTHON_BIN analyze_fem_czm_fatigue_outputs.py \
  --root runs/fem_czm_six_fatigue_300K_smoke \
  --R 0.1 \
  --cycles-max 1e10 \
  --target-crack-extension-um 50
```

## Focused shielded-case validation

Use a dense grid around the prior threshold region. This is longer than a smoke test but much smaller than the full six-case campaign.

```bash
CONDA_ENV=arrhenius-fem-czm \
MODE=pilot \
CASE_FILTER=plastic_shielded_case64_M1 \
KLIST_OVERRIDE="5.0 5.5 6.0 6.25 6.5 6.75 7.0 7.25 7.5 8.0" \
OUTROOT=runs/fem_czm_plastic_shielded_threshold_300K_pilot \
BLOCKS=5000 \
CYCLES_MAX=2e14 \
TARGET_EXT_UM=100 \
MAKE_2D_PLOTS=0 \
bash run_fem_czm_six_fatigue_300K.sh
```

For manuscript-quality local Paris points, increase `TARGET_EXT_UM` to 250 um after the focused validation behaves correctly.


## v5 analysis corrections

- Integrated atlas x-values always use the initial range `(1-R) K_initial`; geometry-evolved final K is reserved for local block plots.
- No-growth points that stop at the numerical block limit are shown as open downward triangles and labeled block-limited bounds.
- True cycle-horizon upper bounds are shown as filled downward triangles.
- `fem_czm_rate_defined_thresholds.csv` reports log-interpolated rate-defined thresholds using measured integrated rates.


## Local Paris-point color convention (v6)

The local `da/dN` versus local `Delta K_J` plot now uses a continuous `viridis`
color scale for accumulated crack extension, `Delta a` in micrometres. Material
class is represented by marker shape. The established output filenames are
retained, and an additional explicit `*_colored_by_extension` alias is written.
No simulations need to be rerun; rerunning the analysis is sufficient.
