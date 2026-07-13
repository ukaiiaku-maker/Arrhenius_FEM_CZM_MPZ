# V1 fully EXP-floor four-class tuning

This calibration searches for new 1-D parameterizations in which **both crack-tip emission and crack opening use EXP-floor barriers**.

## Calibration target

The fit target is the observable prior 1-D response from `forward_input_parameter_KT_overlay`:

- ceramic-like `Kc(T)`
- peak-like `Kc(T)`
- weak-temperature-dependence `Kc(T)`
- DBTT-like `Kc(T)`

The old linearized crack-opening barrier parameters are **not** treated as targets. The prior 2-D toughness points are carried only as validation overlays in the final figure.

## Search architecture

The default production search has three stages:

1. **Broad stage**: 256 Sobol cleavage/state contexts crossed with the 96 supplied EXP-floor emission surfaces = 24,576 candidates.
2. **Local stage**: perturbations around the best broad candidates from each target class.
3. **Final stage**: fine `dK` evaluation of the finalists over the dense 5 K target grid.

The searched cleavage/state parameters are:

- `cleave_G00_eV`
- `cleave_gT_eV_per_K`
- `cleave_sigc0_GPa`
- `cleave_sT_MPa_per_K`
- `cleave_exp_a`
- `cleave_exp_n`
- `cleave_floor_frac`
- `cleave_S_hs_kB`
- `chi_shield`
- `N_sat`

The emission channel uses the 96-surface EXP-floor design from the later fatigue/strength work.

## Production run

```bash
OUT=runs/v1_exp_floor_four_class_tuning \
N_CONTEXTS=256 \
BROAD_DK=0.25 \
LOCAL_DK=0.05 \
FINAL_DK=0.02 \
RESUME=1 \
bash run_v1_exp_floor_four_class_tuning.sh
```

The default loading rate is `KDOT=0.005 MPa sqrt(m)/s`, matching the prior forward-model target calculation.

## Smoke test

```bash
SMOKE=1 \
OUT=runs/v1_exp_floor_four_class_tuning_smoke \
bash run_v1_exp_floor_four_class_tuning.sh
```

## Restart behavior

`RESUME=1` reuses completed broad, local, and final stage arrays/tables in the same `OUT` directory. Use a new output directory when changing search-grid settings or parameter ranges.

## Main outputs

- `recommended_exp_floor_four_class.csv`: best parameterization for each class
- `recommended_curves_dense.csv`: tuned dense-grid response curves
- `exp_floor_four_class_tuning.png`: target, tuned V1, and prior 2-D validation overlay
- `broad_candidate_scores.csv`
- `local_candidate_scores.csv`
- `finalist_scores.csv`
- `run_config.json`

The recommended parameters are intended as seeds for the expensive adaptive-CZM FEM sweep, not as unique inverse solutions.
