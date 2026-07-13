# Three-rate FEM/CZM toughness comparison with V1 analytical predictions

This package plots the four-class FEM/CZM rate sweep at 1x, 10x, and 100x and overlays V1 analytical predictions at the corresponding rates.

It produces two main 2x2 panel figures:

- `CZM_three_rate_first_passage_vs_analytic.png`
- `CZM_three_rate_median_Rcurve_vs_analytic.png`

The first uses the first-passage fracture toughness. The second uses the median of all event-sampled `K_J` values in each case's `R_curve_event_sampled.csv`.

Incomplete runs are retained and plotted as open symbols.

## Run from `Arrhenius_FEM_CZM`

```bash
unzip -o CZM_Three_Rate_Toughness_With_Analytic.zip

ROOT=runs/four_class_exp_floor_CZM_rates_no_branch_500um_theta45 \
OUT=runs/CZM_three_rate_temperature_comparison \
RATES="1 10 100" \
bash run_plot_czm_three_rate_toughness_with_analytic.sh
```

## Analytical curves

The script first tries to compute rate-specific V1 analytical curves by importing:

```text
run_v1_exp_floor_four_class_tuning.py
```

and using:

```text
four_class_exp_floor_exact_model_inputs.csv
```

with:

```text
Kdot = BASE_KDOT * rate_factor
```

where the default is `BASE_KDOT=0.005`.

If the V1 script is not present or cannot be imported, the script falls back to `four_class_analytical_prediction_final.csv`. If that CSV does not contain `rate_factor`, the same analytical curve is reused for all rates and a warning is printed.


## v2 notes

This version fixes two issues:

1. The script now scans per-case folders even when a top-level `four_class_temperature_summary.csv` exists. This prevents a stale summary file from hiding valid `rate_100x` cases.
2. The wrapper runs with `CONDA_ENV=arrhenius-fem-czm` by default, so the V1 analytical computation uses the project environment rather than the base environment. This should avoid the SciPy `_spropack` import error from base Python 3.13.

It also writes `rate_case_availability.csv` to show exactly which class/temperature/rate cases were found and whether each completed 500 µm.


## v3 fix

The script now supplies default V1 metadata columns such as `exp_Tref_K=300` if they are absent from the compact exact-parameter CSV. This fixes the previous analytical-computation error:

```text
WARNING: V1 analytical computation failed: 'exp_Tref_K'
```
