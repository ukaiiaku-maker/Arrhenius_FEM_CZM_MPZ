# Saved R-curve analysis

This package analyzes the saved FEM/CZM seed outputs without rerunning simulations.

It reads, in order of preference:

1. `R_curve_event_sampled.csv`
2. `steps_*K.csv`

and then performs a bin-median R-curve analysis.

## Run

From the `Arrhenius_FEM_CZM` project root:

```bash
unzip -o Saved_RCurve_Analysis.zip

PYTHON_BIN=/opt/homebrew/Caskroom/miniconda/base/envs/arrhenius-fem-czm/bin/python \
ROOT=runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45 \
OUT=runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/Rcurve_analysis \
bash run_analyze_saved_Rcurves.sh
```

Optional controls:

```bash
BIN_UM=25
LATE_WINDOW="700 1000"
MAX_EXT_UM=1000
```

## Main outputs

- `*_binned_seed_Rcurves_fit.png`
- `seed_Rcurve_metrics_and_fits.csv`
- `class_Rcurve_metric_summary_complete_only.csv`
- `class_mean_Rcurve_fits.csv`
- `seed_binned_Rcurves_long.csv`

The analysis is R-curve-like and simulation-facing, not an ASTM-valid fracture test analysis.
