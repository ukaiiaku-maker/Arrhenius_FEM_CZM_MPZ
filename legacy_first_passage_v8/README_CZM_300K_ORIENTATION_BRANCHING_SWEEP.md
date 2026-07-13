# FEM/CZM 300 K R-curve orientation and weak-branching sweep

This package adds a run script and plotter for a one-class 300 K FEM/CZM R-curve campaign.
It does not modify the solver. It templates a completed FEM/CZM command from the existing 500 K replicate campaign or rate sweep.

Default campaign:

- class: `ceramic`
- temperature: 300 K
- orientations: 0, 15, 30, 45 degrees, no branching (`--max-fronts 1`)
- branching comparison: 30 degrees, weak branching (`--max-fronts 3`)
- target crack extension: 1000 µm
- seed: 1201
- output: `runs/czm_Rcurve_300K_orientation_branching_ceramic_v1`

Run from the `Arrhenius_FEM_CZM` project root:

```bash
unzip -o CZM_300K_Orientation_Branching_Rcurve_Sweep.zip

CONDA_ENV=arrhenius-fem-czm \
CLASS=ceramic \
T_K=300 \
THETAS="0 15 30 45" \
BRANCH_THETA=30 \
SEED=1201 \
TARGET_EXT_UM=1000 \
STEPS=80000 \
OUTROOT=runs/czm_Rcurve_300K_orientation_branching_ceramic_v1 \
bash run_czm_300K_orientation_branching_sweep.sh
```

For the weakT alternative, use:

```bash
CLASS=weakT OUTROOT=runs/czm_Rcurve_300K_orientation_branching_weakT_v1 bash run_czm_300K_orientation_branching_sweep.sh
```

Plot after completion:

```bash
PYTHON_BIN=/opt/homebrew/Caskroom/miniconda/base/envs/arrhenius-fem-czm/bin/python \
ROOT=runs/czm_Rcurve_300K_orientation_branching_ceramic_v1 \
bash run_plot_czm_300K_orientation_branching_Rcurves.sh
```

Outputs:

- `plots/orientation_sweep_no_branch_Rcurves.png`
- `plots/branching_comparison_Rcurves.png`
- `plots/orientation_branching_Rcurve_summary.csv`
- `plots/orientation_branching_binned_Rcurves_long.csv`
