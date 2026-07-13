# Four-class FEM/CZM 500 K seeded R-curve campaign

This package runs the four active fully EXP-floor FEM/CZM classes

- `ceramic`
- `peak`
- `weakT`
- `DBTT`

at 500 K for five independent **solver-seed realizations**, each to a projected
crack extension of 1000 micrometers.

## Statistical interpretation

The constitutive model, barrier parameters, mesh, loading rate, and CZM settings
are identical in all five runs. The replicate index changes only the random-number
seed for the cleavage first-passage process.

For each renewal event the model draws

\[
H_n=-\ln U_n,\qquad U_n\sim \mathrm{Uniform}(0,1),
\]

so \(H_n\sim\mathrm{Exp}(1)\). The cumulative physical hazard action is
compared with this latent action. The sampled values are therefore consequences
of the Arrhenius hazard model; they are not five manually chosen threshold
parameters.

The same five seeds are reused across the four classes as common random numbers.
This gives paired stochastic histories and reduces finite-sample noise in
cross-class comparisons. The mesh seed is intentionally held fixed; varying it
would mix physical first-passage scatter with discretization sensitivity.

## Activation volume and Gumbel statistics

For a locally linear barrier

\[
\Delta G(\sigma)\simeq \Delta G_c-V_{\mathrm{eff}}(\sigma-\sigma_c),
\qquad
V_{\mathrm{eff}}=-\frac{\partial \Delta G}{\partial \sigma},
\]

the failure stress under monotonic loading has a minimum-type Gumbel form. Its
local stress scale is

\[
\beta_\sigma=\frac{k_B T}{V_{\mathrm{eff}}},
\qquad
\operatorname{SD}(\sigma_f)=\frac{\pi}{\sqrt 6}\beta_\sigma.
\]

Using the simplest local proportionality \(d\sigma_{\mathrm{eff}}/dK\simeq
\sigma_{\mathrm{eff}}/K\), the corresponding stress-intensity scale is

\[
\beta_K=\frac{k_B T}
{V_{\mathrm{eff}}\,d\sigma_{\mathrm{eff}}/dK}.
\]

The postprocessor reads the solver's `vstar_cleave_b3` diagnostic, computes these
local scales, and compares the predicted Gumbel standard deviation with the
empirical spread of the five extension-aligned curves. Because plastic shielding,
geometry, and state evolve, the full R-curve ensemble is generally a
history-dependent mixture rather than one exact Gumbel distribution. Five runs
are sufficient for a pilot mean and a diagnostic scatter check, not a definitive
distribution fit.

## Install

Copy the five files into the `Arrhenius_FEM_CZM` project root beside:

- `run_four_class_exp_floor_czm_500um_sweep.py`
- `four_class_exp_floor_exact_model_inputs.csv`
- `arrhenius_fracture/`

No installed package file is overwritten.

## Production run

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM

CONDA_ENV=arrhenius-fem-czm \
MAX_JOBS=1 \
OUTROOT=runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45 \
bash run_four_class_czm_500K_5rep_1000um.sh
```

Defaults:

- temperature: 500 K
- solver seeds: `1101 1102 1103 1104 1105`
- target extension: 1000 micrometers
- `dU=2e-7 m`, `dt=8.4 s`
- maximum steps: 50000
- branching disabled, `max_fronts=1`
- snapshots disabled to control disk usage

The campaign is restart-safe. Completed cases are skipped unless `FORCE=1`.

Enable snapshots every 100 micrometers with:

```bash
SAVE_SNAPSHOTS=11 SNAPSHOT_BY_EXT_UM=100 \
bash run_four_class_czm_500K_5rep_1000um.sh
```

## Plot only

```bash
PLOT_ONLY=1 \
OUTROOT=runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45 \
bash run_four_class_czm_500K_5rep_1000um.sh
```

## Main outputs

```text
four_class_500K_individual_and_mean_R_curves.png
four_class_500K_mean_R_curve_comparison.png
four_class_500K_activation_volume_vs_extension.png
four_class_500K_replicate_summary.csv
four_class_500K_activation_volume_gumbel_summary.csv
four_class_500K_R_curve_ensemble_and_gumbel_statistics.csv
```

Each class also receives:

```text
R_curve_replicates_mean_and_gumbel_scatter.png
R_curve_ensemble_and_gumbel_statistics.csv
```

The mean is calculated after interpolation onto a common 5 micrometer extension
grid, so runs with more recorded events do not receive greater statistical
weight.
