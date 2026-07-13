# Mixed-mode FEM/CZM v7: exact-backend full-circle control

## Purpose

V7 fixes the two control defects exposed by the v6 production-backend calibration:

1. `traction_probe_reliable` combined two different questions: whether the local normal/shear traction phase was measurable, and whether a crystallographic cleavage/slip candidate was valid. Consequently, the 0° and +30° phase probes were rejected even though their phase errors were only 0.110° and 0.072°.
2. The v6 coordinate `z = asinh(U_shear/U_open)` forced the opening coefficient to be positive. The reported +30° basis solution actually had a small negative opening coefficient; v6 silently reversed that vector when forming `z`, so the stored coordinate and the applied boundary coefficients were inconsistent.

The -30° result showed a third issue: the linear basis inversion predicted -30°, but the exact adaptive-CZM probe produced only -4.33°. That state is cancellation-sensitive, so the basis matrix is suitable only as an initial guess.

## V7 formulation

Boundary loading is represented on the complete unit circle:

```
U_open  = U_total cos(alpha)
U_shear = U_total sin(alpha)
```

`alpha` is retained as an unwrapped full-circle coordinate. Negative remote opening is not silently changed; it is recorded through:

```
loading_open_is_tensile
loading_open_coeff
loading_shear_coeff
loading_alpha_deg
```

For every requested traction phase, the calibrator:

1. runs exact one-step adaptive-CZM pure-opening and pure-sliding probes;
2. uses the measured response matrix only to initialize `alpha`;
3. runs the exact production backend at that initial state;
4. performs local secant/bracket/root searches directly on the measured phase;
5. accepts only a finite phase probe within `CAL_PSI_TOL_DEG`.

The calibration no longer requires a valid crystallographic candidate. That is a separate physical-run requirement. The files now distinguish:

```
traction_phase_probe_reliable
directional_metrics_reliable
traction_probe_reliable  # both conditions
```

The physical campaign performs a second event-state full-circle control loop and accepts a production point only when both phase control and directional metrics are reliable.

## Installation

Copy the contents of this folder into the root of `Arrhenius_FEM_CZM`. V6 files are not overwritten.

## Verification

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM
chmod +x verify_mixed_mode_fem_czm_v7.sh
chmod +x run_mixed_mode_fem_czm_v7_campaign.sh

CONDA_ENV=arrhenius-fem-czm \
bash verify_mixed_mode_fem_czm_v7.sh
```

Expected ending:

```
Ran 13 tests
OK
MIXED_MODE_V7 verification OK
```

The tests include the exact v6 response matrix, preservation of a negative opening coefficient, full-circle round trips, separate phase reliability, and a synthetic exact-backend root reproducing the reported negative-branch failure.

## Required three-angle preflight

Use new folders:

```bash
rm -rf runs/mixed_mode_fem_czm_v7_preflight_cal
rm -rf runs/mixed_mode_fem_czm_v7_preflight

CONDA_ENV=arrhenius-fem-czm \
PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
CLASSES="ceramic DBTT" \
TARGET_PSI="-30 0 30" \
T_K=500 \
CRYSTAL_THETA_DEG=45 \
TRACTION_PROBE_RADIUS_M=1e-5 \
CAL_PSI_TOL_DEG=0.75 \
EVENT_PSI_TOL_DEG=2 \
MAX_CAL_ROOT_ITERS=20 \
MAX_CONTROL_ITERS=8 \
MAX_JOBS=1 \
RECALIBRATE=1 \
CALROOT=runs/mixed_mode_fem_czm_v7_preflight_cal \
OUTROOT=runs/mixed_mode_fem_czm_v7_preflight \
bash run_mixed_mode_fem_czm_v7_campaign.sh
```

Calibration can be slower than v6 because cancellation-sensitive targets may require several one-step adaptive-CZM probes. Follow the probe count in the printed calibration rows.

## Acceptance criteria

The calibration CSV must show for all three targets:

```
phase_converged = True
first_production_step_verified = True
phase_sample_reliable = True
```

The final campaign status should show:

```
event_phase_control_converged = True
```

and status `event` or `right_censored`. Results labeled `event_phase_mismatch` or `right_censored_phase_mismatch` must not enter the mixed-mode envelope.

## Important warning about negative remote opening

The reported +30° response matrix may require a very small negative remote opening coefficient while maintaining positive local crack-tip normal traction. V7 preserves and audits that state rather than changing it. Before long crack extension, any case with `loading_open_is_tensile=False` should be checked with crack-face contact enabled or with an alternative experimental geometry that generates the same local phase under nonnegative remote opening.

## Principal output files

Calibration:

```
mixed_mode_loading_calibration_v7.csv
mixed_mode_loading_calibration_history_v7.csv
production_backend_basis_v7.json
production_backend_probes/
```

Campaign:

```
campaign_status_v7.csv
mixed_mode_v7_anisotropic_all_cases.csv
<class>/psi_<angle>/mixed_mode_control_history_v7.csv
<class>/psi_<angle>/production_backend_control_final_summary.json
plots_v7/
```
