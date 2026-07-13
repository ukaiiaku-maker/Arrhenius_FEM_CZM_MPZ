# Arrhenius FEM/CZM MPZ v9.1 — three-class 1-D calibration

## Purpose

Version 9.1 uses the moving-process-zone constitutive model introduced in v9.0 to design three unified material parameterizations:

1. **ceramic** — decreasing first-passage toughness with temperature and essentially no rising R-curve at any temperature;
2. **weakT** — intermediate, nearly temperature-independent first-passage toughness with a moderate rising R-curve that reaches a temperature-insensitive developed resistance;
3. **DBTT** — a brittle low-temperature shelf and a tough high-temperature branch, with the transition retained during repeated crack growth.

The peak-shaped class is deliberately excluded.

The fit does not expose the legacy scalar shielding coefficient, scalar emission saturation, per-step emission cap, or stored-energy cleavage offset. Monotonic fracture, fatigue, dwell, mixed mode, and the full FEM/CZM solver continue to use the same fitted material row.

## What is fitted

The calibration is staged to reduce parameter compensation.

### Stage 1: first passage

Fits the elementary crack-opening and emission surfaces:

- cleavage zero-stress barrier and temperature slope;
- cleavage stress scale;
- emission zero-stress barrier and temperature slope;
- emission stress scale.

The DBTT objective additionally requires a low-to-high-temperature switch in the process-zone state at first passage. A direct cleavage-temperature slope alone cannot satisfy the complete DBTT objective.

### Stage 2: repeated-growth/process-zone response

With Stage 1 barriers fixed, fits:

- finite source-site inventory;
- glide barrier;
- trapping and detrapping barriers;
- retained-state recovery barrier;
- source recovery rate;
- source refresh length;
- process-zone length;
- slip-based blunting coefficient.

### Stage 3: joint polish

Jointly adjusts the complete parameter vector against both first passage and repeated growth. This is a local Powell polish initialized from the staged results, not a new unconstrained global search.

## R-curve observables

Each reduced simulation follows 0–1000 micrometres of crack extension. Individual event values are retained, but the objective uses robust summaries:

- `K_init`: first crack-renewal event;
- `early_rise_per_100um`: median of all pairwise slopes within the early window;
- `K_plateau`: median event resistance in the 700–1000 micrometre window;
- `plateau_rise_per_100um`: robust terminal slope in that same window;
- `delta_KR = K_plateau - K_init`;
- a robust saturating R-curve fit, recorded as a diagnostic rather than used as the primary target.

The plateau is accepted only when the terminal slope is also near zero. Therefore a still-rising curve cannot obtain a good score merely by passing through the desired late-window median.

## Class constraints

### Ceramic

- monotonically decreasing `K_init(T)`;
- monotonically decreasing developed resistance;
- `delta_KR` between 0 and 1.5 MPa sqrt(m);
- greater high-temperature plastic activity is allowed, but it may not raise the resistance curve appreciably.

### Weak-temperature/FCC-like

- total first-passage variation no greater than approximately 2.5 MPa sqrt(m);
- total plateau variation no greater than approximately 3 MPa sqrt(m);
- `delta_KR` between 3 and 7 MPa sqrt(m);
- finite early rise followed by a near-zero terminal slope.

### DBTT

- monotonically increasing first-passage and developed resistance;
- low-temperature `delta_KR` between 0 and 1.5 MPa sqrt(m);
- high-temperature `delta_KR` between 7 and 14 MPa sqrt(m);
- large low-to-high-temperature separation in both first passage and plateau;
- weak shielding/blunting on the low-temperature shelf;
- developed shielding or blunting before the high-temperature first event and through its late R-curve.

## Design targets

The editable target table is:

```text
mpz_three_class_design_targets.csv
```

These values are **synthetic design anchors**, not claimed experimental measurements. They define the requested response topology and an initial absolute toughness scale. Modify the CSV before the production fit when a different absolute scale is preferred; no code change is required.

The initial material rows are:

```text
mpz_three_class_initial_guesses.csv
```

They are marked `INITIAL_GUESS_NOT_CALIBRATED` and are intentionally not publication parameterizations.

## Installation and verification

Use the isolated project environment:

```bash
conda activate arrhenius-fem-czm
cd Arrhenius_FEM_CZM_MPZ_v9_1_three_class_tuning
python -m pip install -e .
pytest -q
bash verify_mpz_v9.sh
```

The supplied package passes 142 tests and 16 subtests.

## Recommended run sequence

All stages are restartable and write to separate directories under `runs/mpz_v9_1_three_class_tuning`.

### 0. Evaluate the inherited guesses

```bash
STAGE=smoke bash run_mpz_three_class_tuning.sh
```

This is expected to fail the weak-T and DBTT design objectives. Its purpose is to verify the integration, output tables, and plots.

### 1. Fit first passage

```bash
STAGE=first \
POPSIZE=6 \
MAXITER=25 \
bash run_mpz_three_class_tuning.sh
```

This stage samples 300, 500, 700, 800, 900, 1100, and 1200 K so the DBTT transition is not hidden between sparse anchors.

### 2. Fit repeated growth

```bash
STAGE=rcurve \
POPSIZE=5 \
MAXITER=20 \
bash run_mpz_three_class_tuning.sh
```

This stage uses 300, 700, 900, and 1100 K and a coarse 25 micrometre crack increment for the global search. It follows each candidate to 1000 micrometres.

### 3. Joint polish

```bash
STAGE=joint \
MAXITER=20 \
MAXFEV=1200 \
bash run_mpz_three_class_tuning.sh
```

The joint stage uses a 20 micrometre increment and includes both 800 and 900 K to resolve the transition.

### 4. Numerical convergence

```bash
STAGE=verify bash run_mpz_three_class_tuning.sh
```

The audit evaluates all temperatures from 300–1200 K and varies:

- nominal loading increment `dK`;
- process-zone bin count;
- physical crack increment.

The final table should not be promoted to the FEM/CZM validation stage unless first passage, plateau resistance, and R-curve increment are stable under these changes.

## Main outputs

Each optimization stage writes:

```text
mpz_three_class_parameters.csv
mpz_three_class_metrics.csv
mpz_three_class_event_Rcurves.csv
mpz_three_class_score_components.csv
mpz_three_class_predictions.json
three_class_target_comparison.png
three_class_representative_Rcurves.png
class_results/<class>/...
checkpoints/<class>_best.json
```

The event table includes pre-renewal shielding, blunting, mobile and retained content, source availability, local slip, and cumulative emission. This is necessary to determine whether a target response was obtained through a plausible evolved process zone.

## After convergence

Use the same final parameter table for:

```bash
python run_mpz_fatigue_matrix.py \
  --parameters runs/mpz_v9_1_three_class_tuning/03_joint/mpz_three_class_parameters.csv \
  --classes "ceramic weakT DBTT" \
  --require-fitted \
  --out runs/mpz_v9_1_three_class_fatigue
```

```bash
python run_mpz_dwell.py \
  --parameters runs/mpz_v9_1_three_class_tuning/03_joint/mpz_three_class_parameters.csv \
  --classes "ceramic weakT DBTT" \
  --require-fitted \
  --out runs/mpz_v9_1_three_class_dwell
```

Then use the small FEM/CZM validation matrix before any full production sweep:

```bash
python run_mpz_fem_czm_validation_matrix.py \
  --parameters runs/mpz_v9_1_three_class_tuning/03_joint/mpz_three_class_parameters.csv \
  --classes "ceramic weakT DBTT" \
  --temperatures "300 800 1100" \
  --target-ext-um 100 \
  --require-fitted \
  --out runs/mpz_v9_1_three_class_fem_czm_validation
```

The validation protocol may use one active front to isolate constitutive persistence, but the underlying solver still retains anisotropy, mixed mode, branching, multiple fronts, coalescence, retirement, cyclic mechanics, restart, and snapshots. A branch-enabled check follows only after the single-front constitutive screen passes.
