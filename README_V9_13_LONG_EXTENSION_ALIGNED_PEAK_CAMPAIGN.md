# v9.13 long-extension aligned-peak campaign

This campaign evaluates the 25 candidates that passed the archived 50 micrometre
`y__peak_like_1d` criterion, reruns them with a mechanically calibrated longer
loading map, measures peak stability with crack extension, and writes an exact
temperature-axis-transformed registry whose retained candidates peak near a
common target temperature.

The default campaign is:

- candidates: first 25 rows with `y__peak_like_1d=True` in the prior
  `ranked_candidates.csv`;
- temperatures: 700, 800, 900, 1000, 1100, 1200, and 1300 K;
- target extension: 100 micrometres;
- diagnostic checkpoints: 25, 50, 75, and 100 micrometres;
- common aligned peak: 900 K;
- long-peak estimator: discrete temperature-grid maximum.

## Required long loading map

The previous v10.2.22 loading map covers only approximately 52.08 micrometres
and is rejected for a 100 micrometre campaign. Do not tile or repeat it.

A new reference 2-D case must first reach at least the requested extension. The
extractor reads the completed case's:

- `steps_*K.csv`;
- `stochastic_avalanche_geometry_events.json`;
- `run_args.json`.

For each accepted event it records:

- `KJ/Uapp` as the event-preceding displacement-to-K geometry factor;
- the stored threshold action;
- `event_advance_m` as path translation;
- `x1-x0` as projected R-curve advance.

Example extraction:

```bash
"$CONDA_PREFIX/bin/python" \
  scripts/extract_v10222_long_rcurve_loading_map.py \
    --case-dir /absolute/path/to/completed/reference/T300K_th45_seed3621 \
    --reference-candidate-id v912_targeted_local_peak_013476_0083 \
    --minimum-coverage-um 100 \
    --out runs/v9_13_long_map/v10_2_22_long_rcurve_loading_map.json
```

The extractor fails if accepted steps and geometry events do not match, if the
CRN seed is inconsistent, if projected advance is nonpositive, or if the final
coverage is below the requested minimum.

## Run the complete campaign

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v9_13_dbtt_temperature_shelf

git fetch origin v9.13-long-extension-aligned-peak-campaign
git switch --track origin/v9.13-long-extension-aligned-peak-campaign
"$CONDA_PREFIX/bin/python" -m pip install -e .

SOURCE_RANKING=/absolute/path/to/\
v9_13_autonomous_dbtt_wave1_promoted_50um_v1/ranked_candidates.csv

LOADING_MAP=/absolute/path/to/v10_2_22_long_rcurve_loading_map.json

SOURCE_RANKING="$SOURCE_RANKING" \
LOADING_MAP="$LOADING_MAP" \
OUTROOT=runs/v9_13_peak25_long100um_align900K_v1 \
TARGET_EXT_UM=100 \
TARGET_PEAK_K=900 \
JOBS=4 \
bash scripts/run_v913_long_peak25_alignment.sh
```

The shell runner may extract the map automatically when `LOADING_MAP` is not
set:

```bash
SOURCE_RANKING="$SOURCE_RANKING" \
REFERENCE_CASE_DIR=/absolute/path/to/completed/reference/T300K_th45_seed3621 \
REFERENCE_CANDIDATE_ID=v912_targeted_local_peak_013476_0083 \
OUTROOT=runs/v9_13_peak25_long100um_align900K_v1 \
TARGET_EXT_UM=100 \
TARGET_PEAK_K=900 \
JOBS=4 \
bash scripts/run_v913_long_peak25_alignment.sh
```

## Stages

### 1. Preparation

`scripts/prepare_v913_long_peak25_campaign.py`:

- selects exactly 25 peak-like rows;
- normalizes the `x_raw__*` parameter columns into the active 29-field v9.13
  candidate contract;
- writes `selected_peak25_registry.csv`;
- verifies loading-map coverage before any R-curve is launched;
- records source and loading-map SHA-256 hashes.

### 2. Long autonomous screen

The established resumable `scripts.run_v913_autonomous_dbtt_search` module runs
25 candidates by seven temperatures, or 175 R-curves. It is invoked with
`python -m` to ensure the repository root is available for the existing
`scripts.*` imports.

### 3. Alignment analysis

`scripts/analyze_v913_long_peak_alignment.py` reads the event-resolved case JSON
files and uses strict checkpoint semantics. A checkpoint value is NaN unless the
case actually reached that crack extension.

For each checkpoint it reports:

- peak temperature and peak value;
- peak rise;
- post-peak drop;
- post-peak minimum;
- high-temperature rebound;
- whether the maximum is on a temperature-grid boundary.

Peak-temperature drift is evaluated across 25, 50, 75, and 100 micrometres:

- drift up to 50 K: `stable`;
- drift from 50 to 100 K: `moderate_drift`;
- drift above 100 K: `extension_dependent`.

A candidate is not aligned when any case is incomplete, a checkpoint is not
reached, the long peak lies on the temperature boundary, the post-peak drop is
below 1 MPa sqrt(m), or the peak is extension-dependent.

For each retained candidate, the exact transformation uses

```text
lambda_i = 900 K / T_peak,100um,i
```

and writes the transformed active parameter row. The aligned R-curve table is
derived by mapping each original temperature to `lambda_i*T`; redundant solver
runs are not required.

## Main outputs

Under `OUTROOT/preparation`:

- `selected_peak25_registry.csv`;
- `selected_peak25_registry.ids.txt`;
- `selected_peak25_registry.manifest.json`;
- optionally the extracted long loading map and audit CSV.

Under `OUTROOT/long_screen`:

- resumable per-case JSON files;
- `case_results_checkpoint.csv`;
- `R_curve_events.csv`;
- `ranked_candidates.csv`;
- run contract, progress, and search manifest files.

Under `OUTROOT/alignment`:

- `long_case_checkpoint_table.csv`;
- `checkpoint_peak_metrics.csv`;
- `candidate_long_peak_metrics.csv`;
- `aligned_registry_Tp900K.csv`;
- `aligned_case_results_derived.csv`;
- `alignment_rejections.csv`;
- `refinement_case_plan.csv`;
- `alignment_manifest.json`.

## Verification performed

- Python syntax checks for the new module and scripts;
- shell syntax check;
- unit tests for strict checkpoint handling, interior peak metrics, rebound,
  and peak-drift classification;
- preparation test selecting exactly 25 candidates from the archived ranking;
- fail-closed test rejecting a 52 micrometre map for a 100 micrometre target;
- synthetic end-to-end alignment test: a stable 1000 K peak produced
  `lambda=0.9`, while a boundary peak was rejected;
- extraction regression against the archived v10.2.22 reference case: all map
  arrays were reproduced to floating-point precision.
