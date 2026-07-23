# v9.13 autonomous DBTT candidate search

This campaign searches the existing 4,096-row v9.12 targeted-local registry
with the autonomous 1-D model calibrated to the v10.2.22 two-dimensional
R-curves. It targets only the difficult DBTT/peak response. Ceramic-like and
weak-temperature parameterizations are outside this search.

## Why the existing pool is searched first

The registry is not a new global design. It contains local scrambled-Sobol
families around eight v9.12 parents:

- 2,048 plateau-family rows;
- 1,024 bridge-family rows; and
- 1,024 peak-family rows.

Eleven kinetic and state coordinates vary. The cleavage surface and several
geometry/constitutive quantities remain fixed. This is nevertheless the
correct first pool because it was generated around the prior DBTT and peak
basins. A new reduced-model candidate generator is justified only if the
calibrated autonomous 1-D search shows that this population contains no robust
peak.

## Objective and optimizer

No new DBTT objective is introduced. The driver calls the existing functions
in `augment_mpz_v9_12_directional_peak_targets.py`:

- low-temperature baseline: median response at \(T\leq700\) K;
- high-temperature response: median at \(T\geq1000\) K;
- intermediate peak window: 800–1000 K;
- peak rise relative to the pre-800 K baseline;
- peak drop relative to the final-temperature response; and
- peak prominence: the smaller of the peak rise and peak drop.

The established `peak_like_1d` threshold remains
1 MPa\(\sqrt{\mathrm m}\). For this campaign the response trajectory is
autonomous \(K(25\,\mu\mathrm m,T)\), replacing the obsolete prescribed-
protocol \(\Delta K_\mathrm{micro}(T)\).

The existing ExtraTrees training and acquisition scripts are reused after the
first evaluated wave. The acquisition is configured for 85% predicted
peak-quality candidates and 15% uncertainty/diversity exploration.

## First acquisition wave

The first 128 Sobol rows from each peak parent form a nested, low-discrepancy
256-candidate batch. Each candidate is evaluated at
700, 800, 900, 1000, 1100, and 1200 K to 25 µm. The screening integration uses
the calibrated event exponent \(p=0.95\) and a cleavage-hazard increment of
0.25. This is a five-times-coarser search approximation than the accepted
calibration setting; promoted candidates must be rerun with 0.05.

First check out this branch and refresh the editable installation in the active
conda environment. For a single-branch clone, first add the v9.13 branch to the
set that Git recognizes as remote branches:

```bash
BRANCH=v9.13-dbtt-4096-autonomous-search
git remote set-branches --add origin "$BRANCH"
git fetch origin

if git show-ref --verify --quiet "refs/heads/$BRANCH"; then
  git switch "$BRANCH"
  git merge --ff-only "origin/$BRANCH"
else
  git switch --no-track -c "$BRANCH" \
    "refs/remotes/origin/$BRANCH"
  git branch --set-upstream-to="origin/$BRANCH" "$BRANCH"
fi

"$CONDA_PREFIX/bin/python" -m pip install -e .

test -f scripts/run_v913_autonomous_dbtt_4096_wave1.sh
```

Then run from the repository root:

```bash
REGISTRY=candidates/v9_12_targeted_local_4096_registry.csv \
MAX_JOBS=4 \
OUT=runs/v9_13_autonomous_dbtt_4096_peak_wave1_v1 \
bash scripts/run_v913_autonomous_dbtt_4096_wave1.sh
```

The driver is case-resumable. Rerunning the same command reads completed case
JSON files and schedules only missing candidate/temperature pairs.

Primary outputs are:

- `case_results_checkpoint.csv`: one row per candidate and temperature;
- `R_curve_events.csv`: complete accepted event trajectories;
- `autonomous_dbtt_training_table.csv`: exact existing objective applied to
  \(K(25\,\mu\mathrm m,T)\);
- `ranked_candidates.csv`: evaluated candidates ranked by peak acceptance and
  prominence;
- `promoted_registry.csv`: the leading 48 evaluated rows;
- `autonomous_dbtt_surrogate.joblib`: refitted existing optimizer; and
- `next_active_registry.csv`: the next 256-row active-learning wave.

## Accurate promotion gate

After inspecting wave 1, rerun the promoted rows over the complete temperature
grid to 50 µm with the accepted integration resolution:

```bash
PYTHONPATH=. python -u scripts/run_v913_autonomous_dbtt_search.py \
  --candidate-registry \
    runs/v9_13_autonomous_dbtt_4096_peak_wave1_v1/promoted_registry.csv \
  --base-physics-json mpz_v9_13_v10222_transfer_common_physics.json \
  --loading-map \
    runs/v9_13_v10222_rcurve_targets_v1/v10_2_22_rcurve_loading_map.json \
  --policy-json mpz_v9_12_targeted_local_search_policy.json \
  --families \
  --per-parent 0 \
  --temperatures 300 400 500 600 700 800 900 1000 1100 1200 \
  --checkpoint-um 25 \
  --target-extension-um 50 \
  --translation-action-exponent 0.95 \
  --max-hazard-increment 0.05 \
  --jobs 4 \
  --promote-count 20 \
  --out runs/v9_13_autonomous_dbtt_wave1_promoted_full_v1
```

The final promotion decision must inspect the complete event-resolved
\(K(\Delta a,T)\) curves, not only the scalar peak score. The subsequent 2-D
campaign remains responsible for confirming that the post-peak branch is
plastic failure rather than a lower but still valid cleavage toughness.

## Decision after wave 1

- If robust peaks occur, continue active learning within this 4,096-row pool
  using `next_active_registry.csv`.
- If peaks occur only at 25 µm and disappear by 50 µm, use those failures as
  training data and expand the existing v9.12 local Sobol bounds.
- If no peak-like rows occur, construct a new proposal population with the
  existing v9.12 generator and optimizer. A separate phenomenological 0-D
  fracture model is not required.
