# v9.13 exact DBTT temperature shift

This diagnostic applies only the exact temperature-axis transformation to
`v912_targeted_local_peak_013476_0083`. It does not modify the cleavage shelf or
run any shelf-parameter variants.

The default targets move the current 1000 K peak to 700 K and 1100 K. The
runner evaluates the original candidate at 700--1200 K and each transformed
candidate at `lambda*T`, using the same loading map, CRN thresholds,
`translation_action_exponent=0.95`, and `max_hazard_increment=0.05`.

The requested extension defaults to 50 micrometres. The runner fails before any
simulation if it exceeds the loading-map coverage.

## Run

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v9_13_dbtt_temperature_shelf

git pull --ff-only
"$CONDA_PREFIX/bin/python" -m pip install -e .

REGISTRY=candidates/v9_12_targeted_local_4096_registry.csv

test -f "$REGISTRY" || {
  gzip -dc candidates/v9_12_targeted_local_4096_registry.csv.gz > "$REGISTRY"
}

LOADING_MAP="${LOADING_MAP:-$(
  find /Volumes/Data/Data/Nanopillar_calculation \
    -name v10_2_22_rcurve_loading_map.json \
    -print -quit
)}"

test -n "$LOADING_MAP" && test -f "$LOADING_MAP" || {
  echo "ERROR: calibrated v10.2.22 loading map was not found"
  exit 1
}

"$CONDA_PREFIX/bin/python" -m pytest -q \
  tests/test_dbtt_transform_v913.py \
  tests/test_v913_temperature_shift_entry.py

OUT=runs/v9_13_0083_temperature_shift_only_v1
mkdir -p "$OUT"

nohup /usr/bin/caffeinate -dimsu \
  "$CONDA_PREFIX/bin/python" -u \
  scripts/run_v913_dbtt_temperature_shift.py \
    --candidate-registry "$REGISTRY" \
    --candidate-id v912_targeted_local_peak_013476_0083 \
    --base-physics-json mpz_v9_13_v10222_transfer_common_physics.json \
    --loading-map "$LOADING_MAP" \
    --current-peak-temperature-K 1000 \
    --target-peak-temperatures-K 700 1100 \
    --base-temperatures-K 700 800 900 1000 1100 1200 \
    --target-extension-um 50 \
    --out "$OUT" \
  </dev/null > "$OUT/driver.log" 2>&1 &

PID=$!
echo "$PID" | tee "$OUT/driver.pid"
tail -f "$OUT/driver.log"
```

The run contains 18 R-curves: six original cases and six paired cases for each
of the two target peak temperatures.

## Outputs

- `temperature_scale_identity.csv`: paired K values and maximum differences for
  every event/state field.
- `temperature_scale_events.csv`: baseline and shifted event-resolved R-curves.
- `temperature_scaled_candidates.csv`: the two transformed active parameter
  rows.
- `temperature_shift_manifest.json`: loading-map coverage, numerical settings,
  maximum event-K difference, and pass/fail status.

The process exits with status 2 if the maximum paired event-K difference exceeds
`1e-8 MPa sqrt(m)`.
