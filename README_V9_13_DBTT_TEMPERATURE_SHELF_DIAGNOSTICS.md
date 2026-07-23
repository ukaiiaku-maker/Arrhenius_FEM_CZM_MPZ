# v9.13 DBTT temperature and lower-shelf diagnostics

This branch adds a narrow diagnostic for candidate
`v912_targeted_local_peak_013476_0083`. It does not retrain the surrogate,
change the calibrated shared 1-D reduction, or launch another broad search.

## 1. Exact temperature-axis translation

For a desired temperature scale `lambda = T_peak,target / T_peak,current`, the
candidate transform is

- `Tref_K <- lambda*Tref_K`;
- cleavage and emission `G00_eV <- lambda*G00_eV`;
- `peierls_H0_eV <- lambda*peierls_H0_eV`;
- `taylor_H0_eV <- lambda*taylor_H0_eV`;
- cleavage and emission `sT_GPa_per_K <- sT_GPa_per_K/lambda`;
- all `gT_eV_per_K`, activation entropies, attempt frequencies, EXP-floor
  shapes, fractional floors, source/state parameters, and geometry remain
  unchanged.

The runner recomputes paired R-curves and reports event-by-event differences.
This numerical identity check is retained because the absolute shared
`floor_min_eV` is not a candidate coordinate, even though the active fractional
floors dominate candidate 0083.

## 2. Lower-shelf diagnostics

Two cleavage-only perturbation families are evaluated.

### Anchored linear cleavage pivot

The cleavage zero-stress energy and characteristic stress are tilted about a
high-temperature anchor. At the anchor, the complete cleavage barrier is
identical at every stress. Emission, Peierls transport, Taylor completion,
source density, shielding, backstress, blunting, and geometry are unchanged.
This is the cleanest test of whether the lower shelf can be reduced while
preserving the high-temperature propagation increment exactly.

The transform is rejected if either the cleavage or emission zero-stress energy
or characteristic stress leaves the specified positive domain anywhere in the
evaluation temperature range. This prevents apparent improvements that rely on
constitutive positivity clamps.

### Global cleavage stress-axis scale

The cleavage `sigc0` and `sT` are multiplied by a common factor. This can lower
the opening shelf more strongly, but it no longer preserves the cleavage
barrier at the high-temperature anchor. It is included as a controlled test of
whether the Peierls/Taylor/shielding/backstress increment remains approximately
unchanged when the absolute opening scale is reduced.

## Run

Use a clean checkout of this branch and the calibrated v10.2.22 loading map.
The existing loading map covers 52.08184 micrometres, so the default target is
50 micrometres and the runner fails before simulation if the requested target
exceeds map coverage.

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v9_13_integrated_dbtt_search

git fetch origin v9.13-dbtt-temperature-shelf-diagnostics
git switch --track origin/v9.13-dbtt-temperature-shelf-diagnostics

"$CONDA_PREFIX/bin/python" -m pip install -e .

REGISTRY=candidates/v9_12_targeted_local_4096_registry.csv
if [[ ! -f "$REGISTRY" ]]; then
  gzip -dc candidates/v9_12_targeted_local_4096_registry.csv.gz > "$REGISTRY"
fi

OUT=runs/v9_13_0083_temperature_shelf_diagnostic_v1

"$CONDA_PREFIX/bin/python" -m pytest -q \
  tests/test_dbtt_transform_v913.py

"$CONDA_PREFIX/bin/python" -u \
  scripts/run_v913_dbtt_temperature_shelf_diagnostics.py \
  --candidate-registry "$REGISTRY" \
  --candidate-id v912_targeted_local_peak_013476_0083 \
  --base-physics-json mpz_v9_13_v10222_transfer_common_physics.json \
  --loading-map \
    runs/v9_13_v10222_rcurve_targets_v1/v10_2_22_rcurve_loading_map.json \
  --current-peak-temperature-K 1000 \
  --target-peak-temperatures-K 700 1100 \
  --shelf-temperature-K 700 \
  --anchor-temperature-K 1000 \
  --target-extension-um 50 \
  --out "$OUT"
```

## Primary outputs

- `temperature_scale_identity.csv`: paired original/scaled temperatures and
  maximum event-level differences;
- `temperature_scaled_candidates.csv`: transformed active candidate rows;
- `shelf_scan_aggregate.csv`: lower-shelf toughness, anchor toughness,
  `Delta K_50-first`, backstress, tip radius, and front width for each variant;
- `shelf_scan_temperature_summary.csv`: full temperature response;
- `shelf_scan_events.csv`: complete event-resolved R-curves;
- `shelf_scan_candidates_and_rejections.csv`: accepted transformed rows and
  explicit positive-domain rejection reasons;
- `diagnostic_manifest.json`: loading-map coverage and numerical settings.

A useful preservation criterion is that the anchor
`Delta K_50-first` and maximum backstress remain close to their baseline values
while the 700 K `K_50` decreases. The code reports both ratios directly; it does
not hard-code an acceptance threshold.
