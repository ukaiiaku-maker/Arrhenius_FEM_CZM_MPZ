# v9.13 autonomous 1-D R-curve calibration

This workflow makes the persistent-site 1-D model predict the same
stress-intensity-versus-crack-extension response measured in the
`v10_2_22_top5_dbtt_physical_width_50um_theta45_crn3621_v1` 2-D campaign.
It supersedes interpreting the old prescribed ten-segment protocol as a
complete R-curve prediction.

## Calibration boundary

The five candidate rows are immutable. Their cleavage, emission, Peierls,
Taylor, source-density, correlation, and `c_blunt` values are fingerprinted
when targets are built. A calibration run stops if any of these values has
changed.

Only quantities common to every candidate and temperature may be calibrated:

| Shared quantity | Role | Default bounds |
| --- | --- | ---: |
| `K_geometry_scale` | Reduced displacement-to-\(K\) normalization | 0.95–1.05 |
| `reference_source_area_scale` | 1-D representation of 2-D source area | 0.5–2 |
| `emission_geometry_scale` | Common slip-trace projection scale | 0.5–2 |
| `activation_to_line_content_scale` | Aggregate activation normalization | 0.5–2 |
| `persistent_backstress_scale` | Reduced Taylor-backstress mapping | 0.5–2 |
| `blunting_slip_fraction_scale` | Line-content-to-tip-radius mapping | 0.5–2 |
| `translation_action_exponent` | Event advance within the cleavage clock | 0.25–2 |

The current transfer already matches the archived constitutive state closely,
so the default calibration changes only `translation_action_exponent`. This
prevents a mechanical event-order discrepancy from being absorbed into the
physical candidate parameterizations.

## Target tables

`build_v913_v10222_rcurve_targets.py` extracts two optimization tables:

1. one row for each of 50 candidate/temperature models, containing
   \(K_\mathrm{first}\), \(K(10\,\mu\mathrm m)\),
   \(K(25\,\mu\mathrm m)\), \(K(50\,\mu\mathrm m)\), maximum backstress,
   minimum front width, and maximum tip radius;
2. one row for each of 800 accepted crack events, containing the complete
   \(K(\Delta a)\) curve.

Checkpoint values use the first accepted stochastic crack event at or beyond
the requested extension. They are not interpolated across a finite crack
jump.

The loading map is extracted once from the plastic-free 300 K reference
trajectory. It contains:

- the pre-event displacement-to-\(K\) geometry;
- the common-random-number cleavage thresholds;
- the path advance used for moving-frame state convection; and
- the projected advance used on the R-curve abscissa.

No evaluated case supplies its accepted 2-D \(K\), time-step, or crack-advance
history to the autonomous driver.

## Event operator

The original autonomous trial held the tip fixed throughout a cleavage
first-passage interval and convected the state only after the event. That
created artificial high-temperature retained populations and errors exceeding
20 MPa√m in the late R-curve.

The corrected `hazard_coupled` operator distributes the already sampled event
advance over its cleavage-action interval:

\[
f_a = \left(\frac{H}{H_c}\right)^p ,
\]

where \(H/H_c\) is the completed fraction of the event threshold and \(p\) is
one shared exponent. The emission geometry stays frozen at the pre-event
extension, preserving the 2-D update order while moving-frame convection and
constitutive kinetics are integrated together.

This is still autonomous: the only per-event inputs are the common stochastic
clock and the shared crack-geometry map, both known before evaluating a new
candidate.

## Reproduce the calibration

Build the targets:

```bash
PYTHONPATH=. python -u scripts/build_v913_v10222_rcurve_targets.py \
  --archive /path/to/v10_2_22_top5_dbtt_physical_width_50um_theta45_crn3621_v1.zip \
  --candidate-registry candidates/v9_13_persistent_sites_top5_registry.csv \
  --out runs/v9_13_v10222_rcurve_targets_v1
```

Fit the one shared event exponent on exact candidate/temperature pairs:

```bash
PYTHONPATH=. python -u scripts/calibrate_v913_rcurve_to_v10222_top5.py \
  --mode grid \
  --candidate-registry candidates/v9_13_persistent_sites_top5_registry.csv \
  --base-physics-json mpz_v9_13_v10222_transfer_common_physics.json \
  --targets runs/v9_13_v10222_rcurve_targets_v1/v10_2_22_rcurve_checkpoint_targets.csv \
  --event-targets runs/v9_13_v10222_rcurve_targets_v1/v10_2_22_rcurve_event_targets.csv \
  --target-manifest runs/v9_13_v10222_rcurve_targets_v1/v9_13_v10_2_22_rcurve_target_manifest.json \
  --loading-map runs/v9_13_v10222_rcurve_targets_v1/v10_2_22_rcurve_loading_map.json \
  --train-cases \
    v912_targeted_local_plateau_010759_0403:800 \
    v912_targeted_local_peak_013476_0162:1000 \
    v912_targeted_local_peak_013476_0368:900 \
    v912_targeted_local_peak_005518_0118:1200 \
  --optimize-parameters translation_action_exponent \
  --translation-exponent-grid 0.5 0.65 0.8 0.9 0.95 1.0 \
  --out runs/v9_13_v10222_rcurve_calibration_v1
```

The four training cases include both event-allocation bifurcations and smooth
plastic R-curves. The final phase of the command evaluates the selected
calibration on all 50 cases.

## Accepted transfer

The finite grid selects
`translation_action_exponent = 0.95`; every other shared scale remains exactly
one. With the candidate rows unchanged, the 50-case acceptance sweep gives:

| Comparison | Count | MAE (MPa√m) | RMSE (MPa√m) | Maximum error (MPa√m) |
| --- | ---: | ---: | ---: | ---: |
| Four checkpoints per case | 200 | 0.133 | 0.599 | 8.010 |
| Every accepted event | 800 | 0.131 | 0.367 | 8.010 |
| Checkpoints excluding one phase outlier | 199 | 0.093 | 0.194 | 0.869 |
| Events excluding the same phase outlier | 799 | 0.121 | 0.234 | 0.974 |

The sole large discrepancy is candidate `v912_targeted_local_peak_013476_0162`
at 1000 K. Its plastic-resistance jump occurs one finite stochastic event too
early in 1-D, producing an 8.010 MPa√m error at the 10 µm checkpoint. The same
curve returns to errors of 0.869 MPa√m at 25 µm and 0.783 MPa√m at 50 µm.
This is an event-phase limitation, not a persistent bias in the R-curve.

Per-candidate checkpoint errors over all ten temperatures are:

| Candidate suffix | MAE (MPa√m) | RMSE (MPa√m) | Maximum error (MPa√m) |
| --- | ---: | ---: | ---: |
| `0118` | 0.042 | 0.080 | 0.240 |
| `0162` | 0.305 | 1.289 | 8.010 |
| `0314` | 0.087 | 0.172 | 0.598 |
| `0368` | 0.069 | 0.157 | 0.639 |
| `0403` | 0.162 | 0.263 | 0.573 |

Thus the autonomous 1-D model is quantitatively consistent with the archived
2-D R-curves for this common-random-number campaign, with the explicit
single-event caveat above. The complete target and prediction tables are
versioned at:

- `runs/v9_13_v10222_rcurve_targets_v1/`; and
- `runs/v9_13_v10222_rcurve_alpha0p95_all50_v1/`.

## Use in the next parameter search

After acceptance, freeze the shared loading map, common physics, and calibrated
event exponent. A new broad search may then vary only the existing candidate
parameter fields. Each candidate is evaluated autonomously against the same
temperature grid and common-random-number event sequence. Screening and
promotion should use:

- the complete event-resolved \(K(\Delta a,T)\) response;
- the four checkpoint values;
- backstress, front-width, and tip-radius state diagnostics; and
- explicit failure-mode classification for high-temperature plastic failure.

This ordering separates model calibration from material parameter search: the
1-D/2-D reduction is fixed first, and candidate kinetics are searched only
after the transfer error is quantified.
