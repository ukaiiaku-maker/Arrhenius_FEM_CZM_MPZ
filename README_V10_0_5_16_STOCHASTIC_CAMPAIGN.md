# v10.0.5.16 stochastic PF-parity campaign

## Production runner

```bash
bash run_v10_0_5_16_0118_300_1200K_200um_stochastic_family_campaign.sh
```

The default campaign revision writes to:

```text
runs/v10_0_5_16_0118_300_1200K_200um_stochastic_pf_parity_v2
```

## Temperature-specific stochastic streams

Each temperature receives one explicit 32-bit seed derived from SHA-256 of:

```text
v10.0.5.16 | BASE_HAZARD_SEED | temperature_K
```

The seed map is written before simulation to:

```text
temperature_seed_map.csv
```

The runner fails before starting if any derived seeds collide. The campaign
postprocessor also fails if two temperatures have the same seed or the same
accepted threshold/event-length sequence fingerprint.

## Snapshots

Visual solver plots are enabled by default. The runner requests:

- state snapshots at every accepted cleavage event;
- extension-triggered snapshots every 10 micrometres by default;
- 20 target snapshots arranged in five columns by default;
- serialized signed MPZ snapshots for restart and accounting.

Set `GENERATE_SOLVER_PLOTS=0` only for a diagnostic run that intentionally omits
visual snapshot figures.

## Final analysis

After all temperature cases complete, the runner calls:

```bash
python scripts/analyze_v100516_stochastic_campaign.py \
  --campaign-root "$CAMPAIGN_ROOT" \
  --target-extension-um "$TARGET_EXT_UM"
```

Outputs are placed in `final_analysis/`:

- `R_curves_all_temperatures.png` and `.pdf`;
- individual PNG/PDF R-curves for each temperature;
- `K_vs_temperature_initial_mean_max_end.png` and `.pdf`;
- `event_advance_vs_crack_extension.png` and `.pdf`;
- `R_curve_event_points_all_temperatures.csv`;
- `K_temperature_and_seed_summary.csv`;
- `temperature_seed_map.csv`;
- `analysis_audit.json`.

The toughness metrics are defined as:

- **initial:** `K_J` at the first accepted cleavage event;
- **mean:** arithmetic mean of `K_J` over all accepted cleavage events;
- **maximum:** maximum `K_J` over accepted cleavage events;
- **endpoint:** `K_J` at the final accepted cleavage event.

Only rows with `n_fire > 0` are plotted as R-curve event points, and their count
must equal the accepted adaptive-CZM geometry-event count.

## Refinement metadata lifecycle correction

The completed revision uses a verified captured refinement mesh when the final
runtime mesh object exists but no longer carries the production refinement
annotations. This changes only the post-run audit source; mesh generation,
mechanics, constitutive evolution, and quality gates are unchanged.
