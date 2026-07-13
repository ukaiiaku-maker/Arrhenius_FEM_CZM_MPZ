# Four-class 1x / 10x / 100x FEM-CZM rate sweep

This add-on runs the existing no-branch 500-um four-class R-curve sweep at three
loading-rate factors: 1x, 10x, and 100x.

The nominal displacement increment is held fixed and the nominal physical time
increment is scaled as

    dt(rate_factor) = BASE_DT / rate_factor.

With the default 4-mm specimen height, BASE_DU=2e-7 m and BASE_DT=8.4 s, the
nominal engineering strain rates are approximately:

- 1x:   5.95238e-6 1/s
- 10x:  5.95238e-5 1/s
- 100x: 5.95238e-4 1/s

The sharp-front adaptive event controller refines each proposed step by applying
the same `adaptive_frac` to both dU and dt. Therefore the imposed rate is
preserved under adaptive event refinement.

## Snapshot defaults

Snapshots are enabled by default:

- SAVE_SNAPSHOTS=12
- SNAPSHOT_COLS=6
- SNAPSHOT_BY_EXT_UM=50

For a 500-um target this gives extension-triggered coverage from initiation
through sustained growth. Each case should produce `field_snapshots_<T>K.png`.

## Dynamic timestep preflight

Before the full 120-case campaign, the wrapper by default runs a three-rate
preflight on weakT at 900 K to 25 um extension. Disable only when intentionally
resuming a previously validated campaign:

    PREFLIGHT_RATE_SMOKE=0 bash run_four_class_exp_floor_czm_rate_sweep.sh

## Adaptive audit

After the campaign, `audit_four_class_rate_sweep.py` writes:

- adaptive_timestep_audit.csv
- adaptive_timestep_audit_summary.json

The audit verifies rate preservation and reports:

- nominal and realized rate;
- maximum relative rate error;
- min / p01 / median / p99 adaptive fraction;
- min / median / max accepted dt;
- number and fraction of steps at the minimum adaptive fraction;
- maximum predicted adaptive clock increment;
- fraction of steps above the adaptive target;
- final crack extension;
- snapshot PNG presence.

Minimum-fraction hits are reported but are not automatically classified as
failures because genuinely unstable crack propagation is allowed to reach that
limit by design.

## Rate comparison outputs

The final summarizer writes:

- rate_temperature_summary.csv
- four_class_rate_effect_Kinit_vs_T.png
- four_class_rate_effect_Kprop_vs_T.png

The existing per-rate R-curve outputs and per-case snapshots are retained.
