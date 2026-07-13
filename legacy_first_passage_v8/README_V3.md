# Stateful PD scratch-initiation model v3

This overlay updates the v2 discrete-site model. It contains no AT1 or AT2 calls.

## Install

Unzip over the root of the existing fatigue project after the v1/v2 stateful-PD files are present:

```bash
unzip -o stateful_pd_model_v3_overlay.zip
python -m unittest tests.test_stateful_pd_core
```

Six core tests should pass.

## Recommended verification run

```bash
OUT=runs/sn_stateful_pd_v3_700_pair \
STRESSES="700" \
CASES="no_shield shielded" \
MESH_SEED=1 \
PD_SEED=1 \
SNAPSHOT_EVERY=25 \
PRINT_EVERY=5 \
bash run_sn_stateful_pd_pilot.sh
```

The new spatial controls can be changed through:

```bash
PD_INITIATION_RADIUS=240e-6
PD_INITIATION_TAPER=60e-6
PD_INITIATION_BACK_EXTENT=60e-6
PD_AMPLIFICATION_DAMAGE_SCALE=0.05
```

## Files to inspect first

For each case inspect:

- `pd_initiation_diagnostics_final.png`
- `pd_state_final.npz`
- `summary.json`
- `sn_stateful_pd_history.csv`

Before calibrating event rates, confirm that:

1. coupling-shell candidate sites are zero;
2. the maximum birth intensity lies at or immediately ahead of the scratch root;
3. intact PD amplification is one;
4. bond softening remains zero until a realized stable defect exists.

## Calibration sequence

Do not begin by increasing the hit-memory time over several decades. First verify spatial localization. Then tune, in order:

1. physical site density;
2. crack-opening barrier and stress concentration;
3. hit count and delivery prefactor;
4. stabilization/healing competition;
5. growth and linkage kinetics.

Use multiple PD seeds only after the mean first-event probability is in a useful range.
