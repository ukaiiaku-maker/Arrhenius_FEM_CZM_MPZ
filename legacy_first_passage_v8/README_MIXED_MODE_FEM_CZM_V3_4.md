# Mixed-mode FEM/CZM v3.4 campaign probe fix

This is a runner-only correction for the validated v3.3 mechanics and circular-phase calibration.

## Defect corrected

The v3.3 campaign performed an interface probe with:

```bash
python -m arrhenius_fracture.mixed_mode_first_passage_v3_3 --help
```

The mixed-mode driver requires `--target-mode-phase-deg`, so that probe exited before any physical case was launched. The actual case command already contained the required target argument.

V3.4 imports the module and verifies its `MODEL_ID` using `python -c`; it does not invoke the required CLI during the probe. Regression tests also verify that every physical-case argument list includes `--target-mode-phase-deg`.

## Dependencies

The v3.3 files must already be installed, including:

- `arrhenius_fracture/mixed_mode_first_passage_v3_3.py`
- `calibrate_mixed_mode_loading_v3_3.py`
- `tests/test_mixed_mode_first_passage_v3_3.py`

## Verify

```bash
CONDA_ENV=arrhenius-fem-czm bash verify_mixed_mode_fem_czm_v3_4.sh
```

## Resume the three-angle preflight

The validated calibration can be reused without recalculation:

```bash
CONDA_ENV=arrhenius-fem-czm \
CLASSES="ceramic" \
TARGET_PSI="-30 0 30" \
CALROOT=runs/mixed_mode_fem_czm_v3_3_circular_phase_preflight \
OUTROOT=runs/mixed_mode_fem_czm_v3_4_preflight \
bash run_mixed_mode_fem_czm_v3_4_campaign.sh
```

Set `RECALIBRATE=1` only when the mesh, geometry, or elastic settings have changed.
