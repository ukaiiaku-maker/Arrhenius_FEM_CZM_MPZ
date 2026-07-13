# Mixed-mode FEM/CZM v8 — production-backend probe Boolean fix

V8 retains the v7 exact adaptive-CZM full-circle calibration and anisotropic calibrated-tip mechanics. It corrects a CSV Boolean decoding error that caused every valid phase probe to be discarded.

## Root cause

The production backend writes `traction_phase_probe_reliable=True`. During CSV aggregation, pandas may represent this as floating-point `1.0`. The v7 helper accepted `"1"` but rejected `"1.0"`. Consequently:

- `phase_sample_reliable` was false for every probe;
- the root finder ignored all exact-backend samples;
- 0 and +30 degree targets were rejected despite sub-degree errors;
- the -30 degree search returned its initial basis estimate after 23 unused probes.

V8 decodes bool, NumPy bool, integer, float, and text representations robustly. It also writes `phase_probe_flag_raw` and `phase_probe_flag_decoded` into the calibration history.

## Expected behavior

For each target, valid probes now participate in bracket/secant/root selection. A target passes when the exact adaptive-CZM first-step traction phase is within `CAL_PSI_TOL_DEG`.

The mechanics are otherwise unchanged from v7.
