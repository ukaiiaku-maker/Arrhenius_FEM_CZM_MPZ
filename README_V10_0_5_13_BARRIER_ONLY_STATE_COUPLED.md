# v10.0.5.13.3 barrier-only tip/source-only FEM/CZM campaign

## Purpose

Run the four selected response classes through the existing full 2-D `Arrhenius_FEM_CZM_MPZ` mechanics while changing only the Arrhenius barrier landscape.

The solver remains the inherited full 2-D FEM/CZM stack: plane-strain FEM, cubic anisotropic elasticity, domain J, adaptive-CZM topology and quality transactions, elastic continuum bulk under `tip_only`, the moving crack-tip MPZ source/state evolution, conservative rollback, and the validated 330/240/100 um mesh/J policy. No solver subsystem is removed or replaced.

## Barrier-only transfer

Only cleavage, emission, Peierls, and Taylor barrier surfaces and attempt frequencies are transferred. Candidate source count/density, source refresh, encounter/retention, shielding/back-stress, blunting, MPZ-grid recommendation, Taylor-correlation, and developed-state fields are provenance only.

## Tip/source-only policy

The campaign requires `bulk_plasticity_mode=tip_only`: elastic FEM bulk, moving crack-tip MPZ source interpretation, no uniform bulk mobile/retained state, and restart/completion verification of that policy.

## Rate-preserving macro-ramp

The historical rate uses `dU=2e-7 m`, `dt=8.4 s`. The launcher defaults to `dU=2e-5 m`, `dt=840 s`, preserving exactly the same `dU/dt`. The adaptive-event controller scales both by the same accepted trial fraction near events. The launcher rejects an accidental rate change unless `ALLOW_RATE_CHANGE=1` is explicitly supplied.

The 100x step is a screening acceleration, not assumed convergence. Compare the first DBTT 700 K result against a 50x case (`dU=1e-5`, `dt=420`) using first-event K, emitted/retained state at first event, and the 0--20 um R-curve before promoting 100x to the 40-case sweep.

## Campaign

Four options, 300--1200 K in 100 K increments, 100 um target extension, 45 degree orientation, straight single front, deterministic events, two concurrent cases, restart-safe verification, and live prefixed output.

## Validation

```bash
python -m pip install -e . --no-deps
python -m pytest -q \
  tests/test_v100513_barrier_only.py \
  tests/test_v1005131_preserved_state.py \
  tests/test_v1005132_startup_resolution_warning.py \
  tests/test_v1005133_tip_only_ramp.py \
  tests/test_v1005123_phase_c_repairs.py
bash -n run_v10_0_5_13_3_barrier_only_monotonic.sh
```
