# v10.0.5.13.3 barrier-only tip/source-only FEM/CZM campaign

## Purpose

Run the four selected response classes through the existing full 2-D
`Arrhenius_FEM_CZM_MPZ` mechanics while changing only the Arrhenius barrier
landscape.

The active solver remains the inherited full FEM/CZM stack:

- 2-D plane-strain FEM equilibrium;
- cubic anisotropic elasticity at 45 degrees;
- domain J integral;
- adaptive-CZM topology and geometry-quality transactions;
- elastic continuum bulk under the validated `tip_only` coupling;
- moving crack-tip MPZ source activity, emission, transport, retention,
  clearing, shielding, and crack-advance renewal;
- conservative rollback and existing completion/audit layers;
- 330 um physical refinement support with 240 um cluster J and 100 um local J.

No solver subsystem is removed or replaced.

## Barrier-only transfer contract

The response registry supplies only:

- crack-opening EXP-floor barrier;
- emission EXP-floor barrier;
- Peierls barrier, entropy, shape, and attempt frequency;
- Taylor barrier, entropy, shape, and attempt frequency.

Candidate-specific source count/density, source refresh, encounter/retention,
shielding, back-stress, blunting, MPZ-grid recommendation, Taylor-correlation,
and developed-state fields are recorded as provenance but are not applied.

## Tip/source-only policy

The campaign explicitly requires:

- `bulk_plasticity_mode = tip_only`;
- continuum bulk role = elastic FEM only;
- source interpretation = moving crack-tip MPZ only;
- no uniform bulk mobile/retained state;
- restart/completion verification of the same policy.

## Rate-preserving adaptive macro-ramp

The historical loading rate is defined by:

```text
dU = 2e-7 m
dt = 8.4 s
```

The v10.0.5.13.3 launcher defaults to a 100x numerical macro-step:

```text
dU = 2e-5 m
dt = 840 s
```

The physical rate `dU/dt` is unchanged. The existing adaptive-event controller
multiplies both values by the same accepted trial fraction near a kinetic event.
The launcher rejects an accidental loading-rate change unless
`ALLOW_RATE_CHANGE=1` is explicitly supplied for a rate study.

The 100x step is a screening acceleration, not assumed convergence. The first
DBTT 700 K case must be compared against a 50x case (`dU=1e-5`, `dt=420`) using
first-event K, emitted/retained state at first event, and the 0--20 um R-curve.
The 100x step is promoted to the 40-case sweep only when those results agree
within the campaign tolerance.

## Four barrier options

- `ceramic_primary`: `ceramic_restart02_candidate00`
- `weakT_primary`: `weakT_restart00_candidate00`
- `dbtt_primary`: `DBTT_restart04_candidate03`
- `peak_primary`: `DBTT_restart05_candidate61`

## Screening matrix

The full screening run contains 40 cases:

- four response options;
- 300 through 1200 K at 100 K increments;
- 100 um target extension;
- one straight active front, branching disabled;
- deterministic event statistics and no stochastic emission;
- two concurrent cases by default;
- three snapshots per case, with solver plots disabled to reduce disk use.

Long 500 um calculations should be selected after inspecting the 100 um map.

## Restart behavior

A case is skipped only when the requested target extension and barrier-only
manifest are complete, the barrier option matches, 330 um refinement support
is verified, and both the production and v9.11 integration audits confirm
`tip_only` with the moving crack-tip MPZ and no uniform bulk mobile/retained
state.

Any incompatible or interrupted case is archived by command/manifests/log tail,
removed, and rerun cleanly.

## New-install isolation

The version-specific launcher resolves its root from its own location and
refuses to launch unless `arrhenius_fracture` resolves inside that installation.

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
