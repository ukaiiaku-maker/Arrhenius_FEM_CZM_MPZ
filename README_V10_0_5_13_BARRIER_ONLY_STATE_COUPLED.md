# v10.0.5.13.1 barrier-only, state-coupled FEM/CZM campaign

## Purpose

Run the four selected response classes through the existing full 2-D
`Arrhenius_FEM_CZM_MPZ` mechanics while changing only the Arrhenius barrier
landscape.

The active solver remains the inherited full FEM/CZM stack:

- 2-D plane-strain FEM equilibrium;
- cubic anisotropic elasticity at 45 degrees;
- domain J integral;
- adaptive-CZM topology and geometry-quality transactions;
- state-coupled `bulk_same_pt_km` mobile/retained continuum state;
- moving process-zone transport, retention, release, recovery, shielding, and
  blunting;
- conservative rollback and existing completion/audit layers;
- productionized 330 um physical refinement support with 240 um cluster J and
  100 um local J.

No solver subsystem is removed or replaced.

## Barrier-only transfer contract

The response registry supplies only:

- crack-opening EXP-floor barrier;
- emission EXP-floor barrier;
- Peierls barrier, entropy, shape, and attempt frequency;
- Taylor barrier, entropy, shape, and attempt frequency.

The following candidate-specific reduced-model fields are recorded as
provenance but are not applied:

- source sites or source density;
- source refresh or time recovery;
- encounter/retention coefficients;
- shielding magnitude or clipping;
- back-stress weights;
- blunting coefficient;
- MPZ length/bin recommendation;
- Taylor-correlation closure values;
- developed-state initialization.

v10.0.5.13.1 also avoids assigning nominal common replacements for those
fields.  Their active values are whatever the existing full 2-D solver and
explicit CLI construct.  The campaign fixes only the common validation mode and
resolution: `bulk_same_pt_km`, 100 um MPZ, and 80 bins.

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

A case is skipped only when all of the following are verified:

- requested target extension was reached;
- barrier-only production manifest completed without exception;
- candidate state fields were not applied;
- 330 um physical refinement metadata was verified;
- explicit bulk mobile/retained state was active;
- at least one bulk-state update occurred.

Any other existing case directory is treated as an interrupted/non-resumable
partial case.  Its command, manifests, and final 200 log lines are archived in
`interrupted_case_logs/`, the partial directory is removed, and that case is
rerun from its initial state.  Completed compatible cases are retained.

## New-install isolation

`run_v10_0_5_13_barrier_only_monotonic.sh` resolves its project root from its own
location and refuses to launch unless the imported `arrhenius_fracture` package
resolves inside that same directory.  This prevents the historical shared
namespace from loading `PF-fracture-fatigue` or a deleted legacy worktree.

## Validation

```bash
python -m pip install -e . --no-deps

python -m pytest -q \
  tests/test_v100513_barrier_only.py \
  tests/test_v1005131_preserved_state.py \
  tests/test_v1005123_phase_c_repairs.py

bash -n run_v10_0_5_13_barrier_only_monotonic.sh
```

The complete inherited mechanics/state subset is also exercised by the
release-specific GitHub Actions workflow.
