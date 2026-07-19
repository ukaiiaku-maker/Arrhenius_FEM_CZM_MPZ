# v10.0 PF-equivalent kinetic CZM

Branch: `v10.0-pf-equivalent-kinetic-czm`

Reference: `ukaiiaku-maker/PF-fracture-fatigue`, branch
`v10.1.7.1-final-production-temperature-sweep`.

## Purpose

This branch adds a new selectable front-local state, `kinetic_campaign_czm`,
without changing the package-level `legacy_scalar` or `moving_pz` defaults. It
ports the promoted PF material rows and front-state mechanisms before any new
parameter optimization.

## Implemented

- Exact packaged PF manifests for `ceramic`, `weakT`, and `DBTT`.
- Explicit opening, cleavage, and emission stress channels.
- Active signed elastic shielding only in the cleavage channel.
- Fixed-scale local Taylor back stress only in the emission channel.
- Exact bounded finite-source emission.
- Exponential source refresh from crack advance only.
- Continuous cleavage action and MPZ translation inside one 5 um checkpoint.
- Strang-split plastic/cleavage integration with action and translation limits.
- One-checkpoint firing with `dt_consumed_s` and `dt_unused_s`.
- Developed-state cumulative and residence diagnostics.
- Trial cohesive state with `abrupt` and `clock_linear` mappings.
- Unilateral compressive contact retained at fully damaged interfaces.
- Full trial transaction snapshot and rollback, including the complete adaptive
  advance log and front-local kinetic state.
- Predictor-corrector progressive mechanics coupling.
- Guarded transformation of the actual production `sharp_front.run_2d` at three
  exact source anchors. The run aborts if an anchor changes or is non-unique.
- Progressive single-front Mode-I entry point that retains the v9.18.5.6
  triangle/area quality gate, edge-aware insertion, pre-refined corridor, and
  component-wise incremental x anchoring without importing the old v9.17
  renew-then-open controller.
- Isolated prescribed-K Stage-A trace exporter and fail-closed comparison.
- Stage-B abrupt Mode-I entry point using the mature FEM/CZM loop.
- No-artificial-controls audit.

## Validation status and blocked promotion

The progressive loop is implemented but has not yet passed a local repository
pytest run, a one-segment nonlinear FEM smoke, rollback/veto exercise, or penalty
convergence. It is therefore an experimental Stage-C gate, not production code.

Accordingly:

- no three-class 100 um progressive matrix is marked production-ready;
- no 500 um nine-temperature production runner is marked production-ready;
- a progressive run must write
  `full_progressive_trial_loop_active=true` and pass the no-artificial-controls
  audit;
- the first progressive run must remain one front, branching off, monotonic
  Mode I, one material, one temperature, and one 5 um segment.

## Checkout

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v9_18_5_8_state_coupled_material

git fetch origin \
  refs/heads/v10.0-pf-equivalent-kinetic-czm:refs/remotes/origin/v10.0-pf-equivalent-kinetic-czm

git worktree add \
  ../Arrhenius_FEM_CZM_MPZ_v10_0_pf_equivalent_kinetic_czm \
  -b v10.0-pf-equivalent-kinetic-czm \
  origin/v10.0-pf-equivalent-kinetic-czm

cd ../Arrhenius_FEM_CZM_MPZ_v10_0_pf_equivalent_kinetic_czm
```

## Foundation tests

```bash
conda activate arrhenius-fem-czm
CONDA_ENV=arrhenius-fem-czm bash run_v10_0_foundation_tests.sh
```

This compiles the transformed actual `sharp_front.run_2d`, runs the focused
unit tests, and writes isolated CZM state histories for all three material
classes. It does not certify PF parity without PF reference histories.

## Stage-A parity comparison

Generate a CZM trace:

```bash
python -m arrhenius_fracture.kinetic_campaign_parity \
  --material weakT \
  --history-csv prescribed_K_history.csv \
  --out runs/v10_stageA_weakT
```

Then compare against a PF v10.1.7.1 reference JSON generated from exactly the
same K/T/dt history:

```bash
python -m arrhenius_fracture.kinetic_campaign_parity \
  --material weakT \
  --history-csv prescribed_K_history.csv \
  --reference-json PF_REFERENCE_RECORDS.json \
  --out runs/v10_stageA_weakT_compared
```

The comparison checks clock action, micro advance, source budget, mobile and
retained populations, Taylor back stress, active shielding, blunting, and
cumulative state diagnostics.

## Stage-B abrupt regression

The Stage-B entry point keeps abrupt cohesive insertion but uses the new kinetic
state and PF parameters. It is a regression/parity gate, not the final model.

```bash
python -m arrhenius_fracture.mode_i_first_passage_v10_0 \
  --v10-material-class weakT \
  --czm-opening-coupling abrupt \
  --mode 2d \
  --temperatures 700 \
  --crystal-aniso --crystal-compete --no-crystal-branch --max-fronts 1 \
  --crack-backend adaptive_czm \
  --out runs/v10_stageB_weakT_700K
```

After completion:

```bash
python audit_v10_0_no_artificial_controls.py \
  runs/v10_stageB_weakT_700K
```

Do not use `--require-progressive` for Stage B.

## Stage-C progressive one-segment smoke

Do not increase the target until this exact gate passes:

```bash
ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM=5 \
ARRHENIUS_PREFINED_MODE_I_CORRIDOR=1 \
ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY=0.035 \
ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO=0.08 \
ARRHENIUS_MAX_TIP_H_OVER_DA=0.75 \
ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE=0.02 \
python -m arrhenius_fracture.mode_i_first_passage_v10_0_progressive \
  --v10-material-class weakT \
  --czm-opening-coupling clock_linear \
  --mode 2d \
  --temperatures 700 \
  --steps 50000 \
  --nx 36 --ny 72 \
  --tip-h-fine 1e-6 --tip-ratio 1.20 \
  --dU 2e-7 --dt 8.4 \
  --da-phys-um 5 \
  --target-crack-extension-um 5 \
  --crystal-aniso --crystal-compete --no-crystal-branch --max-fronts 1 \
  --crack-backend adaptive_czm \
  --mpz-length-um 100 --mpz-n-bins 200 \
  --out runs/v10_stageC_progressive_weakT_700K_5um_v1
```

Then require the progressive audit:

```bash
python audit_v10_0_no_artificial_controls.py \
  runs/v10_stageC_progressive_weakT_700K_5um_v1 \
  --require-progressive
```

The smoke is not passed unless all of the following are true:

- exactly one trial insertion and one committed event;
- monotonically increasing `B` and cohesive damage;
- total micro advance equals 5 um;
- `mpz_advance_on_commit_m=0`;
- no second topology event consumes `dt_unused_s`;
- no quality veto or unsupported cohesive endpoint;
- no non-finite mechanics state;
- the progressive and no-artificial-controls audits both pass.

## Non-negotiable implementation audit

The branch does not add:

- temperature-dependent source count;
- temperature-dependent shielding coefficient;
- empirical `N_max`;
- per-step emission cap;
- temporal source recycling;
- stored-energy subtraction from cleavage;
- independent cohesive failure criterion;
- smeared variational fracture;
- wake shielding as the primary toughening mechanism;
- barrier re-fitting.

## Promotion order

1. Focused unit tests and production-loop transform compilation.
2. PF/CZM isolated front-state parity.
3. Abrupt CZM regression.
4. Progressive one-segment smoke and rollback/veto audit.
5. Three-level normal/tangential penalty convergence.
6. Three-class, three-temperature, 100 um matrix.
7. Nine-temperature, 500 um production matrix.

No later stage should be launched merely because initiation magnitudes differ.
The developed state and normalized propagation response must also pass their
matched audits.
