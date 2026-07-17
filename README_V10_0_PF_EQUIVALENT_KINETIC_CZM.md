# v10.0 PF-equivalent kinetic CZM

Branch: `v10.0-pf-equivalent-kinetic-czm`

Reference: `ukaiiaku-maker/PF-fracture-fatigue`, branch
`v10.1.7.1-final-production-temperature-sweep`.

## Purpose

This branch adds a new selectable front-local state, `kinetic_campaign_czm`,
without changing the package-level `legacy_scalar` or `moving_pz` defaults.  It
ports the promoted PF material rows and the front-state mechanisms before any
new parameter optimization.

## Implemented foundation

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
- Predictor-corrector stepper API for progressive mechanics coupling.
- Isolated prescribed-K Stage-A trace exporter and fail-closed comparison.
- Stage-B abrupt Mode-I entry point using the mature FEM/CZM loop.
- No-artificial-controls audit.

## Deliberately blocked

The existing production `sharp_front.run_2d` loop advances topology only after
a renewal fires.  The required progressive lifecycle inserts an intact trial
interface before action develops, solves mechanics with that interface, and
then updates `damage=B` transactionally.  That dedicated single-front 2-D loop
is not yet promoted in this branch.

Accordingly:

- `mode_i_first_passage_v10_0 --czm-opening-coupling clock_linear` exits with an
  explicit error instead of silently using the old renew-then-open loop.
- No three-class 100 um progressive matrix runner is marked production-ready.
- No 500 um nine-temperature production runner is marked production-ready.
- The no-artificial-controls audit can require
  `full_progressive_trial_loop_active=true` and will fail until that gate is
  implemented and validated.

## Foundation checkout

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

This runs focused unit tests and writes isolated CZM state histories for all
three material classes.  It does not certify PF parity without PF reference
histories.

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
state and PF parameters.  It is a regression/parity gate, not the final model.

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

After a completed Stage-B run:

```bash
python audit_v10_0_no_artificial_controls.py \
  runs/v10_stageB_weakT_700K
```

Do not use `--require-progressive` for Stage B.  That option is reserved for the
future dedicated progressive loop and must fail until it is actually active.

## Non-negotiable implementation audit

The branch does not add:

- temperature-dependent source count;
- temperature-dependent shielding coefficient;
- empirical `N_max`;
- per-step emission cap;
- temporal source recycling;
- stored-energy subtraction from cleavage;
- independent cohesive failure criterion;
- AT2 fracture;
- wake shielding as the primary toughening mechanism;
- barrier re-fitting.

## Promotion order

1. Focused unit tests.
2. PF/CZM isolated front-state parity.
3. Abrupt CZM regression.
4. Dedicated progressive single-segment 2-D loop.
5. Penalty convergence.
6. Three-class, three-temperature, 100 um matrix.
7. Nine-temperature, 500 um production matrix.

No later stage should be launched merely because initiation magnitudes differ.
The developed state and normalized propagation response must also pass their
matched audits.
