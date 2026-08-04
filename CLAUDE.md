# CLAUDE.md

## Governing documents

Read this file and `FEM_CZM_HANDOFF.md` completely before editing source code.
`FEM_CZM_HANDOFF.md` is the authoritative scientific and numerical contract.

Maintain `CLAUDE_PROGRESS.md` as the durable state record for this workspace.
Update it after every completed milestone, before every commit, and before the
conversation context becomes limited.

## Authoritative workspace

Work only in:

`/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude`

Use only the Conda environment:

`arrhenius-fem-czm-claude`

The following existing repositories are read-only references:

- `/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_17_paper7_self_contained`
- `/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1`

Write run outputs only under this workspace's top-level `runs/`, unless a
specific diagnostic explicitly requires a new absolute output location.
Never alter or overwrite reference run directories.

Expected remote baseline:

- repository: `ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ`
- source branch: `v10.0.5.18.4.0-theta0-pf-parity`
- verified source commit at workspace creation:
  `293491157063484bb10df6adc479f847a6a08ba6`
- development branch:
  `claude/v10.0.5.18.4.0-j-controlled-loading`

## Primary objective

Develop the theta-zero FEM/CZM model until it can reproduce the qualitative
PF sharp-interface Arrhenius-hazard behavior for the four immutable production
parameterizations.

The immediate blocking task is not a broad parameter campaign. It is to replace
or supplement the slow fixed remote-displacement ramp with a transactional
PF-reference driving-force controller based on `J_target(t)` or `KJ_target(t)`.

First establish Peak, 1000 K, theta=0, seed=8666 prefracture parity. Then
qualify 20–50 micrometers of crack growth, then 400 micrometers, then 1000
micrometers. Only after those gates pass should work expand to the temperature
sweep and four parameterizations.

## Non-negotiable physics

- Fracture initiation and renewal remain Arrhenius hazard controlled.
- Do not add an independent athermal or empirical cohesive initiation criterion.
- The cohesive representation is geometry/separation numerics, not competing fracture physics.
- Preserve the stochastic cleavage threshold, emission thresholds, RNG state,
  event-length draw, selected direction, exact event endpoint, and physical
  event identity through all trial solves and retries.
- Numerical subsegments may represent one physical event but may not become
  multiple physical events.
- Geometry commit is atomic. Failure restores all constitutive, hazard,
  cohesive, mesh, MPZ, and RNG state.
- Do not weaken triangle-quality or local-area-ratio gates.
- Do not tune the four material rows independently.
- Do not replace or renormalize the audited PF active kernel.
- Preserve the distinction between moving tip MPZ state and full-field bulk
  Peierls–Taylor state.
- Rejected controller iterations must be completely transactional.
- Do not tune Arrhenius parameters to compensate for a defective J-integral,
  loading controller, MPZ translation, state transfer, or remesher.

## Working method

Proceed beyond planning into implementation, tests, focused local runs,
diagnosis, and evidence-based correction.

For each defect:

1. identify the earliest incorrect state transition;
2. add a focused regression where practical;
3. implement the smallest physics-neutral correction;
4. run the narrow test;
5. run the related inherited gates;
6. run a real production-entry smoke;
7. commit the passing milestone;
8. update `CLAUDE_PROGRESS.md`.

Do not stop after an audit or plan. Continue automatically to the next
qualified milestone.

Stop and ask for user input only when:

1. the handoff does not resolve a genuine physics ambiguity;
2. an operation would be destructive outside this isolated workspace;
3. credentials, pushing, merging, or remote mutation are required;
4. the same failure remains after two distinct evidence-based corrections;
5. a long production campaign is ready to launch but has not yet been approved.

Do not push, merge, rebase, force-reset, make PR #64 ready, or remove its draft
status. Local commits on the development branch are authorized.

## Immediate sequence

### Phase 0 — provenance and baseline

Report and record:

```bash
pwd
git rev-parse --show-toplevel
git branch --show-current
git rev-parse HEAD
git status --short
git remote -v
conda run -n arrhenius-fem-czm-claude \
  python -c "import sys, arrhenius_fracture; print(sys.executable); print(arrhenius_fracture.__file__)"
```

Verify the exact PF reference paths, parameter row, seed, kernel SHA-256,
archived slow-ramp FEM result, package version, and production entry point.

Run the focused tests from the handoff and reproduce the low-J slow-ramp rate.

### Phase 1 — PF driving-force trajectory

Read the prefracture PF `steps_1000K.csv` and create an audited target history.
Determine whether `J_target(t)` or `KJ_target(t)` is numerically more robust.

The target reader must preserve physical time intervals and provenance.

### Phase 2 — transactional controller

Implement a safeguarded controller with:

- elastic square-root predictor where appropriate;
- bounded secant/proportional correction;
- convergence tolerances and iteration limits;
- immutable accepted state during trial iterations;
- exact rollback of plastic, MPZ, hazard, threshold, event draw, cohesive,
  geometry, and RNG state after rejected trials;
- stochastic-event localization within a target interval;
- no assignment of a crossing event to the interval endpoint;
- no threshold redraw during retries.

### Phase 3 — prefracture qualification

For Peak, 1000 K, theta=0, seed=8666, compare against PF:

- target and achieved J/KJ versus physical time;
- B;
- N_em;
- bulk plastic work;
- first-passage time;
- first-passage J and KJ;
- controller iteration count and rejected-trial state equality.

Do not begin long crack growth until this passes.

### Phase 4 — short growth

Run 20 or 50 micrometers and audit:

- exact event length and endpoint;
- atomic geometry transaction;
- MPZ translation;
- B renewal;
- N_em evolution;
- J-contour recentering;
- cohesive work and energy balance;
- post-event continuation.

### Phase 5 — 400 and 1000 micrometers

Proceed only after short-growth acceptance. Preserve restart safety,
provenance, monitoring, and prior outputs.

### Phase 6 — temperature and four-class parity

Use the four immutable registry rows and common numerical treatment.
Do not use parameter-specific numerical workarounds.

## Git discipline

Keep commits narrow and descriptive. Before every commit:

```bash
git diff --check
git status --short
```

Suggested milestones:

- `test: reproduce theta0 slow driving-force ramp`
- `feat: add PF driving-force trajectory reader`
- `feat: add transactional J-controlled loading`
- `test: preserve state across rejected J-controller trials`
- `test: localize hazard events within J-target intervals`
- `feat: add theta0 prefracture parity audit`
- `fix: qualify theta0 short cohesive growth`
- `feat: add restart-safe theta0 parity campaign`

Never commit top-level `runs/`, large checkpoints, figures, videos, or archived
reference data.

## Durable progress record

`CLAUDE_PROGRESS.md` must always contain:

- date and local time;
- repository path;
- branch and HEAD;
- environment and import path;
- current phase;
- exact files changed;
- commits created;
- tests run and results;
- simulations run and results;
- latest accepted physical state;
- latest numerical blocker;
- failed approaches that should not be repeated;
- exact next command;
- exact next scientific acceptance gate.

Before context becomes limited, finish the current atomic edit if safe, run the
narrowest relevant check, commit passing work where appropriate, update the
progress file, and provide a concise continuation handoff.
