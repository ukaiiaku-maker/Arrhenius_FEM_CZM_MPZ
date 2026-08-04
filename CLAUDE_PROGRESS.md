# Claude progress

- Updated: 2026-08-03 22:03 PDT
- Repository: /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude
- Branch: claude/v10.0.5.18.4.0-j-controlled-loading
- HEAD: ad06e4c4197f768e3e3f2b182228afafce2b703e at start of session
  (HEAD~1 = 293491157063484bb10df6adc479f847a6a08ba6, the verified source
  commit on v10.0.5.18.4.0-theta0-pf-parity; ad06e4c is only the local
  workspace-handoff commit adding CLAUDE.md/FEM_CZM_HANDOFF.md/etc.). A new
  commit for the Phase 1 trajectory reader follows this update — see
  "Commits created" below for the exact hash once made.
- Remote: origin = https://github.com/ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ.git
- Conda env: arrhenius-fem-czm-claude
  (python: /opt/homebrew/Caskroom/miniconda/base/envs/arrhenius-fem-czm-claude/bin/python)
- Package import path: arrhenius_fracture/__init__.py in this workspace
- Package version (authoritative, via `importlib.metadata.version('arrhenius-fem-czm')`): 10.0.5.18.4.0
  - NOTE: `arrhenius_fracture.__version__` (a hardcoded string at
    arrhenius_fracture/__init__.py:91) still reads '10.0.2' — stale dead
    code left over from an old commit ("Expose v10.0.2 package version"),
    never bumped through later pyproject.toml version bumps. Confirmed this
    predates the workspace (not introduced here; `git show ad06e4c` does not
    touch this file). Not physics-relevant; every real version gate
    (preflight/launcher scripts) checks `importlib.metadata`, not this
    attribute. Left unfixed — out of scope for the current objective.

## Current phase

Phase 0 (provenance/baseline) is **complete**. Phase 1 (PF driving-force
trajectory reader) is **complete**. Not yet started: Phase 2 (transactional
J-controlled loading controller).

## Phase 0 — completed findings

- PF reference case directory exists and is intact:
  `/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/steps_1000K.csv`
  — 1652 data rows, columns include `KJ_Pa_sqrtm`, `J_effective_direct_J_per_m2`,
  `J_signed_direct_J_per_m2`, `B`, `N_em`, `crack_extension_m`, `dt_cur_s`
  (constant 8.4 s through the prefracture regime),
  `W_bulk_plastic_cumulative_J_per_m`, etc. SHA-256 of this exact file:
  `666d839ae8ef4267ae4ac7ab52a4220de676525bfb70917744e25a119c291b0c`.
- PF kernel family.json SHA-256 matches exactly:
  `a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a`
  at `.../runs/v10_2_28_kernel_cache/1447653d.../family.json`.
- **PF prefracture/first-passage row located precisely**: row index 176
  (0-based; PF `step`=177) is the first row with `crack_extension_m>0`.
  At that row: KJ=53.9019 MPa√m, J≈6531 J/m² (matches FEM_CZM_HANDOFF.md §5
  "J≈6531 J/m², KJ≈54 MPa√m" first-passage figure). Rows 0–175 (176 rows,
  physical time 8.4s to 1478.4s) are the exact PF prefracture regime.
- Archived slow-ramp FEM baseline directory confirmed at:
  `/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_17_paper7_self_contained/runs/v10_0_5_18_4_0_peak1000_theta0_seed8666_pf_full_field_parity_v1_slow_displacement_ramp_baseline_20260803_202328/`
  — launcher.pid (31068) is no longer alive (deliberately stopped, not
  crashed — no traceback/error in logs). The true log tail goes further
  than the CLAUDE.md snapshot (step 3600): it reaches step 9400,
  KJ=11.367 MPa√m, B=0, N_em=0, a=0.500mm (still no growth) before being
  stopped. This confirms the "slow driving-force ramp" diagnosis in
  FEM_CZM_HANDOFF.md §5.
- **Independently reproduced the same slow ramp in this workspace**: an
  800-step local smoke run (see command below) produced KJ values
  bit-identical to the archived baseline at every printed step (step 200
  KJ=2.215, 400=3.111, 600=3.786, 800=4.346 MPa√m) — full determinism
  confirmed (same seed 8666, same kernel, same parameter row). Command used:
  ```
  conda run -n arrhenius-fem-czm-claude env \
    STEPS=800 PRINT_EVERY=200 SAVE_SNAPSHOTS=0 TARGET_EXT_UM=1000 \
    GENERATE_SOLVER_PLOTS=0 \
    CAMPAIGN_ROOT="$(pwd)/runs/phase0_slow_ramp_reproduction_smoke_20260803" \
    bash run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh
  ```
  Output retained at `runs/phase0_slow_ramp_reproduction_smoke_20260803/`
  (not committed; `runs/` is out of scope for git per CLAUDE.md).
- All files named in FEM_CZM_HANDOFF.md §4 exist in this workspace.
- Production entry point confirmed:
  `python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production`
  invoked via `run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`
  (STEPS/CAMPAIGN_ROOT/etc. are env-overridable; checks
  package_version==10.0.5.18.4.0 but does NOT hard-check branch name, so it
  runs fine on the dev branch).
  `scripts/preflight_and_launch_v10051840_theta0_full_field.sh` DOES
  hard-check `EXPECTED_BRANCH=v10.0.5.18.4.0-theta0-pf-parity`, which our
  dev branch will never match — do not weaken that guard; use the plain
  production launcher or pytest directly for dev-branch verification instead.
- Focused regression suite (Phase 0 required set) — 19/19 passed:
  tests/test_theta0_pf_parity_v10051840.py,
  tests/test_theta0_pf_full_field_parity_v10051840.py,
  tests/test_theta0_pf_full_field_metadata_v10051840.py,
  tests/test_corridor_refinement_metadata_v10051840.py,
  tests/test_atomic_path_corridor_v10051839.py,
  tests/test_atomic_path_corridor_node_compaction_v10051839.py.

## Phase 1 — completed: PF driving-force trajectory reader

New module: `arrhenius_fracture/pf_theta0_driving_force_trajectory_v10051840.py`
- `load_pf_theta0_driving_force_trajectory(steps_csv_path)` — reads a PF
  `steps_*.csv`, validates required columns, reconstructs physical time by
  cumulative-summing the per-row `dt_cur_s` (never inferred from step count
  alone), records SHA-256 provenance of the source file, and locates the
  first row with nonzero `crack_extension_m` as `first_passage_row_index`
  (+ step/time/J/KJ at that row). Raises on missing columns, non-positive
  dt, or non-increasing step — fails closed rather than silently continuing
  on malformed reference data.
- `prefracture_slice(trajectory)` — returns the strict prefix of rows
  before first passage. Per FEM_CZM_HANDOFF.md §6 ("no assignment of a
  crossing event to the interval endpoint"), the first-passage row itself
  is deliberately excluded from the slice; callers read its J/KJ from
  `first_passage_*` fields instead.
- `interpolate_target(trajectory, physical_time_s, channel="J"|"KJ")` —
  bounded linear interpolation; raises ValueError on any request outside
  `[time[0], time[-1]]` (extrapolation forbidden, matching the
  `kernel_extrapolation_allowed=False` convention already used for the PF
  signed shielding kernel elsewhere in this codebase).

Regression test: `tests/test_pf_theta0_driving_force_trajectory_v10051840.py`
— 7/7 passing, including one integration-style test guarded by
`pytest.mark.skipif` against the real PF reference CSV (verifies exact
SHA-256, exact first-passage row 176, KJ=53.9019 MPa√m at that row).

This reader is a pure function of the CSV; it does not touch FEM state,
hazard state, RNG, or mesh. It is the input Phase 2's transactional
controller will consume to compute `J_target(t)` at each candidate solver
step.

## Files changed this session

- New: `arrhenius_fracture/pf_theta0_driving_force_trajectory_v10051840.py`
- New: `tests/test_pf_theta0_driving_force_trajectory_v10051840.py`
- Modified: `CLAUDE_PROGRESS.md` (this file)
- Untracked, not committed: `runs/phase0_slow_ramp_reproduction_smoke_20260803/`
  (local smoke-run output; `runs/` stays out of git per CLAUDE.md)

## Commits created

Pending this update — the plan is to commit the two new files together as
`feat: add PF driving-force trajectory reader` (one of the milestones
CLAUDE.md suggests verbatim), including this CLAUDE_PROGRESS.md update.
Check `git log --oneline -3` for the actual resulting hash.

## Tests run

- 19/19 Phase 0 focused tests (see list above).
- 7/7 new Phase 1 trajectory-reader tests.
- 26/26 combined (Phase 0 set + Phase 1 test) run together — no
  interaction/import-order issues.

## Simulations run

- Read-only inspection of the archived slow-ramp baseline (no new run).
- One completed 800-step local smoke run (see Phase 0 section) — confirmed
  bit-identical KJ trajectory vs. the archived baseline; not run to first
  passage (would take ~145,000 steps at the current fixed-displacement
  ramp rate, which is precisely the defect Phase 2 will fix).

## Latest accepted physical state

No FEM/CZM solver source code has been modified yet — only a new read-only
trajectory-reader module was added. The archived + freshly reproduced
slow-ramp evidence (step 9400, KJ=11.4 MPa√m, B=0, N_em=0, no growth) is
the accepted Phase 0 baseline confirming the reported defect. The PF
prefracture target trajectory (176 rows, up to KJ=53.9 MPa√m) is now
available in-code as the Phase 2 controller's reference input.

## Latest numerical blocker

The fixed remote-displacement ramp (dU/dt matched to the PF rate) produces
a far slower J/KJ trajectory than the PF reference (~145,000 steps
extrapolated to first passage). This is a controller problem, not a
hazard/kernel/mesh problem — none of B, N_em, or crack extension have
moved yet at step 9400. Root cause per FEM_CZM_HANDOFF.md §5: FEM and PF
geometries have different compliance, so equal dU/dt does not produce
equal dJ/dt. Not yet fixed — that is Phase 2's job. No FEM/CZM solver code
has been touched to address it yet; only the read-only target-trajectory
input Phase 2 will consume has been built and tested.

## Failed approaches that should not be repeated

- Do not attempt to "fix" the slow ramp by increasing dU/dt or otherwise
  tuning displacement rate directly — this was the approach that produced
  the archived baseline problem in the first place, and CLAUDE.md forbids
  tuning around a defective loading controller. The correct fix is the
  transactional J_target(t)/KJ_target(t) controller (Phase 2).
- Do not run `scripts/preflight_and_launch_v10051840_theta0_full_field.sh`
  directly on the dev branch — it will hard-fail on the branch-name check
  by design (source-repo gate, not applicable here). Use the plain
  production launcher or pytest instead for local dev-branch verification.
- Do not run heredoc-piped `conda run -n ... python3 - <<'PY' ... PY`
  invocations expecting captured stdout when redirected to a file — output
  was silently empty (buffering interaction between `conda run` and stdin
  heredocs). Write the script to a real file and run
  `conda run -n ... python3 /path/to/script.py` instead; that works
  reliably.

## Exact next command

1. Design and implement Phase 2: the transactional trial-solve controller.
   Concretely: a function that, given the current accepted FEM/hazard/MPZ
   state and a target `(physical_time_s, KJ_target)` pair from
   `prefracture_slice(...)`, (a) predicts a trial boundary displacement via
   the elastic square-root predictor
   `U_new = U_old * sqrt(KJ_target / max(KJ_old, KJ_min))` (note: the
   handoff's predictor formula is written in terms of J; since
   `KJ ∝ sqrt(J)`, the equivalent KJ-based predictor is
   `U_new = U_old * (KJ_target / max(KJ_old, KJ_min))`, i.e. linear in KJ,
   not sqrt — verify which the codebase's existing J vs. KJ solve path
   expects before choosing), (b) solves the FEM system at the trial U,
   (c) compares achieved KJ against target within a tolerance, (d) applies
   a safeguarded secant/proportional correction and retries if outside
   tolerance, (e) on acceptance, integrates hazard/plasticity state over
   the actual elapsed physical time and localizes any stochastic crossing
   strictly inside the interval, (f) on rejection, restores exact prior
   state (bulk plastic strain, forest density, MPZ populations, B,
   thresholds, event-length draw, geometry, cohesive state, RNG) with zero
   net change. Before writing solver-loop-integration code, first locate
   where the existing production entry
   (`mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production.py`)
   currently applies the fixed `dU` ramp per step, and how it already
   captures/restores state on a rejected trial solve (the atomic cohesive
   commit/rollback machinery in `cohesive_trial_state.py` /
   `KineticCZMTransactionSnapshot` likely has the pattern to extend, per
   FEM_CZM_HANDOFF.md §6 — read that before designing a new snapshot
   mechanism from scratch).
2. Add a focused regression proving state is bit-identical before and
   after a rejected controller trial (no hazard/RNG/MPZ/plastic drift).
3. Wire the controller into a new or flagged production-entry variant
   (do not silently change the existing archived-baseline entry point's
   default behavior without a flag, so the Phase 0 slow-ramp reproduction
   remains re-runnable for comparison).

## Exact next scientific acceptance gate

Phase 3 prefracture parity for Peak / 1000 K / theta=0 / seed=8666: FEM
reaches PF's first-passage J≈6531 J/m² (KJ≈53.9 MPa√m) in a practical
number of controller-driven steps (not ~145,000), with B(t), N_em(t), and
bulk plastic work tracking the PF trajectory shape, before any
crack-growth qualification (Phase 4) is attempted.

## Continuation handoff (for a fresh session or Codex)

If context is exhausted before Phase 2 lands, everything needed to resume
is in this file plus the two new Phase-1 files. A fresh agent should: (1)
read this file, CLAUDE.md, and FEM_CZM_HANDOFF.md in full; (2) re-run
`git log --oneline -5` and `git status --short` to confirm no uncommitted
drift; (3) re-run the 26-test combined focused suite listed above to
re-confirm the accepted baseline before writing new code; (4) proceed
directly to the "Exact next command" section above — Phase 0 and Phase 1
do not need to be repeated.
