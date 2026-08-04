# Claude progress

- Updated: 2026-08-04 (session continuation — real-engine controller
  audit output added; see `FEM_PF_PARITY_SCORECARD.md` for the scientific
  status ledger, which now organizes ongoing work per explicit user
  instruction. This file remains the narrative/commit/environment record.)
- Repository: /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude
- Branch: claude/v10.0.5.18.4.0-j-controlled-loading
- Development workflow: all source-code development is version-controlled
  in Git on this branch. Local commits are authorized and ongoing; nothing
  has been pushed yet (push is a separate, explicit, user-approved step —
  see "Exact next command" below for the precise push command). The local
  `arrhenius-fem-czm-claude` Conda environment runs the editable install
  from this exact local Git checkout — confirmed via
  `importlib.metadata.version('arrhenius-fem-czm')`.
- Completed commits on this branch this session, oldest first:
  - `7562911` — `feat: add PF driving-force trajectory reader` (Phase 1)
  - `3cc67bb` — `feat: add safeguarded J/KJ-target loading controller` (Phase 2a)
  - `ad2c446` — `feat: wire default-off PF KJ-target controller` (Phase 2b checkpoint)
  - `040610f` — `docs: add FEM/PF parity scorecard and PF artifact regeneration contract`
  - `99f6e5c` — `test: add real-engine controller audit output and integration tests`
  (HEAD before this session's work: `ad06e4c`, the workspace-handoff
  commit; `ad06e4c~1` = `293491157063484bb10df6adc479f847a6a08ba6`, the
  verified source commit on `v10.0.5.18.4.0-theta0-pf-parity`.)
- Remote: origin = https://github.com/ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ.git
- Conda env: arrhenius-fem-czm-claude
  (python: /opt/homebrew/Caskroom/miniconda/base/envs/arrhenius-fem-czm-claude/bin/python)
- Package version (authoritative, via `importlib.metadata.version('arrhenius-fem-czm')`): 10.0.5.18.4.0
  - NOTE: `arrhenius_fracture.__version__` (hardcoded at
    `arrhenius_fracture/__init__.py:91`) still reads `'10.0.2'` — stale
    dead code predating this workspace, not physics-relevant, not fixed
    (out of scope). Every real version gate uses `importlib.metadata`.

## EXTERNAL BLOCKER — the immutable PF reference artifacts are not
## currently trustworthy at their documented path. Read this before any
## flag-ON validation or Phase 3 work.

The read-only PF reference repository
(`/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1`,
outside this workspace) is being actively written to by a concurrent,
unrelated process — not this session, not this workspace. This session
never wrote to that repo (only `ls`/`shasum`/`wc`/read-only Python CSV
reads); `git grep`/`git status` in *this* workspace show no trace of it.

**Exact expected values (from FEM_CZM_HANDOFF.md / CLAUDE.md, verified
earlier this session before the external changes began):**

```text
PF case directory:
  PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/
    v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1/
    v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/steps_1000K.csv
  steps_1000K.csv SHA-256:
    666d839ae8ef4267ae4ac7ab52a4220de676525bfb70917744e25a119c291b0c

PF kernel:
  PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/v10_2_28_kernel_cache/
    1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json
  family.json SHA-256:
    a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a
```

**Current observed state (as of 2026-08-04 08:36-09:09 PDT):**
- The `steps_1000K.csv` at the exact path above no longer exists. It
  existed and was used successfully twice earlier this session (the Phase
  0 baseline reproduction smoke, and a flag-off regression re-run
  immediately after the `sharp_front.py` edit) — both are valid, already
  gathered evidence; only *later* re-runs against this path are blocked.
- `runs/` in that repo now contains newly created files from an apparently
  live, unrelated "v10.4.x bulk plasticity campaign"
  (`LATEST_v10_4_4/5/6_bulk_plasticity_campaign.txt`,
  `v10_4_6_theta0_rate1x_four_class_500um_bulk_plasticity_base3621_v1/`).
- The kernel file still exists at its documented path, but its **content
  has changed**: observed SHA-256 is now
  `d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3`
  (≠ the expected `a85b57ad9e...` above). The hash-named cache directory
  is being reused/overwritten by the concurrent process.

**Why this matters / what it does not affect:**
- The existing SHA-256 guards in `run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`
  and `scripts/preflight_and_launch_v10051840_theta0_full_field.sh`
  correctly fail closed on this mismatch — confirmed this checkpoint did
  **not** weaken either guard (see "Verification evidence" below). Nothing
  has silently run against wrong physics.
- Nothing in this workspace's own Git history is affected. Phases 0/1/2a
  are fully committed and independent of the external repo's live state.
  Phase 1's real-reference integration test is
  `pytest.mark.skipif`-guarded on the file's existence, so it now skips
  (not fails) — confirmed in this checkpoint's test run (36 passed, 1
  skipped).
- Only **flag-ON physical qualification** of the Phase 2b controller (and
  all of Phase 3) is blocked — the flag-off/default production path does
  not depend on this file at all for its own correctness, only the
  launcher script's own pre-existing precondition check and its final
  `compare_v10051840_theta0_pf_parity.py` call require the file to exist
  to run *at all* (flag on or off) via that particular script.

**Requirement before any further flag-ON work**: the authoritative PF
`steps_1000K.csv` and kernel `family.json` must be recovered or
regenerated, then copied into a **dedicated, immutable local reference
directory inside this workspace's own control** (e.g.
`reference_inputs/pf_v10_4_1_theta0_peak1000K_seed8666/`, not committed —
large reference data — but stable and outside the live, externally-mutable
`PF-fracture-fatigue_.../runs/` tree), and independently re-verified
against the two SHA-256 values above before being trusted again. Do not:
- copy/"rescue" data from the live external repo without explicit
  direction (it's someone else's active work area);
- accept a same-named-but-different-hash artifact as a substitute;
- retry against the live path repeatedly hoping it stabilizes.

## Current phase

Phase 0 (provenance/baseline): **complete**. Phase 1 (PF driving-force
trajectory reader): **complete, committed at `7562911`**. Phase 2a
(standalone controller decision logic): **complete, committed at
`3cc67bb`**. Phase 2b (wiring into `sharp_front.run_2d` + launcher
plumbing): **default-off integration complete and checkpointed at
`ad2c446`; controller audit-output + real-engine integration tests added
at `99f6e5c`** (real mesh, real solve, synthetic target — proves the
controller tracks a target KJ(t) correctly against the actual FEM engine,
not just synthetic callbacks). **Flag-ON qualification against the real,
audited PF reference data is still pending recovery of the immutable PF
reference artifacts** (see EXTERNAL BLOCKER below and
`PF_REFERENCE_REGENERATION_CONTRACT.md`). Do **not** describe the
controller as physically qualified against PF — only its wiring,
default-off safety, and real-engine numerical correctness against a
synthetic target have been verified. This project's scientific status is
now tracked primarily in **`FEM_PF_PARITY_SCORECARD.md`** (gate ladder +
per-quantity table); consult that file first for "what's actually been
proven," this file for "how/when/which commit."

## Phase 0 — completed findings

- PF reference case verified intact (before the external changes): 1652
  data rows, columns include `KJ_Pa_sqrtm`, `J_effective_direct_J_per_m2`,
  `B`, `N_em`, `crack_extension_m`, `dt_cur_s` (constant 8.4 s through
  prefracture). SHA-256 recorded above.
- PF kernel SHA-256 verified matching (before the external changes),
  recorded above.
- PF prefracture/first-passage row located precisely: row index 176
  (0-based; PF `step`=177) is the first row with `crack_extension_m>0`.
  At that row: KJ=53.9019 MPa√m, J≈6531 J/m² (matches FEM_CZM_HANDOFF.md §5
  "J≈6531 J/m², KJ≈54 MPa√m"). Rows 0-175 (176 rows, physical time 8.4s to
  1478.4s) are the exact PF prefracture regime.
- Archived slow-ramp FEM baseline confirmed at
  `Arrhenius_FEM_CZM_MPZ_v10_0_5_17_paper7_self_contained/runs/v10_0_5_18_4_0_peak1000_theta0_seed8666_pf_full_field_parity_v1_slow_displacement_ramp_baseline_20260803_202328/`
  (a different, unaffected reference repo). Its log reaches step 9400,
  KJ=11.367 MPa√m, B=0, N_em=0, a=0.500mm (still no growth) before being
  deliberately stopped — confirms the "slow driving-force ramp" diagnosis
  in FEM_CZM_HANDOFF.md §5.
- Independently reproduced the same slow ramp in this workspace: an
  800-step local smoke matched the archived baseline bit-for-bit at every
  printed step (step 200 KJ=2.215, 400=3.111, 600=3.786, 800=4.346 MPa√m).
- All files named in FEM_CZM_HANDOFF.md §4 exist in this workspace.
  Production entry point:
  `python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production`,
  invoked via `run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`.
  `scripts/preflight_and_launch_v10051840_theta0_full_field.sh` hard-checks
  `EXPECTED_BRANCH=v10.0.5.18.4.0-theta0-pf-parity` — do not use it on this
  dev branch; use the plain launcher or pytest directly instead.
- Focused regression suite (Phase 0 required set): 19/19 passed.

## Phase 1 — completed: PF driving-force trajectory reader (`7562911`)

`arrhenius_fracture/pf_theta0_driving_force_trajectory_v10051840.py`:
- `load_pf_theta0_driving_force_trajectory(steps_csv_path)` — reads a PF
  `steps_*.csv`, validates required columns, reconstructs physical time by
  cumulative-summing per-row `dt_cur_s` (never from step count alone),
  records SHA-256 provenance, locates `first_passage_row_index` (first
  `crack_extension_m>0`). Fails closed on missing columns/non-positive
  dt/non-increasing step.
- `prefracture_slice(trajectory)` — strict prefix before first passage;
  the crossing row itself is excluded (FEM_CZM_HANDOFF.md §6: never assign
  a crossing event to an interval endpoint).
- `interpolate_target(trajectory, physical_time_s, channel="J"|"KJ")` —
  bounded linear interpolation; raises outside `[time[0], time[-1]]`
  (no extrapolation, matching the kernel's own no-extrapolation contract).

Test: `tests/test_pf_theta0_driving_force_trajectory_v10051840.py` — 7
tests. One (`test_real_pf_reference_prefracture_first_passage_matches_handoff`)
is `pytest.mark.skipif`-guarded on the real file's existence — now skips
due to the EXTERNAL BLOCKER above; the other 6 (synthetic-CSV) still pass.

Pure function of the CSV; touches no FEM/hazard/RNG/mesh state.

## sharp_front.py architecture (researched this session — expensive to
## re-derive; read this before touching sharp_front.py again)

`arrhenius_fracture/sharp_front.py::run_2d(args)` is the **single shared
step loop** for the entire `mode_i_first_passage_*` family. ~15 other
modules `inspect.getsource(sharp_front.run_2d)` and AST/text-patch it
rather than calling a different loop. **Treat `run_2d`'s source text as
load-bearing for other modules.** Before this session's edit, `git grep`
confirmed the only two external text anchors
(`kinetic_progressive_2d_v10.py`/`v1002.py`'s `adaptive_target = min(...)`
at line 2435, and an `info = eng.step(KJ, T, dt_cur)` anchor at ~3071) both
fall **outside** the edited region (~2451-2665 pre-edit numbering), so the
edit should not break those patchers.

Production entry chain for this workspace's target case:
`mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production.py:46`
→ `..._parity.py:298 main()` → `_base.main`
(`mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor`) → ... →
`sharp_front.py:1095 run_2d`. Every wrapper `main()` in this chain forwards
the raw CLI argv list mostly verbatim (only lightweight string-based
`_option_value`/`_require_*` checks, not real argparse) down to
`sharp_front.py`'s own `_build_parser().parse_args(argv)` (its `main()`,
line ~4454) — so a new flag added to `_build_parser()` flows through
cleanly as long as it's present in the CLI args list somewhere upstream.

Inside `run_2d` (theta=0, non-deflect, non-fatigue path — the
Peak/1000K/seed8666 case), pre-edit structure:
- Outer loop `while step < args.steps:`, inner trial loop `while True:`,
  which already snapshots `u_saved/ep_saved/rho_saved/Uapp_saved` before
  every trial and restores them on rejection.
- The fixed ramp (`dU_step = cfg.loading.dU_top * trial_frac`) was applied
  unconditionally, then the mechanics/plasticity solve ran, then KJ was
  computed via `compute_J_integral`, then the **only existing accept/reject**
  check (`predicted_clock > adaptive_target`, a hazard-clock/CFL safety
  check, not a KJ-target check) decided retry vs. accept.
- **Key finding**: front-engine hazard/RNG/MPZ state (`B`, `N_em`,
  `_hazard_rng`, event-length draw, etc.) is only touched **after** the
  trial loop breaks (step commit). It is not read or written during trial
  evaluation for the non-deflect path — so a KJ-target retry loop
  gets full transactionality for free from the existing
  `u_saved/ep_saved/rho_saved/Uapp_saved` scaffolding; no additional
  front-engine snapshot was needed for this integration.
- `compute_J_integral` (`j_integral.py:50`) and `solve_dirichlet`
  (`fem.py:173`) are pure/side-effect-free given their inputs — a trial
  evaluation (assemble → solve → update_plasticity → compute_J_integral)
  is exactly the `solve_trial(Uapp) -> KJ` contract Phase 2a's controller
  expects.
- Full front-engine transactional snapshot/restore already exists (for a
  different purpose) at `persistent_site_moving_tip_v100515.py:71-115`
  and subclasses — not needed for this integration, but the reference
  implementation if a future change moves hazard integration earlier into
  the trial loop.
- No pre-existing KJ/J-target retry loop existed anywhere; closest prior
  art (`mixed_mode_first_passage_v8.py`'s safeguarded secant for
  mode-mixity angle) wraps a whole `run_2d` call from outside, a different
  pattern than the per-step inner-loop retry needed here.

## Phase 2a — completed: standalone J/KJ-target controller (`3cc67bb`)

`arrhenius_fracture/pf_theta0_j_controlled_loading_v10051840.py`:
- `elastic_predictor(Uapp_old, achieved_old, target, power, achieved_min)`
  — `U_new = U_old * (target/max(achieved_old, achieved_min))**(1/power)`.
  `channel_power_for("J")==2.0` reduces to the handoff's literal sqrt
  formula; `channel_power_for("KJ")==1.0` gives the equivalent linear form.
- `solve_to_target(...)` — safeguarded predictor+secant loop: predicts,
  calls the caller's `solve_trial` callback, checks tolerance, falls back
  to a bounded secant correction (`max_growth_ratio`/`max_shrink_ratio`)
  capped at `max_iterations`. Never raises on non-convergence — returns
  `converged=False` with the best-effort last trial.
- Returns `JControlledStepResult` with full per-trial history
  (`TrialRecord`, `accepted` flag) for the event-resolved audit table
  FEM_CZM_HANDOFF.md §8 asks for.
- Deliberately touches no FEM/mesh/plasticity/MPZ/hazard/RNG state — pure
  scalar control logic, dependency-injected via `solve_trial`.

Test: `tests/test_pf_theta0_j_controlled_loading_v10051840.py` — 11/11
passing (predictor formula, exact/nonlinear convergence, safeguard
clamping, graceful non-convergence, achieved_min floor, rejected-trial
bookkeeping, config validation).

## Phase 2b — default-off integration: CHECKPOINTED, flag-ON qualification
## PENDING (do not describe as qualified)

### What changed (this checkpoint)

`arrhenius_fracture/sharp_front.py`:
1. New setup block before `while step < args.steps:`: if
   `args.pf_kj_target_csv` is set, loads the PF trajectory via
   `load_pf_theta0_driving_force_trajectory` + `prefracture_slice`, builds
   a `JControllerConfig`, initializes `physical_time_accepted = 0.0`,
   `KJ_accepted = 0.0`. If unset (default `None`), `pf_target_trajectory
   stays None` and nothing else changes.
2. Three new argparse flags appended at the very end of `_build_parser()`
   (right before its `return p`, so no existing argument-block text was
   touched): `--pf-kj-target-csv` (default `None`), `--pf-kj-target-rel-tol`
   (default 1e-3), `--pf-kj-target-max-iterations` (default 25).
3. Trial-loop restructure: `pf_target_active = pf_target_trajectory is not
   None and not deflect and not fatigue_mode`. When `pf_target_active`,
   looks up `KJ_target` at `physical_time_accepted + dt_cur`, defines a
   `_pf_solve_trial(Uapp)` closure resetting `u/ep_gp/rho_gp` from
   `u_saved/ep_saved/rho_saved` each call and re-solving (same
   assemble/solve/update_plasticity/`compute_J_integral` sequence as the
   original code), calls `solve_to_target(...)`, and sets
   `Uapp/KJ = result.accepted_Uapp/achieved` before falling through to the
   **unchanged** `predicted_clock = eng.predict_clock_increment(...)` call
   and the **completely unchanged** hazard-clock accept/reject block. When
   `pf_target_active` is False, the `else:` branch is the **original code,
   byte-identical, re-indented one level** — confirmed no `git diff`
   changes to its actual statements, only indentation.
   - First-step-from-rest handling: when `Uapp_saved <= 0` (step 0 only),
     seeds the controller with one real solve at
     `max(cfg.loading.dU_top, 1e-12)` rather than a degenerate `(0, 0)`
     pair (the predictor is undefined at `Uapp_old=0`, per
     `test_achieved_min_floor_prevents_division_by_zero_from_rest_state`).
     Later steps seed from `(Uapp_accepted, KJ_accepted)` with no extra solve.
4. Post-accept commit: added `if pf_target_trajectory is not None:
   KJ_accepted = float(KJ); physical_time_accepted += dt_cur` right after
   the existing `Uapp_accepted = Uapp` — additive, no-op when unset.

`run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`:
- Added `PF_KJ_TARGET_CSV=${PF_KJ_TARGET_CSV:-}` (default empty) and a
  `PF_KJ_TARGET_FLAG` computed the same way as the pre-existing
  `PLOT_FLAG` pattern, appended to the production-entry invocation via
  `${PF_KJ_TARGET_FLAG:+"$PF_KJ_TARGET_FLAG" "$PF_KJ_TARGET_CSV"}`.
  Verified this exact two-word optional-flag Bash idiom works correctly
  under the system's Bash 3.2 (both empty and non-empty cases tested
  standalone) before using it, given this script's own prior history of a
  Bash-3.2-specific bug (`2934911 Make optional plot flag compatible with
  macOS Bash 3.2`). Did **not** touch the pre-existing `FAMILY_SHA_REQUIRED`
  kernel-hash guard or its comparison logic anywhere in this diff.

### Verification evidence for this checkpoint

- `py_compile` on `sharp_front.py`: succeeds.
- `bash -n` on the launcher script: succeeds.
- Focused suite: **36 passed, 1 skipped** (the 1 skip is the
  `skipif`-guarded real-PF-reference integration test, correctly skipping
  now that the external file is gone — not a new failure). `git diff
  --check`: clean.
- Default-off (flag absent) bit-identical evidence: an 800-step smoke run
  performed **immediately after the `sharp_front.py` edit, before the
  launcher-script edit**, reproduced the exact same KJ values as the Phase
  0 baseline at every printed step (200/400/600/800 → 2.215/3.111/3.786/
  4.346 MPa√m) — proves the `sharp_front.py` restructuring alone is a
  behavioral no-op when `--pf-kj-target-csv` is absent.
- A **third** attempted re-verification (after also editing the launcher
  script) hit the EXTERNAL BLOCKER above — the launcher's own
  pre-existing precondition check (present before any of this session's
  edits, used later in the same script for the `compare_v10051840_theta0_pf_parity.py`
  call) requires `steps_1000K.csv` to exist at the documented path
  **regardless of the new flag**, and that file is now gone externally.
  This is not a regression from this checkpoint's edits — it is the same
  external blocker, confirmed by inspecting the failing script's own
  precondition-check output.
- Given the above, the launcher-script diff was instead verified by (a)
  direct code review (the new lines are purely additive and default to
  producing zero extra CLI tokens when `PF_KJ_TARGET_CSV` is unset — see
  diff), and (b) a standalone Bash 3.2 test of the exact
  `${VAR:+"$A" "$B"}` expansion pattern confirming it contributes zero
  words when empty and exactly two correctly-quoted words when set. The
  composition of (evidence for `sharp_front.py` alone) + (proof the
  launcher change is a no-op when unset) covers the "flag-off remains
  bit-identical" requirement even though a fresh end-to-end run via the
  launcher could not be performed a third time this session.
- `git grep` confirms: no hardcoded fallback/replacement PF path exists
  anywhere in the new code (the CSV path is only ever taken from
  `args.pf_kj_target_csv` / `$PF_KJ_TARGET_CSV`, both user/launcher
  supplied); the `FAMILY_SHA_REQUIRED` guard lines are untouched by the
  diff; `runs/` remains gitignored (new run artifacts were not staged).

### What is explicitly NOT yet verified (do not claim otherwise)

- No flag-ON run against real PF reference data has completed. The one
  attempt failed at the launcher's own precondition check due to the
  EXTERNAL BLOCKER, before any controller code executed.
- No evidence yet that the controller converges within
  `max_iterations`, tracks a real PF KJ(t) trajectory sanely, or reaches
  first passage faster than the fixed ramp, against real (non-synthetic)
  FEM state.
- UPDATE (`99f6e5c`): a real-integration test now exists
  (`tests/test_pf_theta0_j_controlled_loading_real_engine_v10051840.py`,
  4 tests) — but against a synthetic target trajectory and the minimal
  `legacy_scalar` front engine, not the real PF reference data or
  production front-engine configuration. See
  `FEM_PF_PARITY_SCORECARD.md` Gate 1 row for exactly what this does and
  does not prove.

## Files changed this session (cumulative)

- New, committed at `7562911`: `arrhenius_fracture/pf_theta0_driving_force_trajectory_v10051840.py`,
  `tests/test_pf_theta0_driving_force_trajectory_v10051840.py`
- New, committed at `3cc67bb`: `arrhenius_fracture/pf_theta0_j_controlled_loading_v10051840.py`,
  `tests/test_pf_theta0_j_controlled_loading_v10051840.py`
- Modified, committed at `ad2c446`: `arrhenius_fracture/sharp_front.py`,
  `run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`,
  `CLAUDE_PROGRESS.md`
- New, committed at `040610f`: `FEM_PF_PARITY_SCORECARD.md`,
  `PF_REFERENCE_REGENERATION_CONTRACT.md`
- Modified/new, committed at `99f6e5c`: `arrhenius_fracture/sharp_front.py`
  (controller audit-output block), `tests/test_pf_theta0_j_controlled_loading_real_engine_v10051840.py`
- Untracked, not committed (correctly gitignored under `runs/`):
  various `runs/*_20260804/` smoke/regression-check directories, several
  of which failed at the launcher's PF-file precondition check once the
  EXTERNAL BLOCKER hit (documented, not a code defect)
- Untracked, not part of this project's code, not touched:
  `.claude/settings.json`, `.claude/settings.local.json` (local harness
  config, not FEM/CZM source)

## Tests run (cumulative)

- 19/19 Phase 0 focused tests.
- 7/7 Phase 1 trajectory-reader tests (6 synthetic + 1 skipif-guarded real).
- 11/11 Phase 2a controller tests.
- 4/4 real-engine integration tests (`99f6e5c`) — real mesh/solve/plasticity
  against a synthetic target, no PF artifact needed.
- Latest combined focused-suite run: 40 passed, 1 skipped (the
  skipif-guarded real-PF-reference test correctly skips per the EXTERNAL
  BLOCKER — expected, not a failure).

## Latest accepted physical state

No FEM/CZM physics has changed. The default (flag-off) production path is
byte-identical in behavior to before this session (proven via bit-identical
smoke reproduction). A new, entirely opt-in code path exists and is
believed correct by construction (code review + the sharp_front.py
architecture notes above + Phase 2a's unit tests) but has not yet been
exercised against real FEM state end-to-end.

## Latest numerical blocker

Two, layered:
1. **Physics/controller** (unchanged from before): the fixed
   remote-displacement ramp produces a far slower J/KJ trajectory than
   the PF reference (~145,000 steps extrapolated to first passage). The
   fix (the Phase 2b controller) is now wired in behind a flag but not
   yet physically validated — see item 2.
2. **External/operational** (new this session): the immutable PF
   reference artifacts this fix needs to validate against are not
   currently trustworthy at their documented path (see EXTERNAL BLOCKER).
   This blocks flag-ON validation and all of Phase 3 until resolved.

## Failed approaches that should not be repeated

- Do not "fix" the slow ramp by tuning displacement rate directly —
  CLAUDE.md forbids tuning around a defective loading controller.
- Do not run `scripts/preflight_and_launch_v10051840_theta0_full_field.sh`
  on this dev branch — it hard-fails on a branch-name check by design.
- Do not use heredoc-piped `conda run -n ... python3 - <<'PY'` expecting
  captured stdout when redirected to a file — silently empty output was
  observed. Write the script to a real file and run
  `conda run -n ... python3 /path/to/script.py` instead.
- Do not retry the flag-ON smoke against the live external PF repo path
  hoping it stabilizes — see EXTERNAL BLOCKER.

## Exact next command

1. **Do not push yet.** After this checkpoint commit lands, the exact push
   command (to run only when explicitly approved) is:
   ```
   git push -u origin claude/v10.0.5.18.4.0-j-controlled-loading
   ```
   Do not modify PR #64, merge, or change its draft status.
2. Obtain or regenerate the authoritative PF `steps_1000K.csv` and kernel
   `family.json`, copy them into a dedicated immutable local reference
   directory inside this workspace (outside the live, externally-mutable
   `PF-fracture-fatigue_.../runs/` tree), and independently re-verify both
   SHA-256 values recorded in the EXTERNAL BLOCKER section above before
   trusting them again.
3. Only then: run the flag-ON smoke
   (`PF_KJ_TARGET_CSV=<frozen local copy>/steps_1000K.csv`, short `STEPS`,
   e.g. 50), inspect that KJ tracks the PF target sanely with a reasonable
   iteration count and no exceptions, add a dedicated regression test for
   the real-integration path, then commit as
   `feat: qualify PF KJ-target controller against frozen reference data`
   (or similar) — only at that point is it fair to describe the
   controller as physically qualified.
4. Proceed to Phase 3 (prefracture parity) only after step 3 passes.

## Exact next scientific acceptance gate

Phase 3 prefracture parity for Peak / 1000 K / theta=0 / seed=8666: FEM
reaches PF's first-passage J≈6531 J/m² (KJ≈53.9 MPa√m) in a practical
number of controller-driven steps (not ~145,000), with B(t), N_em(t), and
bulk plastic work tracking the PF trajectory shape, before any
crack-growth qualification (Phase 4) is attempted. Blocked until the
EXTERNAL BLOCKER above is resolved.

## Continuation handoff (for a fresh session or Codex)

Everything needed to resume is in this file plus the four new Phase
1/2a/2b files. A fresh agent should: (1) read this file, CLAUDE.md, and
FEM_CZM_HANDOFF.md in full; (2) re-run `git log --oneline -5` and
`git status --short`/`git diff --stat` to confirm no uncommitted drift
from what's recorded here; (3) re-run the focused suite (`36 passed, 1
skipped` is the expected current baseline) to re-confirm before writing
new code; (4) check whether the EXTERNAL BLOCKER has been resolved — if
not, do not attempt flag-ON work; if resolved, follow "Exact next command"
above starting at step 2. Phase 0, Phase 1, Phase 2a, and the default-off
Phase 2b wiring do not need to be repeated or re-derived — in particular,
do not re-run the `sharp_front.py` architecture research (Explore-agent
pass); it is fully recorded above.
