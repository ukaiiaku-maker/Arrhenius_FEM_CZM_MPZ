# Claude progress

- Updated: 2026-08-03 22:40 PDT
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
trajectory reader) is **complete**. Phase 2a (the controller's decision
logic: predictor + safeguarded secant + tolerance/iteration bookkeeping,
as a standalone FEM-independent module) is **complete**. Phase 2b (wiring
that controller into the real `sharp_front.run_2d` step loop) is **not
started** — this is the next task and is the highest-risk remaining step;
see "sharp_front.py architecture" and "Exact next command" below before
touching that file.

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

## sharp_front.py architecture (researched this session — expensive to
## re-derive; read this before touching sharp_front.py)

`arrhenius_fracture/sharp_front.py::run_2d(args)` is the **single shared
step loop** for the entire `mode_i_first_passage_*` family. ~15 other
modules (`kinetic_progressive_2d_v10/v1002/v1003.py`,
`mode_i_first_passage_v9_18_5.py`,
`mode_i_first_passage_v10_0_5_3_fatigue_audited.py`,
`mode_i_first_passage_v10_0_5_5_stochastic_vhcf.py`, etc.) do not call a
different loop — they `inspect.getsource(sharp_front.run_2d)` and
AST/text-patch it, or monkey-patch pieces of it, to build progressively
specialized variants. **Treat `run_2d`'s source text as load-bearing for
other modules, not just for its own behavior** — an edit that changes
matched text patterns could silently break unrelated wrapper modules.

The production entry chain for this workspace's target case is:
`mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production.py:46`
→ `..._parity.py:298 main()` (monkey-patches policy/defaults) →
`_base.main` (`mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor`) →
... → `sharp_front.py:1095 run_2d`.

Inside `run_2d`, the per-step structure (theta=0, non-deflect, non-fatigue
path, which is what the Peak/1000K/seed8666 case uses):
- Outer loop: `while step < args.steps:` at **`sharp_front.py:2451`**.
- Inner trial loop: `while True:` at **`sharp_front.py:2456`**.
  - **Already snapshots** `u_saved = u.copy(); ep_saved = ep_gp.copy();
    rho_saved = rho_gp.copy(); Uapp_saved = Uapp_accepted` at
    **lines 2457-2460**.
  - Fixed ramp applied at **line 2472**: `dU_step = cfg.loading.dU_top *
    trial_frac` (production forces `--dU 2.0e-7`), then `Uapp = Uapp_saved
    + dU_step` (line 2473).
  - FEM solve: `assemble_mechanics`/`solve_dirichlet`/`update_plasticity`
    inside an `n_stagger` loop at **lines 2477-2492**. `update_plasticity`
    DOES mutate `ep_gp`/`rho_gp` even during a "trial" — but this is safe
    because of the snapshot above.
  - KJ computed (non-deflect path) at **lines 2554-2557**:
    `J, KJ, _ = compute_J_integral(mesh, u, sigma_gp, psi_gp, d,
    np.array([a_tip, 0.0]), np.array([1.0, 0.0]), mat, ell=..., ...)`.
  - **Only existing accept/reject** at **lines 2581-2586**: rejects only if
    `predicted_clock > adaptive_target` (a hazard-clock/CFL safety check,
    NOT a KJ-target check); on reject it restores `u, ep_gp, rho_gp, Uapp`
    from the saved copies and shrinks `trial_frac`; on accept it `break`s.
  - **Key finding**: front-engine hazard/RNG/MPZ state (`B`, `N_em`,
    `_hazard_rng`, `_emission_rngs`, event-length draw, etc.) is only
    touched **after** this trial loop breaks (step commit, further down in
    `run_2d`, not yet located precisely). It is NOT read or written during
    trial evaluation for the non-deflect path. This means: **a KJ-target
    retry loop that replaces/extends the lines-2581-2586 accept condition
    gets full transactionality for free from the existing
    u_saved/ep_saved/rho_saved/Uapp_saved scaffolding** — no additional
    front-engine snapshot is needed for the retry loop itself, because
    front-engine state simply isn't mutated yet at that point. (A full
    front-engine `_capture_state()`/`_restore_state()` — see below — would
    only become necessary if a future change moves hazard integration
    earlier into the trial loop; it is not necessary for the currently
    planned integration.)

Full front-engine transactional snapshot/restore already exists (for a
different purpose — internal microstepping retries — but is the reference
implementation if ever needed): `_capture_state()`/`_restore_state()`
rooted at `persistent_site_moving_tip_v100515.py:71-115`, subclassed at
`persistent_site_stochastic_tip_v100516.py:188` (adds hazard RNG
`bit_generator.state`, threshold/action/event-index, event-length draw),
`persistent_site_stochastic_emission_v100518.py:331` (adds per-slip-system
emission RNGs + thresholds), etc. Do **not** use
`arrhenius_fracture/cohesive_trial_state.py::KineticCZMTransactionSnapshot`
for this — it is gated on `state_model == 'kinetic_campaign_czm'` and
`ARRHENIUS_CZM_OPENING_COUPLING == 'clock_linear'`, neither of which the
theta0 production entry sets (it uses `front_state_model = 'moving_pz'`
and the `atomic_path_corridor_czm` backend) — that transaction class is
simply not wired into this production path.

`compute_J_integral` (`j_integral.py:50`) is pure/read-only given
`(mesh, u, sigma_gp, psi_gp, d, ...)`. `solve_dirichlet` (`fem.py:173`)
returns a new array and doesn't mutate inputs. So a trial evaluation
(assemble → solve → optionally update_plasticity → compute_J_integral) is
exactly what `solve_trial(Uapp) -> KJ` in the new controller module
expects, matching the existing trial-loop body almost verbatim.

There is no pre-existing KJ/J-target retry loop anywhere in the codebase.
The closest prior art is `mixed_mode_first_passage_v8.py`'s safeguarded
secant controller for mode-mixity phase angle (`safeguarded_alpha_update`,
line ~89), but that wraps an entire `sf.run_2d(args)` call from outside
via monkey-patched `femmod.solve_dirichlet`/`jimod.compute_J_integral`,
not a per-step inner-loop retry — architecturally a different pattern
than what's needed here (we need to seek a target within one step, not
across a whole run).

## Phase 2a — completed: standalone J/KJ-target controller decision logic

New module: `arrhenius_fracture/pf_theta0_j_controlled_loading_v10051840.py`
- `elastic_predictor(Uapp_old, achieved_old, target, power, achieved_min)` —
  `U_new = U_old * (target/max(achieved_old, achieved_min))**(1/power)`.
  `channel_power_for("J")==2.0` reduces this to the handoff's literal
  square-root formula; `channel_power_for("KJ")==1.0` gives the equivalent
  linear form (since KJ ~ sqrt(J) ~ U in the elastic prefracture regime).
- `solve_to_target(Uapp_initial, achieved_initial, target, solve_trial,
  config)` — the safeguarded predictor+secant loop: predicts, calls the
  caller's `solve_trial` callback, checks tolerance, and on rejection
  falls back to a secant correction (bounded by `max_growth_ratio`/
  `max_shrink_ratio` per iteration — the "bounded secant/proportional
  correction" the handoff requires), capped at `max_iterations`. Never
  raises on non-convergence — returns `converged=False` with the
  best-effort last trial so the caller can decide (e.g. right-censor the
  step, matching the `control_state: right_censored_endpoint` convention
  already used elsewhere in this codebase's run outputs).
- Returns `JControlledStepResult` with the full `trials` history
  (`TrialRecord` per iteration, `accepted` flag) so a caller/audit can
  record every rejected trial's `(Uapp, achieved)` pair — feeding directly
  into the event-resolved audit table FEM_CZM_HANDOFF.md section 8 asks
  for (`target J`, `achieved J`, `controller iteration count`, etc.).
- **Deliberately does not touch FEM state, mesh, plasticity, MPZ, hazard,
  or RNG** — it is pure scalar control logic, dependency-injected via
  `solve_trial`. This was a deliberate scope decision this session: wiring
  it into the monolithic, AST-patched `run_2d` is materially riskier and
  deserves its own dedicated step with a real production smoke test, not
  a rushed edit bundled with writing the algorithm itself.

Regression test: `tests/test_pf_theta0_j_controlled_loading_v10051840.py`
— 11/11 passing. Covers: predictor formula reduction to the handoff's
literal sqrt form (J channel) and its linear KJ-channel equivalent; exact
1-iteration convergence for a synthetic linear response; secant-refined
convergence for a synthetic nonlinear (softening) response; safeguard
clamping of an oversized predictor jump; graceful (non-exception)
non-convergence reporting under a pathological constant-response
callback; zero-call short-circuit when already within tolerance; the
achieved_min floor preventing division-by-zero from a Uapp=0/achieved=0
rest state; rejected-trials bookkeeping; and config validation.

## Files changed this session

- New: `arrhenius_fracture/pf_theta0_driving_force_trajectory_v10051840.py`
  (Phase 1, committed at 7562911)
- New: `tests/test_pf_theta0_driving_force_trajectory_v10051840.py`
  (Phase 1, committed at 7562911)
- New: `arrhenius_fracture/pf_theta0_j_controlled_loading_v10051840.py`
  (Phase 2a, not yet committed as of this update)
- New: `tests/test_pf_theta0_j_controlled_loading_v10051840.py`
  (Phase 2a, not yet committed as of this update)
- Modified: `CLAUDE_PROGRESS.md` (this file)
- Untracked, not committed: `runs/phase0_slow_ramp_reproduction_smoke_20260803/`
  (local smoke-run output; `runs/` stays out of git per CLAUDE.md)

## Commits created

- `7562911` — `feat: add PF driving-force trajectory reader` (Phase 1).
- Pending this update: a second commit for the two Phase 2a files, planned
  message `feat: add safeguarded J/KJ-target loading controller`. Check
  `git log --oneline -3` for the actual resulting hash after this update
  is committed.

## Tests run

- 19/19 Phase 0 focused tests (see list above).
- 7/7 Phase 1 trajectory-reader tests.
- 11/11 Phase 2a controller tests.
- 37/37 combined (Phase 0 set + both new test files) run together — no
  interaction/import-order issues.

## Simulations run

- Read-only inspection of the archived slow-ramp baseline (no new run).
- One completed 800-step local smoke run (see Phase 0 section) — confirmed
  bit-identical KJ trajectory vs. the archived baseline; not run to first
  passage (would take ~145,000 steps at the current fixed-displacement
  ramp rate, which is precisely the defect Phase 2 will fix).

## Latest accepted physical state

**`sharp_front.py` (the real FEM/CZM solver loop) has still not been
modified.** Both new modules this session (Phase 1 reader, Phase 2a
controller logic) are pure, read-only/side-effect-free additions that sit
beside the solver, not inside it. The archived + freshly reproduced
slow-ramp evidence (step 9400, KJ=11.4 MPa√m, B=0, N_em=0, no growth) is
still the accepted Phase 0 baseline. The PF prefracture target trajectory
(176 rows, up to KJ=53.9 MPa√m) and a tested, FEM-independent
predictor+secant target-seeking algorithm are now both available in-code,
ready to be wired into `run_2d`.

## Latest numerical blocker

The fixed remote-displacement ramp (dU/dt matched to the PF rate) produces
a far slower J/KJ trajectory than the PF reference (~145,000 steps
extrapolated to first passage). This is a controller problem, not a
hazard/kernel/mesh problem — none of B, N_em, or crack extension have
moved yet at step 9400. Root cause per FEM_CZM_HANDOFF.md §5: FEM and PF
geometries have different compliance, so equal dU/dt does not produce
equal dJ/dt. **Still not fixed** — the fix requires editing
`sharp_front.py::run_2d`'s trial loop (lines ~2456-2586, see architecture
notes above), which has not happened yet this session. Only the two
supporting modules it will call (trajectory reader, controller) have been
built and tested in isolation.

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

Phase 2b: wire `pf_theta0_j_controlled_loading_v10051840.solve_to_target`
into `sharp_front.py::run_2d`'s trial loop. This is the highest-risk step
remaining (editing a monolithic function that ~15 other modules
AST/text-patch — see "sharp_front.py architecture" above) — read that
section in full before starting. Concretely:

1. Do **not** edit `run_2d` in place as the first move. First write a
   small standalone script/test that imports `sharp_front`, builds the
   real args for the Peak/1000K/theta0/seed8666 case (same as
   `run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`'s CLI
   flags), and confirms you can call `assemble_mechanics` +
   `solve_dirichlet` + `compute_J_integral` directly (the same calls at
   `sharp_front.py:2477-2492,2554-2557`) to get a `solve_trial(Uapp)->KJ`
   closure working in isolation against the real mesh/materials, before
   touching the loop itself.
2. In `run_2d`, locate the exact accept/reject block at
   **lines 2581-2586** (`if adaptive_events and predicted_clock >
   adaptive_target ...`). Add a **new, separately flagged** code path
   (e.g. gated on a new `args.j_target_controlled` / an env var, default
   OFF) that, instead of shrinking `trial_frac` to satisfy the hazard
   clock, calls `solve_to_target` with `solve_trial` wrapping the existing
   assemble/solve/(optionally update_plasticity)/compute_J_integral block
   at a candidate `Uapp`, targeting `interpolate_target(prefracture_slice,
   physical_time_s, channel="KJ")` for the current step's PF-reference
   physical time. Reuse the loop's own `u_saved/ep_saved/rho_saved/
   Uapp_saved` restore-on-reject pattern for every rejected controller
   trial (confirmed transactional for free — see architecture notes).
   Do NOT call `update_plasticity` on rejected trials if avoidable, or if
   it must be called (n_stagger loop structure), rely on the existing
   `ep_saved`/`rho_saved` restore to undo it exactly as the current
   hazard-clock path already does.
3. Keep the existing `--dU`-fixed-ramp default path completely unchanged
   (default OFF for the new flag) so the Phase 0 slow-ramp reproduction
   stays re-runnable as a regression baseline for comparison.
4. Add a focused regression proving: (a) `ep_gp`/`rho_gp`/`u`/`Uapp` are
   bit-identical before and after a rejected controller trial in a real
   (small) FEM setup, and (b) a short controller-driven run reaches a
   target KJ materially faster (fewer accepted steps) than the fixed-ramp
   baseline for the same physical-time interval.
5. Run the Phase 0+1+2a 37-test focused suite plus the new regression,
   then a real production-entry smoke (short STEPS, as in the Phase 0
   smoke command) with the new flag enabled, comparing against the
   archived/reproduced slow-ramp baseline's early KJ values as a sanity
   check that non-controller behavior is unaffected when the flag is off.
6. Commit as `feat: add transactional J-controlled loading`, matching
   CLAUDE.md's suggested milestone name.

## Exact next scientific acceptance gate

Phase 3 prefracture parity for Peak / 1000 K / theta=0 / seed=8666: FEM
reaches PF's first-passage J≈6531 J/m² (KJ≈53.9 MPa√m) in a practical
number of controller-driven steps (not ~145,000), with B(t), N_em(t), and
bulk plastic work tracking the PF trajectory shape, before any
crack-growth qualification (Phase 4) is attempted.

## Continuation handoff (for a fresh session or Codex)

If context is exhausted before Phase 2b lands, everything needed to
resume is in this file plus the four new Phase 1/2a files — re-deriving
the "sharp_front.py architecture" section above cost a full Explore-agent
pass this session; do not repeat that research, just read it here. A
fresh agent should: (1) read this file, CLAUDE.md, and
FEM_CZM_HANDOFF.md in full; (2) re-run `git log --oneline -5` and
`git status --short` to confirm no uncommitted drift; (3) re-run the
37-test combined focused suite listed above to re-confirm the accepted
baseline before writing new code; (4) proceed directly to the "Exact next
command" section above (Phase 2b: wiring into `sharp_front.run_2d`) —
Phase 0, Phase 1, and Phase 2a do not need to be repeated.
