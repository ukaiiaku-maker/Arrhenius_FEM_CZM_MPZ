# FEM/CZM ↔ PF Parity Scorecard

This file, not the current coding task, organizes the work. Every session
should update it before ending. See CLAUDE.md and FEM_CZM_HANDOFF.md for
the governing contract, and CLAUDE_PROGRESS.md for narrative session-level
detail (commits, exact commands, architecture notes). This file is the
scientific status ledger.

Reference case for all rows below unless stated otherwise: Peak
parameterization (`v913_paper_peak01_0242980_persistent_sites`), 1000 K,
theta=0, hazard seed 8666, PF run
`v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1`.

## Status summary (2026-08-04, updated)

**Controller infrastructure complete and now validated end-to-end against
the real FEM engine with a synthetic target; physical FEM/PF parity still
untested.**

- The PF driving-force trajectory reader and the safeguarded J/KJ-target
  controller exist, are unit-tested in isolation (18 synthetic tests
  total across both modules), and are wired into `sharp_front.run_2d`
  behind a default-off flag. The default (fixed-ramp) path is proven
  bit-identical to its pre-existing behavior.
- **New this update**: a real-engine integration test (`sharp_front.main`
  called directly, real mesh/assembly/solve/plasticity, no PF kernel
  needed) drove the controller against a small synthetic target KJ(t)
  trajectory and confirmed, against the real solver: achieved KJ tracks
  target within tolerance at every step (2-5 iterations, well under the
  25-iteration cap); physical time is inherited from the target CSV's
  `dt` column, not solver-step count; rejected trials are recorded and
  precede the accepted value; repeated identical runs are bit-identical
  (no cross-trial or cross-run state leakage). This is genuine Gate 1
  evidence — **but against a synthetic target and the minimal
  `legacy_scalar` front engine, not the real PF-audited configuration**,
  so Gate 1 remains "partial," not "pass," until the same is demonstrated
  against the real PF reference trajectory and the production front-engine
  configuration (persistent-site MPZ, signed kernel, full-field bulk PT).
- **No real flag-ON run against the audited PF reference data has
  completed.** Those artifacts (`steps_1000K.csv`, kernel `family.json`)
  are currently unavailable at their documented path — see
  CLAUDE_PROGRESS.md's "EXTERNAL BLOCKER" section and
  `PF_REFERENCE_REGENERATION_CONTRACT.md`.
- **Therefore every row below Gate 1 is still "not tested" or "blocked."**
  No claim of FEM/PF physical parity has been demonstrated at Gate 2 or
  above.
- Important reminder from the governing instructions: matching a
  prescribed PF `KJ(t)` does **not** by itself prove parity. It only
  removes loading-rate/boundary-condition differences so that Gate 2
  (internal-state comparison at matched driving force) becomes a fair
  test. Do not stop at "the controller runs."

## Gate ladder status

| Gate | Description | Status | Notes |
|---|---|---|---|
| 1 | Controlled prefracture mechanics (controller numerics) | **partial** | Predictor/secant logic validated by 11 synthetic-callback unit tests. Now ALSO validated by 4 real-engine integration tests (`99f6e5c`): real mesh/solve/plasticity, synthetic target, achieved KJ tracks target within tolerance, practical iteration counts (2-5 of 25), physical time correctly inherited from target `dt`, bit-identical across repeats. Remaining gap to "pass": same demonstration against the real PF reference trajectory and the production front-engine configuration (persistent-site MPZ + signed kernel + full-field bulk PT), not yet possible — see EXTERNAL BLOCKER. |
| 2 | Prefracture internal-state parity (B, N_em, plastic work, MPZ) | **blocked** | Requires a completed Gate-1 real run against real PF reference data. |
| 3 | Stochastic first-passage parity | **blocked** | Same. |
| 4 | First-event FEM/CZM transaction audit | **blocked** | Same; also requires the controller audit-output fields (FEM_CZM_HANDOFF.md §8) to be implemented — not yet done. |
| 5 | Short-growth parity (20-50 µm) | **not tested** | Not reached. |
| 6 | Intermediate/long growth (400 µm, 1000 µm) | **not tested** | Not reached. |
| 7 | Four immutable material classes | **not tested** | Not reached. |
| 8 | Temperature dependence | **not tested** | Not reached. |

## Parity quantity table

| Quantity | PF reference | FEM/CZM result | Metric | Tolerance | Status | Dominant discrepancy | Likely subsystem | Supporting run / commit | Next experiment |
|---|---|---|---|---|---|---|---|---|---|
| Driving force J(t), KJ(t) | `steps_1000K.csv` (176 prefracture rows before first passage; captured this session, live copy currently unavailable — see blocker) | Fixed-ramp baseline reaches only KJ=11.37 MPa√m by step 9400 (extrapolated ~145,000 steps to first passage). Against a **synthetic** target with the real engine: achieved KJ matches target within ~2e-3 relative at every step (`99f6e5c`). Against the **real PF target**: none yet. | RMSE / relative error of achieved vs. target KJ(t) | ≤ controller `rel_tol`=1e-3 by construction, confirmed achievable against synthetic target | **partial** (synthetic target only) / **blocked** (real PF target) | Real PF target not yet run | loading/controller | Phase 0 archived baseline; `ad2c446` (wiring); `99f6e5c` (real-engine synthetic-target evidence) | Recover frozen reference (see regeneration contract), then repeat this exact test against the real PF row-176 prefracture trajectory |
| Bulk plasticity Wp(t), ep, rho, active volume | Not yet extracted from PF steps table (columns exist: `W_bulk_plastic_cumulative_J_per_m` etc.) | Not tested | Time series comparison at matched physical time | not yet defined | not tested | — | bulk plasticity / Peierls-Taylor kinetics | none | After Gate 1: extract PF `W_bulk_plastic_cumulative_J_per_m` and FEM equivalent at matched t |
| Moving-tip MPZ populations/translation | Not yet extracted | Not tested | Population time series | not yet defined | not tested | — | moving process zone | none | Same |
| Emission N_em(t), hazard/action | PF `N_em` column exists in steps table | Not tested | Time series | not yet defined | not tested | — | emission hazard | none | Same |
| Cleavage B(t), threshold, accumulated action | PF `B` column exists in steps table | Not tested | Time series | not yet defined | not tested | — | cleavage hazard | none | Same |
| First passage: time, J, KJ, threshold identity | J≈6531 J/m², KJ≈53.9019 MPa√m at PF row 176 (step 177, t=1486.8s) — captured exactly this session | Not tested (fixed ramp never reaches it in a practical step count; controller not yet run) | Relative error in J_FP, KJ_FP; first-passage time comparison | ε_J,FP tolerance not yet defined (handoff suggests "comparable", no numeric bound given yet) | **blocked** | n/a | loading/controller, then cleavage hazard | Phase 0 (fixed-ramp extrapolation); Phase 1 (`7562911`, exact PF row located) | Gate 3 once Gate 1/2 pass |
| Event proposal: length factor, physical length, endpoint | PF event-length draw mechanism (threshold_scaled, factor range 0.5-4.0 per production env vars) | Not tested | Exact identity where seeded identically; distribution otherwise | not yet defined | not tested | — | event-length law | none | Gate 4 first-event audit |
| Energy gate: available vs. required | Not yet extracted | Not tested | — | not yet defined | not tested | — | energy gate | none | Gate 4 |
| Cohesive transaction: trial/commit/veto/rollback identity | n/a (PF has no cohesive geometry; sharp-wake backend) | FEM atomic path-corridor CZM backend exists and is exercised by existing focused tests (unrelated to the new controller) | Existing gate-quality regression only | existing v3.9 tolerances (triangle quality ≥0.035, area ratio ≥0.08) | not tested (for controller-driven steps specifically) | — | cohesive representation / geometry-remeshing | `test_atomic_path_corridor_v10051839.py` (pre-existing, still passing) | Gate 4 audit under controller-driven loading specifically |
| Geometry: tip position, path, mesh quality | n/a | Existing mesh-quality gates pass in isolation (pre-existing tests) | — | existing | not tested (under controller) | — | geometry/remeshing | pre-existing tests | Gate 4/5 |
| Renewal: new threshold, B reset, MPZ update | Not yet extracted | Not tested | — | not yet defined | not tested | — | cleavage hazard / MPZ | none | Gate 4 |
| Post-event continuation | n/a | Not tested | Does the next accepted step succeed after a commit? | pass/fail | not tested | — | state transfer / restart | none | Gate 4 |
| Short growth (20-50 µm) event history | PF event statistics not yet extracted | Not tested | Event count, length distribution, spacing, J/KJ at events | not yet defined | not tested | — | multiple (see gate 5 checklist) | none | After Gate 4 passes |
| Intermediate growth (400 µm) | Not yet extracted | Not tested | Developed statistics (distributions, not per-event identity) | not yet defined | not tested | — | multiple | none | After Gate 5 |
| Long growth (1000 µm) | Not yet extracted | Not tested | Same | not yet defined | not tested | — | multiple | none | After Gate 6a |
| Four-class response (Peak/DBTT/weakT/ceramic) | PF ordering/characteristic behavior not yet extracted for classes other than Peak | Not tested | Ordering + characteristic shape, not per-value fit | not yet defined | not tested | — | material-row transfer / common numerics | none | After Gate 6 |
| Temperature response (transition shape/ordering) | PF temperature sweep not yet extracted | Not tested | Transition location/shape, ordering | not yet defined | not tested | — | multiple | none | After Gate 7 |

## Session reporting log

Newest entry first. Use the required form from the governing instructions.

### 2026-08-04 (continued — real-engine controller audit)

```text
Highest completed gate:      None (Gate 1 remains partial).
Current gate:                Gate 1 (controlled prefracture mechanics).
PF reference:                Not used this step -- a synthetic 10-row
                              target trajectory (KJ linear 0.1-5.0 MPa√m,
                              dt=8.4s/row) substituted because the real PF
                              artifacts remain unavailable (EXTERNAL
                              BLOCKER). This was a deliberate substitution
                              to test controller MECHANICS against the
                              real solver, not a substitution for physical
                              parity evidence.
FEM/CZM result:               Real sharp_front.main() run (legacy_scalar
                              front engine, 6x10 mesh, no kernel), 3 steps:
                              achieved KJ = 0.09999518, 0.64443629,
                              1.18783401 MPa√m vs targets 0.1, 0.644444,
                              1.188889 MPa√m -- all within ~2e-3 relative.
                              Iterations: 5, 3, 2 (cap 25). Physical time
                              8.4/16.8/25.2 s, matching cumulative dt
                              exactly (not solver-step count).
Agreement:                    Achieved vs. synthetic target: agrees within
                              ~5e-5 to 2e-3 relative at each step (well
                              inside the controller's own rel_tol=1e-3
                              acceptance criterion, confirmed operating
                              correctly against the real engine).
Largest discrepancy:          None at this synthetic-target scale -- this
                              step tests controller/solver mechanics, not
                              FEM/PF physics, so "discrepancy" is not yet
                              a meaningful question here.
Likely subsystem:             loading/controller (mechanics only; no
                              physics comparison performed).
Evidence:                     4 new integration tests in
                              tests/test_pf_theta0_j_controlled_loading_real_engine_v10051840.py:
                              tolerance tracking, rejected-trial recording
                              and precedence, cross-run determinism
                              (bit-identical), and flag-off no-op (no
                              audit file written).
Change made:                  Added a controller audit-output block to
                              sharp_front.run_2d (pf_controller_audit list,
                              written to
                              pf_kj_target_controller_audit_<T>.json when
                              the flag is set) recording per-step target
                              KJ, achieved KJ, physical time, converged
                              flag, iteration count, and every rejected
                              trial's (Uapp, achieved) pair -- the
                              event-resolved audit trail FEM_CZM_HANDOFF.md
                              section 8 requires.
Physics changed?:             No -- purely additive audit output; default
                              path unaffected (confirmed no audit file is
                              written when the flag is absent).
Regression result:            40 passed, 1 skipped (full focused suite +
                              new real-engine tests).
Real-run result:               4/4 new real-engine tests passed.
Next discriminating experiment: Once the frozen PF reference bundle is
                              available (see PF_REFERENCE_REGENERATION_CONTRACT.md),
                              repeat this exact real-engine test structure
                              but with the real PF row-176 prefracture
                              trajectory as the target and the production
                              front-engine/kernel configuration, to move
                              Gate 1 from "partial" to "pass" or to
                              surface the first real discrepancy.
```

### 2026-08-04 (earlier this session — infrastructure/blocker discovery)

```text
Highest completed gate:      None (Gate 1 is partial: controller numerics
                              validated synthetically; not yet validated
                              against the real FEM engine).
Current gate:                Gate 1 (controlled prefracture mechanics).
PF reference:                steps_1000K.csv statistics captured this
                              session (176 prefracture rows, first
                              passage at row 176/step 177, KJ=53.9019
                              MPa√m, J≈6531 J/m²) before the live copy
                              became unavailable.
FEM/CZM result:               Fixed-ramp baseline (pre-existing defect):
                              KJ=11.37 MPa√m at step 9400, extrapolated
                              ~145,000 steps to first passage. Controller-
                              driven result: none yet (blocked).
Agreement:                    Not assessable yet — no controller-driven
                              real run has completed.
Largest discrepancy:          The loading/controller gap itself (the
                              original defect this whole effort targets):
                              fixed dU/dt produces a far slower dKJ/dt
                              than the PF reference because FEM and PF
                              geometries have different compliance.
Likely subsystem:             loading/controller.
Evidence:                     Phase 0 archived + freshly reproduced
                              slow-ramp smoke (bit-identical at steps
                              200/400/600/800); PF row-176 first-passage
                              values located exactly from the (now
                              externally unavailable) reference CSV.
Change made:                  Wired a safeguarded predictor+secant
                              KJ-target controller into sharp_front.run_2d
                              behind a default-off --pf-kj-target-csv
                              flag (commit ad2c446). Composes with the
                              existing hazard-clock accept/reject check
                              without modifying it.
Physics changed?:             No — default path proven bit-identical;
                              the new path has not yet been exercised
                              against real FEM state.
Regression result:            36 passed, 1 skipped (skip = real-reference
                              integration test, correctly skipping since
                              the file is externally unavailable).
Real-run result:               None yet — blocked by missing/mutated PF
                              reference artifacts at their documented path
                              (see PF_REFERENCE_REGENERATION_CONTRACT.md).
Next discriminating experiment: Recover or regenerate the frozen PF
                              reference bundle (checksummed against
                              666d839a.../a85b57ad9e...), then run the
                              flag-ON controller against it for a short
                              prefracture interval and inspect whether
                              achieved KJ(t) tracks the PF target within
                              tolerance with a practical iteration count
                              (Gate 1), before attempting any Gate 2
                              internal-state comparison.
```
