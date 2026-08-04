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

## Status summary (2026-08-04)

**Controller infrastructure complete; physical parity untested.**

- The PF driving-force trajectory reader and the safeguarded J/KJ-target
  controller exist, are unit-tested in isolation (18 synthetic tests
  total across both modules), and are wired into `sharp_front.run_2d`
  behind a default-off flag. The default (fixed-ramp) path is proven
  bit-identical to its pre-existing behavior.
- **No real flag-ON run against FEM/CZM state has completed.** The PF
  reference artifacts (`steps_1000K.csv`, kernel `family.json`) that this
  requires are currently unavailable at their documented path — see
  CLAUDE_PROGRESS.md's "EXTERNAL BLOCKER" section and
  `PF_REFERENCE_REGENERATION_CONTRACT.md` for full detail and the recovery
  plan.
- **Therefore every row below is "not tested" or "blocked."** No claim of
  FEM/PF physical parity — controlled-driving-force or native-loading — has
  been demonstrated at any gate above Gate 0 (numerical controller
  correctness, tested only against synthetic callbacks, not the real FEM
  engine).
- Important reminder from the governing instructions: matching a
  prescribed PF `KJ(t)` does **not** by itself prove parity. It only
  removes loading-rate/boundary-condition differences so that Gate 2
  (internal-state comparison at matched driving force) becomes a fair
  test. Do not stop at "the controller runs."

## Gate ladder status

| Gate | Description | Status | Notes |
|---|---|---|---|
| 1 | Controlled prefracture mechanics (controller numerics) | **partial** | Predictor/secant logic validated by 11 synthetic-callback unit tests + code review of the `sharp_front.py` wiring (transactionality argued from the existing `u_saved/ep_saved/rho_saved` scaffolding, not yet independently measured against a real FEM engine run). No real-engine integration test exists yet. |
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
| Driving force J(t), KJ(t) | `steps_1000K.csv` (176 prefracture rows before first passage; captured this session, live copy currently unavailable — see blocker) | Fixed-ramp baseline reaches only KJ=11.37 MPa√m by step 9400 (extrapolated ~145,000 steps to first passage). Controller-driven result: none yet. | RMSE / relative error of achieved vs. target KJ(t) | not yet defined (suggest ≤ controller `rel_tol`=1e-3 by construction once running) | **blocked** | n/a — no controller run yet | loading/controller | Phase 0 archived baseline; `ad2c446` (default-off wiring) | Recover frozen reference (see regeneration contract), then run Gate 1 |
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

### 2026-08-04 (this session)

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
