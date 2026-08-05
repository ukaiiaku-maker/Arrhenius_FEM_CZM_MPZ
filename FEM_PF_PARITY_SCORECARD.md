# FEM/CZM ↔ PF Physical-Correspondence Scorecard

This file, not the current coding task, organizes the work. Every session
should update it before ending. See CLAUDE.md and FEM_CZM_HANDOFF.md for
the governing contract, and CLAUDE_PROGRESS.md for narrative session-level
detail (commits, exact commands, architecture notes). This file is the
scientific status ledger.

**Naming note (2026-08-04, corrected):** earlier entries in this file used
the word "parity" in the strict sense of near-exact numerical agreement.
That standard has been corrected — see "Acceptance philosophy" immediately
below. The FEM/CZM model is not required to numerically match the PF model
step by step or reproduce identical stochastic histories; it must exhibit
**physical correspondence** / **cross-model consistency** with PF. Read
older log entries below with that correction in mind rather than at face
value — they are kept as an honest historical record, not rewritten.

Reference case for all rows below unless stated otherwise: Peak
parameterization (`v913_paper_peak01_0242980_persistent_sites`), 1000 K,
theta=0, hazard seed 8666.

**Principal PF reference bank (2026-08-04, superseding the row above as the
default source of PF data):** the PF final v10.2.30 four-class campaign,
`PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1`
(3 rates x 4 classes x 12 temperatures = 144 cases, all verified
`complete_target_extension`). See `PF_FINAL_CAMPAIGN_REFERENCE_SUMMARY.md`
for the full characterization and `PF_FINAL_CAMPAIGN_REFERENCE_MANIFEST.json`
for per-case detail. The older `v10_4_1_..._selective_reuse_base3621_v1`
case referenced in earlier log entries below remains a valid provenance
record for what it audited at the time, but is no longer the primary
source — its own reference artifacts are unavailable at their documented
path (see CLAUDE_PROGRESS.md "EXTERNAL BLOCKER").

## Acceptance philosophy (governs everything below)

The PF and FEM/CZM formulations have different discretizations, stress
fields, crack representations, process-zone mappings, cohesive mechanics,
and numerical regularizations. Differences of order 10-20%, and
potentially somewhat larger for individual quantities, may be
scientifically acceptable when both models exhibit the same physical
mechanisms and broad response trends. **A result should not fail merely
because FEM differs from PF by 20%.**

The primary questions, in priority order:

1. Do both models predict the same fracture regime?
2. Do both exhibit or fail to exhibit R-curve toughening consistently?
3. Are first-passage or effective K_IC values comparable in scale?
4. Are the temperature trends and DBTT-like behavior consistent?
5. Is the ordering of Peak, DBTT, weak-T, and ceramic parameterizations
   preserved?
6. Do both models show comparable competition among plasticity, emission,
   shielding, and cleavage?
7. Are crack-extension events and developed-growth statistics physically
   comparable, even when individual events differ?
8. Can discrepancies be explained by the known differences between PF and
   FEM/CZM representations?

### Quantitative comparison bands (working diagnostic bands, not universal
### fixed tolerances — apply tighter or looser per quantity)

```text
Excellent correspondence:  <10% difference
Good correspondence:       10-20%
Acceptable model-scale:    20-30%, with consistent mechanism/trend
Significant discrepancy:   >30%, requiring explanation
Qualitative failure:       wrong trend, wrong regime, reversed ordering,
                           absent/present R-curve incorrectly, or a
                           physically inconsistent mechanism
```

Individual events are stochastic, so agreement may need to be assessed
over a regime of moderate crack extension rather than point-by-point.

Guidance by quantity:
- K_IC, first-passage K_J, and developed R-curve level: ~20% initial
  target (Good-to-Acceptable band).
- R-curve shape, temperature dependence, material-class ordering: judged
  primarily by trend and mechanism, not a percentage.
- Individual stochastic event times/lengths: not required to match.
- Ensemble means, distributions, correlations: compare these when
  stochastic variability matters, not single realizations.
- Internal quantities (B, N_em, MPZ populations, plastic work): diagnostic
  indicators. They need not agree numerically if the different
  representations nevertheless produce consistent macroscopic behavior.

### Role of controlled-K_J comparisons

The controlled-K_J(t) run (the transactional PF-reference loading
controller wired this session) is a **diagnostic experiment, not the
final definition of success**. Matching a prescribed PF K_J(t) does not
itself prove physical correspondence — it removes loading-rate/boundary-
condition differences so that comparing internal responses under the same
driving-force history becomes a fair, isolating test. It should determine
whether major differences arise from:

- external loading and boundary conditions;
- J-integral extraction;
- bulk plasticity;
- MPZ representation;
- emission kinetics;
- cleavage hazard;
- cohesive dissipation;
- crack-geometry evolution.

If FEM and PF internal trajectories differ under the same K_J(t),
determine whether the difference:

1. represents a genuine implementation error;
2. follows naturally from the different model formulations;
3. changes the macroscopic fracture response materially;
4. can be resolved only through ensemble/statistical comparison.

**Do not modify the FEM/CZM model merely to force internal variables to
overlay the PF curves**, and do not retune parameters merely to reduce a
numerical percentage difference when the mechanisms and macroscopic
conclusions already agree.

## Four categories of comparison

```text
Numerical correctness:
  Does each implementation satisfy its own equations and transactional
  contract? (Tight tolerances are appropriate here — this is not a
  cross-model comparison at all.)

Controlled-driving-force comparison:
  How do the internal responses differ under the same J/K_J history?
  (Diagnostic; see "Role of controlled-K_J comparisons" above.)

Physical correspondence:
  Do the models predict comparable toughness, R-curve behavior,
  temperature trends, material-class ordering, and mechanisms?
  (Judged against the comparison bands above.)

Statistical correspondence:
  Are event distributions and developed-growth statistics comparable
  within stochastic variability? (Ensemble/distributional, not
  event-by-event.)
```

For every comparison recorded in this file, report:
- absolute and relative differences;
- qualitative agreement (trend/regime/ordering);
- whether the same mechanism is active in both models;
- whether the discrepancy affects the macroscopic conclusion;
- whether ensemble sampling is needed to judge it fairly;
- whether the difference is acceptable, unexplained, or clearly erroneous.

## Status summary (2026-08-04, updated — PF campaign pivot + first real controlled-K_J diagnostic)

**Controller infrastructure complete and validated end-to-end against the
real FEM engine, first with a synthetic target (Numerical correctness
category), and now once against a REAL PF K_J(t) trajectory from the
principal 144-case reference bank (still Numerical correctness only — see
below for exactly why this is not yet Physical-correspondence evidence).
No Controlled-driving-force, Physical-correspondence, or Statistical-
correspondence comparison exists yet.**

- The PF driving-force trajectory reader and the safeguarded J/K_J-target
  controller exist, are unit-tested in isolation, wired into
  `sharp_front.run_2d` behind a default-off flag, and validated against
  the real FEM engine (real mesh/assembly/solve/plasticity, a fail-closed
  stochastic-identity check across RNG/threshold/event-length/B/N_em/MPZ
  state, and a direct-step mechanical-state-equality proof). The default
  (fixed-ramp) path is proven bit-identical to its pre-existing behavior.
- **New this entry — first real controlled-K_J run against actual PF
  data**: `sharp_front.main()` invoked directly (bypassing the
  kernel-gated production wrapper — no kernel, no persistent-site MPZ
  engine, no full-field bulk Peierls-Taylor closure; `front_state_model
  legacy_scalar` only) with the production-scale mesh (`nx=36, ny=72,
  tip_h_fine=1e-6, tip_ratio=1.20, da_phys=5e-6` → 1124 nodes,
  `hbar_tip=1.3846e-6 m`, matching the PF case's own
  `stage3_case_status.json` almost exactly) and the correct Peak-class
  cleavage/emission barrier row installed via `sharp_front.py`'s own
  `--cleave-*`/`--emit-*` CLI flags, targeting the frozen
  `reference_inputs/pf_final_v10_2_30/rate1x/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/steps_1000K.csv`
  trajectory. Result: 40/40 controlled steps converged (3-9 controller
  iterations, mean 5), achieved vs. target K_J agreed to ~7e-6 relative at
  the last step (13.587140 vs. 13.587047 MPa√m), `B` stayed exactly `0.0`
  throughout, `N_em` grew steadily (0.62 → 13.53), no premature failure or
  ligament severing (unlike two earlier attempts — a tiny 6x10 toy mesh,
  and the production mesh with generic/default barrier parameters — both
  of which severed the ligament within 6 steps at K_J~1-2 MPa√m; see
  `CLAUDE_PROGRESS.md` for the full failure record). Raw evidence:
  `runs/anchor_diagnostics/peak_1000K_rate1x_barriers/` (gitignored). Now
  reproducible via `scripts/run_pf_anchor_controlled_kj_diagnostic_v10051840.py`
  (see "Reproducible diagnostic script" below).
- **Why this is still Gate 1 / Numerical correctness, not Physical
  correspondence**: this run proves the controller can track a
  *prescribed* PF K_J(t) accurately with the real solver (a numerical
  root-finding/tracking question). It does **not** yet test whether
  FEM/CZM independently *predicts* a comparable K_J(t), first-passage
  toughness, or R-curve — that requires either (a) carrying this same
  configuration to first passage and treating the result as a reduced
  (kernel-free, MPZ-free) diagnostic, or (b) a native-loading run with no
  PF target at all (see "Recommended next experiment" below). Do not cite
  the ~7e-6 tracking agreement as FEM/PF physical correspondence — it is a
  statement about the controller's own root-finder, not about the two
  models' physics agreeing.
- **No Controlled-driving-force, Physical-correspondence, or Statistical-
  correspondence evidence exists yet for any anchor case.** See "Kernel
  provenance limitation" below for what specifically remains blocked
  (production/kernel-enabled physical correspondence) versus what does not
  (kernel-free controlled diagnostics, documented-kernel native-loading
  experiments, campaign-level macroscopic PF comparisons).
- A confirmed-and-fixed AST-patcher anchor regression (see
  CLAUDE_PROGRESS.md) is Numerical-correctness/engineering-practice
  evidence, not a physics finding — recorded in the session log below for
  completeness but not part of the physical-correspondence assessment.

## PF campaign reference bank (144-case characterization)

Full detail lives in `PF_FINAL_CAMPAIGN_REFERENCE_SUMMARY.md` and
`PF_FINAL_CAMPAIGN_REFERENCE_MANIFEST.json` — not reproduced in full here.
Key findings this scorecard uses as the target physical signature to
reproduce:

- **Peak**: broad, gradual transition ~800-950 K (23 → 57 MPa√m plateau).
- **DBTT**: sharp, narrow transition — flat ~21-23 MPa√m through 950 K,
  jumping to ~50-56 MPa√m exactly between 950 K and 1000 K.
- **weak-T**: essentially flat 19-21 MPa√m across the entire 300-1300 K
  range — no transition.
- **ceramic**: monotonically decreasing 13.6 → 9.3 MPa√m, no transition,
  lowest of the four classes throughout.
- **Ordering** `ceramic < weak-T < {Peak, DBTT}` at nearly every
  (rate, temperature).
- **Rate sensitivity**: weak at low T (300-600 K, ~1.1x spread,
  cleavage-dominated) and strong above ~800 K (1.7-2.35x spread,
  thermally-activated plasticity/emission regime) — slower loading gives
  higher first-passage toughness once that regime is active.

This is the physical-correspondence target signature for Gates 3, 7, and 8
below (ordering, transition shape/location, rate trend) — not a set of
numbers FEM must overlay point-by-point.

## Kernel provenance limitation (read before any kernel-enabled run)

Every PF case in the final campaign (and the older v10.4.1 case) refers to
the same live kernel-cache path
(`.../runs/v10_2_28_kernel_cache/1447653d.../family.json`), which was
observed to change content during concurrent, unrelated work in that
repository. No PF case's own audit JSON records an independent byte-level
kernel hash from generation time — only shape/policy metadata — so the
current live `family.json` (frozen into each `reference_inputs/.../family.json`
at time of use, with its own hash recorded) cannot be proven to be the
exact historical kernel bytes that generated that case's `steps_*.csv`.
See `CLAUDE_PROGRESS.md`'s "EXTERNAL BLOCKER" section and
`PF_REFERENCE_REGENERATION_CONTRACT.md` for the full incident record. **Do
not weaken the kernel hash guard** (`pf_theta0_frozen_reference_v10051840.py`)
and **do not describe the live `family.json` as verified-identical to the
original campaign kernel.**

What this blocks: describing any kernel-enabled FEM run (the full
production entry point, which hard-checks `REFERENCE_KERNEL_SHA256`) as
exact production replication of a specific historical PF case.

What this does **not** block:

- kernel-free controlled-K_J diagnostics (the Gate 1 run above uses no
  kernel at all — `legacy_scalar` front engine, no signed shielding);
- native-loading FEM experiments using a clearly documented, independently
  frozen or currently-live kernel snapshot (with its hash recorded and the
  caveat stated, not claimed as historically verified);
- physical-correspondence comparisons at the level of macroscopic PF
  outputs (first-passage K_J, R-curve level, temperature trend, ordering),
  provided the provenance limitation is stated explicitly;
- campaign-level comparison against the 144-case manifest's own recorded
  outputs (those are frozen PF *results*, not something this workspace
  regenerates from the kernel on read).

## Gate ladder status

| Gate | Description | Category | Status | Notes |
|---|---|---|---|---|
| 1 | Controlled prefracture mechanics (controller numerics) | Numerical correctness | **partial** | Predictor/secant logic validated by 11 synthetic-callback unit tests, real-engine tests against a synthetic target, AND (new) a real 40-step run against the actual PF Peak/1000K/rate1x/seed8666 K_J(t) trajectory with production mesh + correct barrier row (achieved-vs-target agreement ~7e-6 relative, all steps converged) — but still without the signed kernel, persistent-site MPZ engine, or full-field bulk Peierls-Taylor closure, and only 40 of ~163+ steps toward first passage. Still "partial," not "pass," until it is carried to first passage (or a native-loading run substitutes) under the full production configuration. This gate legitimately keeps a tight tolerance (the controller's own `rel_tol`) since it is not a cross-model comparison. |
| 2 | Prefracture internal-state comparison (B, N_em, plastic work, MPZ) | Controlled-driving-force comparison | **blocked** | Requires Gate 1 against real PF data. Acceptance is mechanism/trend-based, not numerical equality — see "Revised Gate 2 interpretation" below. |
| 3 | First-passage comparison (time, J, K_J, mechanism) | Physical correspondence (+ statistical, since stochastic) | **blocked** | Judge K_J,FP scale against the ~20% band; do not require exact time/value match. |
| 4 | First-event FEM/CZM transaction audit | Numerical correctness (of the FEM implementation's own contract) | **blocked** | Not a PF comparison — does the FEM transaction (hazard crossing → event proposal → energy gate → cohesive commit → renewal) satisfy its own stated contract. Requires the Gate 3/4 event-audit schema (defined, unpopulated — see `pf_theta0_first_passage_event_audit_v10051840.py`) to be wired to a real event. |
| 5 | Short-growth statistics (20-50 µm) | Statistical correspondence | **not tested** | Event-by-event identity not required; compare distributions/statistics over the interval. |
| 6 | Intermediate/long growth (400 µm, 1000 µm) | Statistical correspondence | **not tested** | Same, over a longer interval — developed R-curve level, not per-event match. |
| 7 | Four immutable material classes | Physical correspondence | **not tested** | Ordering and characteristic shape, not per-value fit. |
| 8 | Temperature dependence | Physical correspondence | **not tested** | Transition location/shape and ordering, judged by trend. |

### Revised Gate 1 interpretation

Gate 1 still requires the controller itself to track the requested K_J(t)
accurately. This is a **numerical-controller** requirement (does the
controller solve its own well-posed root-finding problem correctly) and
legitimately retains a tight tolerance — it says nothing about FEM/PF
physical correspondence by itself.

### Revised Gate 2 interpretation

Gate 2 does **not** require B(t), N_em(t), plastic work, or MPZ state to
match PF exactly. Compare them to identify mechanistic similarities and
differences. A Gate 2 result may **pass** when:

- the same mechanisms activate in the same general loading regime;
- the evolution is qualitatively consistent;
- differences have a defensible physical or representational explanation;
- the resulting first-passage toughness and subsequent growth behavior
  remain comparable.

## Recommended Stage-A anchor cases (from the 144-case characterization)

Nine cases, `rate1x`, chosen to bracket each class's characteristic
behavior (see `PF_FINAL_CAMPAIGN_REFERENCE_SUMMARY.md` for the rationale
behind each). Only the first has any evidence yet (Gate 1 partial, per
above) — the other eight are not yet started.

| Class | T (K) | Rationale | Status |
|---|---|---|---|
| Peak | 1000 | Plasticity plateau; pre-existing reference case | Gate 1 partial (40/40 real controlled-K_J steps against actual PF data; not yet first passage; no kernel/MPZ/bulk-PT) |
| Peak | 300 | Below transition, cleavage-dominated (~24 MPa√m) | not started |
| DBTT | 900 | Below the sharp transition (~22 MPa√m) | not started |
| DBTT | 1000 | Just above the sharp transition (~56 MPa√m) | not started |
| DBTT | 1100 | Above transition, on the plateau (~56 MPa√m) | not started |
| weak-T | 300 | Confirms flatness at low T (~19-20 MPa√m) | not started |
| weak-T | 1300 | Confirms flatness at high T (~20 MPa√m) | not started |
| ceramic | 300 | Confirms monotonic-decrease trend at low T (~13.6 MPa√m) | not started |
| ceramic | 1300 | Confirms monotonic-decrease trend at high T (~9.3 MPa√m) | not started |

Per governing instructions: do not process all nine before learning from
the first. Resolve the earliest scientifically important discrepancy in
Peak/1000K before freezing the remaining eight anchors.

## Physical-correspondence quantity table

Columns: PF reference value/trajectory; FEM/CZM value/trajectory;
comparison metric; comparison band applied (from the bands above, not a
fixed tolerance); status (not tested / excellent / good / acceptable /
significant-discrepancy / qualitative-failure / blocked); whether the
same mechanism is active in both; dominant discrepancy if any; likely
subsystem; supporting run/commit; exact next experiment.

| Quantity | PF reference | FEM/CZM result | Metric | Comparison band | Status | Same mechanism? | Dominant discrepancy | Likely subsystem | Supporting run / commit | Next experiment |
|---|---|---|---|---|---|---|---|---|---|---|
| Driving force J(t), K_J(t) (controller tracking, Numerical correctness only — not a PF comparison) | n/a for this row — see "Driving force under real PF target" below for the actual cross-model row | Synthetic-target real-engine test: achieved K_J matches target within ~2e-3 relative at every step | Relative error of achieved vs. target K_J(t) | N/A (numerical-controller tolerance, not a comparison band) | excellent (controller numerics only) | n/a | Real PF target not yet run | loading/controller | `ad2c446` (wiring); `99f6e5c`, `12b02e6` (real-engine evidence) | Recover frozen reference, then populate the row below |
| Driving force under real PF target (actual Controlled-driving-force-comparison row, OLDER v10.4.1 case) | `steps_1000K.csv` (176 prefracture rows before first passage; captured earlier this session, live copy currently unavailable) | Fixed-ramp baseline (pre-existing defect): reaches only K_J=11.37 MPa√m by step 9400 (extrapolated ~145,000 steps to first passage) | Achieved vs. target K_J(t) once the controller runs against real PF data | ~20% initial target once populated | **blocked** (superseded by the v10.2.30 campaign row below as the active anchor) | n/a yet | n/a yet | loading/controller | Phase 0 archived baseline | Superseded — use the v10.2.30 campaign anchor row below |
| Controller tracking of real PF K_J(t) (v10.2.30 campaign anchor, Peak/1000K/rate1x/seed8666, first 40 steps) | `reference_inputs/pf_final_v10_2_30/rate1x/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/steps_1000K.csv`, first 40 rows, KJ target rising to 13.587047 MPa√m | `sharp_front.main()` direct invocation, production mesh (1124 nodes), real Peak barrier row, no kernel/MPZ/bulk-PT: KJ achieved 13.587140 MPa√m at step 40, `B`=0.0 throughout, `N_em` 0.62→13.53 | Relative error of achieved vs. target K_J(t) | N/A (numerical-controller tolerance, `rel_tol=1e-3`) — **not yet a Physical-correspondence row**; see Status summary | excellent (controller numerics); **not applicable** as a physical-correspondence measurement yet | n/a (this row is Gate 1, not a cross-model comparison) | none — controller tracked the prescribed target within tolerance every step | loading/controller (Gate 1 only) | `runs/anchor_diagnostics/peak_1000K_rate1x_barriers/` (gitignored); reproducible via `scripts/run_pf_anchor_controlled_kj_diagnostic_v10051840.py` | Extend to first passage (~step 164, KJ≈53.9) under this same kernel-free configuration (Path A), or run native-loading with no PF target (Path B/native) to get an actual FEM-predicted K_J,FP for comparison against this same PF value |
| Bulk plasticity Wp(t), ep, rho, active volume | PF `W_bulk_plastic_cumulative_J_per_m` etc. exist in steps table | Not tested | Trend/mechanism comparison at matched physical time under Controlled-driving-force run; NOT required to overlay numerically | Diagnostic only — no fixed band; judge by mechanism | not tested | n/a | — | bulk plasticity / Peierls-Taylor kinetics | none | After Gate 1: extract PF `W_bulk_plastic_cumulative_J_per_m` and FEM equivalent, compare trend not value |
| Moving-tip MPZ populations/translation | Not yet extracted | Not tested | Population trend comparison | Diagnostic only | not tested | n/a | — | moving process zone | none | Same |
| Emission N_em(t), hazard/action | PF `N_em` column exists | Not tested | Trend/regime comparison | Diagnostic only | not tested | n/a | — | emission hazard | none | Same |
| Cleavage B(t), threshold, accumulated action | PF `B` column exists | Not tested | Trend/regime comparison | Diagnostic only | not tested | n/a | — | cleavage hazard | none | Same |
| First-passage: time, J, K_J | J≈6531 J/m², K_J≈53.9019 MPa√m at PF row 176 (step 177, t=1486.8s) — captured this session | Not tested (fixed ramp never reaches it in a practical step count; controller not yet run against real data) | Relative difference in K_J,FP scale (not exact time/value match) | ~20% initial target (Good-to-Acceptable) | **blocked** | n/a | n/a | loading/controller, then cleavage hazard | Phase 0 (fixed-ramp extrapolation); Phase 1 (`7562911`, exact PF row located) | Gate 3 once Gate 1/2 pass |
| Event proposal: length factor, physical length, endpoint | PF event-length draw (threshold_scaled, factor range 0.5-4.0) | Not tested | Distributional comparison over moderate extension, not per-event identity | Statistical correspondence — no fixed band, assess distribution overlap | not tested | n/a | — | event-length law | none | Gate 5/6, ensemble-level |
| Energy gate: available vs. required | Not yet extracted | Not tested | Trend/mechanism (does the gate activate in the same regime) | Diagnostic only | not tested | n/a | — | energy gate | none | Gate 4 |
| Cohesive transaction: trial/commit/veto/rollback identity | n/a (PF uses sharp-wake, no cohesive geometry) | FEM atomic path-corridor CZM backend exercised by pre-existing focused tests | Numerical correctness of the FEM implementation's own contract — not a PF comparison | N/A (own-contract check) | not tested (for controller-driven steps specifically) | n/a | — | cohesive representation / geometry-remeshing | `test_atomic_path_corridor_v10051839.py` (pre-existing, passing) | Gate 4 audit under controller-driven loading |
| Geometry: tip position, path, mesh quality | n/a | Existing mesh-quality gates pass in isolation | Numerical correctness (own-contract) | N/A | not tested (under controller) | n/a | — | geometry/remeshing | pre-existing tests | Gate 4/5 |
| Renewal: new threshold, B reset, MPZ update | Not yet extracted | Not tested | Mechanism comparison | Diagnostic only | not tested | n/a | — | cleavage hazard / MPZ | none | Gate 4 |
| Post-event continuation | n/a | Not tested | Does the next accepted step succeed after a commit (own-contract, pass/fail) | N/A | not tested | n/a | — | state transfer / restart | none | Gate 4 |
| Short growth (20-50 µm) event history | PF event statistics not yet extracted | Not tested | Ensemble statistics: event count, length distribution, spacing, J/K_J at events | Statistical correspondence — distributional overlap | not tested | n/a | — | multiple | none | After Gate 4 passes |
| Intermediate growth (400 µm) | Not yet extracted | Not tested | Developed R-curve level and shape, not per-event identity | ~20-30% band on level; trend on shape | not tested | n/a | — | multiple | none | After Gate 5 |
| Long growth (1000 µm) | Not yet extracted | Not tested | Same | Same | not tested | n/a | — | multiple | none | After Gate 6a |
| Four-class response (Peak/DBTT/weakT/ceramic) | PF ordering/characteristic behavior not yet extracted beyond Peak | Not tested | Ordering + characteristic shape, NOT per-value fit | Qualitative (ordering preserved = pass; reversed ordering = qualitative failure) | not tested | n/a | — | material-row transfer / common numerics | none | After Gate 6 |
| Temperature response (transition shape/ordering) | PF temperature sweep not yet extracted | Not tested | Transition location/shape and ordering, judged by trend | Qualitative | not tested | n/a | — | multiple | none | After Gate 7 |

## Reproducible diagnostic script

The first real controlled-K_J run (above) was performed ad hoc via a
one-off Python invocation. It is now reproducible via
`scripts/run_pf_anchor_controlled_kj_diagnostic_v10051840.py`, which:

- resolves the frozen anchor reference through
  `pf_theta0_frozen_reference_v10051840.resolve_frozen_pf_reference`
  (fails closed if the frozen `steps_1000K.csv`/`family.json` are missing
  or their SHA-256 no longer matches the recorded contract — never falls
  back to the live, externally-mutable PF `runs/` tree);
- assembles the exact `sharp_front` CLI argv (mesh, Peak barrier row,
  controller tolerance/iteration cap, `--pf-kj-target-csv`) from a small
  anchor registry (currently one entry: `peak_1000K_rate1x_seed8666_v10230`);
  new anchors are added by extending the registry, not by hand-editing argv;
- writes a `provenance_manifest.json` into the output directory before
  running (git commit, package version, resolved frozen-reference hashes,
  full assembled argv, anchor id, timestamp) so every run's inputs are
  independently auditable after the fact;
- supports `--dry-run` to validate/print the assembled configuration and
  provenance without invoking the FEM solver at all (used by this
  script's own tests — no long calculation is launched by testing
  configuration generation).

## Recommended next experiment (Path A vs. Path B)

Per the governing instructions' two candidate paths:

**Path A — extend the current kernel-free controlled diagnostic toward
first passage** (more steps of the same Peak/1000K/rate1x configuration
above, same mesh/barrier row, no kernel/MPZ/bulk-PT). **Path B — resolve
kernel provenance** well enough to run the full production
kernel-enabled configuration.

**Recommendation: Path A first.** Reasoning:

- Path A requires no new infrastructure — the script above already runs
  it; extending it is changing `--steps` and rerunning, then inspecting
  `B(t)`, `N_em(t)`, achieved-vs-target K_J, and whether prefracture
  evolution stays stable up to (and past) the PF first-passage row
  (~step 164, K_J≈53.9 MPa√m).
- It gives a concrete, fast (40 steps completed in the first run; ~4x more
  to reach first passage) answer to a real scientific question the
  governing instructions prioritize: does this reduced (kernel-free,
  MPZ-free) configuration remain numerically stable and qualitatively
  sensible all the way to first passage, or does it break down (as the
  two earlier failed configurations did, just faster)? That is
  information Path B cannot produce regardless of how it resolves.
- Path B (kernel reconstruction) is valuable but is an *infrastructure/
  provenance* task, not a physics result — per the governing instructions'
  own emphasis ("prefer obtaining a controlled first-passage and
  native-loading short-growth result over extended infrastructure work"),
  and per this workspace's own repeated experience this session that
  provenance work here has historically taken multiple sessions
  (`PF_REFERENCE_REGENERATION_CONTRACT.md` documents an exhaustive,
  still-unresolved search). Attempting it now would delay the first
  concrete physical-correspondence datapoint.
- Path A's result also directly informs whether Path B is worth pursuing
  soon: if the kernel-free diagnostic already shows credible prefracture
  behavior and a first-passage K_J broadly comparable to PF's 53.9
  MPa√m, that is itself useful reduced-model evidence even before the
  kernel question is resolved. If it instead reveals a new instability,
  that is a physics/numerics bug worth fixing before spending effort on
  kernel provenance at all.

This recommendation does not preclude Path B — it should still be pursued
in parallel or immediately after, since the kernel-enabled production
configuration is required for a claim of full physical correspondence
(MPZ shielding and bulk Peierls-Taylor plasticity are both currently
absent from the Path A diagnostic and are named non-negotiable mechanisms
in `FEM_CZM_HANDOFF.md`).

## Task 4 preparation: native-loading comparison (not yet launched)

The final scientific test is a native-loading FEM run (no PF K_J(t)
target at all) for Peak/1000K, independently predicting initiation K_J,
initiation J, fracture regime, first crack-extension event, and an
initial 20-50 µm R-curve — compared against the PF v10.2.30 Peak/1000K
case using the physical-correspondence bands above, not exact equality.

**Not yet launched this session** — per the governing instructions, do
not launch a large native-loading run until the controlled-K_J/first-
passage workflow above is physically credible (Path A's outcome). The
preparation step completed this session is: the mesh/barrier-row
configuration validated in the Gate 1 diagnostic (production-scale mesh,
correct Peak row, `sharp_front.main()` direct entry) is the same
configuration a native-loading run would reuse, minus `--pf-kj-target-csv`
and with the pre-existing fixed-ramp or a to-be-designed rate schedule in
its place — no new CLI surface is required, only a decision on the
opening-rate schedule once Path A's first-passage step count is known
(so the native run's own rate can be chosen to reach first passage in a
comparable, practical number of steps rather than repeating the original
145,000-step slow-ramp problem this whole effort was created to fix).

## Revised long-term acceptance

The main FEM/CZM objective should be considered successful when it
independently produces:

- physically stable first-passage fracture;
- energy-consistent cohesive crack growth;
- moving-tip plasticity and shielding;
- R-curve development comparable to PF where PF predicts it;
- comparable K_IC or effective toughness scale, initially targeting
  roughly 20%;
- consistent temperature-dependent trends;
- consistent material-class ordering;
- stable 20-50 µm, 400 µm, and eventually 1000 µm growth;
- robust restart and transactional behavior.

**Exact event-by-event or seed-by-seed identity is neither expected nor
required.** Do not retune parameters merely to reduce a numerical
percentage difference when the mechanisms and macroscopic conclusions
already agree.

## Session reporting log

Newest entry first. Use the required form from the governing instructions.
Entries before 2026-08-04's correction use the older, stricter "parity"
framing in places — read them through the "Acceptance philosophy" above,
not at face value.

### 2026-08-04 (continued — scorecard update + reusable diagnostic script)

```text
Highest completed gate:      None (Gate 1 still partial -- unchanged this
                              entry; no new physics run, only
                              documentation and a reproducibility script
                              wrapping the already-reported Gate 1
                              result).
Current gate:                Gate 1 (Numerical correctness).
Real PF reference available?: Yes, for the Peak/1000K/rate1x/seed8666
                              anchor from the v10.2.30 final campaign
                              (frozen and hash-verified). Still no
                              signed-kernel/MPZ/bulk-PT production
                              configuration exercised for any anchor.
Production configuration exercised?: No.
Largest discrepancy:         Not a new physics finding this entry --
                              this entry converted the ad hoc Gate 1
                              diagnostic into a reproducible script and
                              updated the scorecard's documentation of
                              the campaign pivot, the Gate 1 result, the
                              kernel-provenance limitation, and the nine
                              Stage-A anchors. See prior entry below for
                              the actual numerical result being
                              documented (40/40 converged steps, ~7e-6
                              relative K_J tracking, B=0, N_em 0.62->13.53).
Subsystem implicated:        Documentation (this file) + a new
                              reproducibility/provenance script; no
                              change to sharp_front.py or any physics
                              module.
Evidence:                    scripts/run_pf_anchor_controlled_kj_diagnostic_v10051840.py
                              reproduces the exact CLI argv of the
                              original ad hoc run (verified against
                              runs/anchor_diagnostics/peak_1000K_rate1x_barriers/run_args.json)
                              and was confirmed to (a) resolve the real
                              frozen anchor and produce byte-identical
                              argv via --dry-run, and (b) have that argv
                              accepted by sharp_front's own
                              _build_parser() without error. 13 new focused
                              tests (tests/test_run_pf_anchor_controlled_kj_diagnostic_v10051840.py)
                              cover: anchor registry lookup and
                              unknown-anchor failure, argv assembly
                              content and defaults, argv acceptance by
                              the real sharp_front parser, fail-closed
                              resolution on hash mismatch/missing files,
                              successful resolution on matching hashes,
                              provenance manifest contents, dry-run
                              never invoking the solver, fail-closed
                              ordering (bad reference caught before the
                              output directory is even created), refusal
                              to overwrite a non-empty output directory,
                              and an end-to-end --dry-run CLI smoke
                              against the real committed frozen anchor.
Change made:                 New file scripts/run_pf_anchor_controlled_kj_diagnostic_v10051840.py
                              and tests/test_run_pf_anchor_controlled_kj_diagnostic_v10051840.py.
                              FEM_PF_PARITY_SCORECARD.md updated: PF
                              campaign reference-bank pivot recorded,
                              Gate 1 real-diagnostic result documented
                              with an explicit Numerical-correctness-vs-
                              Physical-correspondence distinction, kernel-
                              provenance limitation section added, nine
                              Stage-A anchors tabulated with status, a
                              Path A/Path B next-experiment recommendation
                              (Path A: extend the current kernel-free
                              diagnostic toward first passage) added, and
                              a Task 4 native-loading preparation note
                              added (not launched -- explicitly deferred
                              until Path A's outcome is known).
Physics changed?:            No.
Regression result:           Focused suite (all PF-controller-family
                              tests + the new file): 74 passed, 3
                              skipped. Full tests/ suite: 740 passed, 3
                              skipped, 8 failed -- the same 8 failures as
                              the previously-documented pre-existing
                              baseline (5 stale package-version-string
                              assertions, 1 unrelated fatigue-anchor
                              break, 2 more from the same baseline
                              check); no new failures.
Real-run result:              13/13 new tests pass; no new FEM solve was
                              launched this entry (only --dry-run
                              config/provenance exercised, per
                              instruction not to launch a long
                              calculation while testing script
                              construction).
Next discriminating experiment: Path A (recommended above) -- rerun
                              scripts/run_pf_anchor_controlled_kj_diagnostic_v10051840.py
                              for the peak_1000K_rate1x_seed8666_v10230
                              anchor with --steps increased toward this
                              case's real first-passage row (~step 164,
                              K_J≈53.9 MPa√m), and inspect whether B(t),
                              N_em(t), and controller convergence remain
                              stable all the way there under this
                              reduced (kernel-free, MPZ-free)
                              configuration.
```

### 2026-08-04 (continued — acceptance standard corrected to physical correspondence)

```text
Highest completed gate:      None (Gate 1 still partial; nothing above
                              Gate 1 has been attempted, so the stricter-
                              vs-looser standard hasn't yet been exercised
                              against real data).
Current gate:                Gate 1 (Numerical correctness).
Real PF reference available?: No -- still blocked.
Production configuration exercised?: No.
Largest discrepancy:         None assessed this entry -- this was a
                              documentation/framework correction, not a
                              new experiment. Corrected FEM_PF_PARITY_SCORECARD.md
                              to replace strict numerical-parity language
                              with a physical-correspondence standard
                              (comparison bands: <10% excellent, 10-20%
                              good, 20-30% acceptable with consistent
                              mechanism, >30% significant discrepancy
                              requiring explanation, qualitative failure
                              for wrong trend/regime/ordering). Added the
                              four-category taxonomy (Numerical
                              correctness / Controlled-driving-force
                              comparison / Physical correspondence /
                              Statistical correspondence) and revised
                              Gate 1/2 interpretation text.
Subsystem implicated:        Documentation/acceptance-criteria framework,
                              not a physics or code subsystem.
Evidence:                    N/A (no new experiment this entry).
Commit(s):                   (see CLAUDE_PROGRESS.md for the exact hash
                              of this documentation update)
Exact next scientific experiment: Unchanged from the prior entry --
                              recover/regenerate the frozen PF bundle,
                              then run Gate 1 for real against the
                              production front-engine configuration.
                              Once real PF data is available, apply THIS
                              corrected standard (bands + mechanism
                              trends, not exact equality) when judging
                              Gate 2 and above -- do not hold internal
                              state variables to numerical-equality
                              standards the governing instructions
                              explicitly reject.
```

### 2026-08-04 (continued — stochastic-identity fingerprint + anchor regression found and fixed)

```text
Highest completed gate:      None (Gate 1 still partial).
Current gate:                Gate 1 (controlled prefracture mechanics).
Real PF reference available?: No -- still blocked, see EXTERNAL BLOCKER
                              and PF_REFERENCE_REGENERATION_CONTRACT.md
                              (now a complete, executable request
                              identifying the best-candidate PF commit).
Production configuration exercised?: No -- only the minimal
                              legacy_scalar front engine, via a synthetic
                              target. A fail-closed frozen-reference
                              interface (arrhenius_fracture/
                              pf_theta0_frozen_reference_v10051840.py)
                              and launcher plumbing are now ready for the
                              moment real artifacts arrive.
Largest discrepancy:         Not a physics discrepancy this entry -- a
                              numerical/engineering regression: the
                              earlier default-off wiring (commit ad2c446)
                              silently broke 7 tests in unrelated modules
                              by re-indenting text two AST/text-patch
                              anchors depended on matching verbatim. Not
                              caught by the narrower focused-suite runs
                              used throughout the session; only surfaced
                              running the full `tests/` suite.
Likely subsystem:            loading/controller wiring mechanics (not
                              physics) -- and test-coverage practice
                              (relying on a narrow suite instead of the
                              full one after editing shared code).
Evidence:                    A from-scratch AST-based anchor audit
                              across every module that patches
                              sharp_front.run_2d, cross-checked against a
                              temporary git worktree at the pre-session
                              baseline commit to separate the 7 real
                              regressions from unrelated pre-existing
                              failures. Fixed by no longer wrapping the
                              original mechanics/KJ block in any new
                              conditional; the controller now determines
                              its target Uapp via a small, purely local
                              search before that block, which runs
                              completely unconditionally and unchanged.
Change made:                 c998041 (the anchor fix) plus, earlier in
                              this entry's work: 12b02e6 (stochastic-
                              identity fingerprint + fail-closed
                              transactionality check across RNG/
                              threshold/event-length/B/N_em/MPZ state,
                              not just mechanics), f6fad45 (fail-closed
                              frozen-reference interface + completed PF
                              regeneration request), 4bb2d09 (Gate 3/4
                              event-audit schema, unpopulated).
Physics changed?:            No.
Regression result:           Full suite: 725 passed, 3 skipped, 8
                              failed -- all 8 independently confirmed
                              pre-existing at the baseline commit, none
                              caused by this session.
Real-run result:              9/9 real-engine tests pass (audit output,
                              rejected-trial recording, determinism,
                              flag-off no-op, direct-step state equality,
                              stochastic-identity preservation across
                              rejected trials, fail-closed tamper
                              detection).
Next discriminating experiment: Recover/regenerate the frozen PF
                              reference bundle per the now-complete
                              regeneration request, freeze it into
                              reference_inputs/, then run Gate 1 for
                              real: point --pf-kj-target-csv and the
                              production launcher's PF_FROZEN_REFERENCE_DIR
                              at it and inspect the controller audit +
                              stochastic-identity output against the
                              production front-engine configuration.
```

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
