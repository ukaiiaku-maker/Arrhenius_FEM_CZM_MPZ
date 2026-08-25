# oneD V2 to cyclic-fatigue integration handoff

## 1. Purpose

Carry the completed one-dimensional V2 screening result into the PF/sharp-front fatigue code **without replacing or weakening the v9.14 reversible-dislocation physics**.

This is not a direct copy of the monotonic reduced driver. The fatigue model has cyclic waveform state, signed transport, physical surface return, stochastic first-passage thresholds, cycle acceleration, post-first-passage energy gating, and atomic crack/MPZ transactions that must remain authoritative.

## 2. Source result to import

Local predictive-program provenance:

- branch: `codex/oneD-v2-terminal-predictive-program`
- reported HEAD: `5045fee`
- producer-code commit in the final manifest: `cddc51605aee93d8ecbaa0dda76c12085e08f9de`
- final result status: `COMPLETE`

Final shared **screening** rows:

- Peak: `v913_zeroD_sobol_0242980`
- DBTT: `v913_zeroD_sobol_0202500`
- weak-T: `oneD_v2_focused_weak_T_0016`
- ceramic-like: `oneD_v2_focused_ceramic_like_0018`

The exact parameter vectors are provided in `data/oneD_v2_terminal_predictive_program/oneD_v2_new_four_class_registry.csv` and in the fatigue handoff registry.

## 3. What may be transferred

Transfer only the shared material-physics coordinates and their exact full-precision values:

- cleavage barrier and temperature/stress-shape parameters;
- emission barrier and temperature/stress-shape parameters;
- Peierls and Taylor barrier/entropy/shape parameters;
- initial source density and Taylor correlation coordinates;
- blunting coefficient and source/process-zone geometric reference values where already part of the shared material row.

Keep the new registry versioned and separate from the historical canonical registry. Do not overwrite the v9.13 production registry until the fatigue integration and cyclic validation pass.

## 4. What must not be transferred as material parameters

The oneD V2 program contains backend-reduction coordinates used to make monotonic PF- and FEM/CZM-consistent reduced models practical. These are not material parameters and must not replace the fatigue code's native cyclic algorithms:

- reduced event-length or translation closures;
- reduced hazard-progress or reload thresholds;
- reduced renewal/avalanche grouping rules;
- straight-path mechanics-map assumptions;
- exact-oracle cache or surrogate settings;
- PF-versus-FEM/CZM backend-specific lifecycle corrections.

The fatigue implementation must retain its existing cyclic waveform, persistent-site MPZ, stochastic threshold stream, threshold-correlated event-length law, energy gate, sharp-wake geometry, and high-cycle acceleration architecture.

## 5. Non-negotiable v9.14 reversible-dislocation physics

Preserve the v9.14 minimal-reversible-fatigue and v4 physical-return semantics.

### 5.1 Separate opening and transport channels

- Cleavage remains opening-only and uses the existing parent law.
- New dislocation emission remains opening-only and uses the existing parent law.
- Already-mobile signed line content uses the **unclipped** transport drive `K_transport = K_applied - K_shield`.
- Signed GND/internal stress is added to the transport stress.
- Signed Peierls mobility determines transport direction and speed.
- Do not reintroduce `max(K_transport, 0)` into the mobile transport channel.
- Do not create a reverse-nucleation channel unless separately authorized.

### 5.2 Physical return rule

A left-boundary mobile outflow reverses the emission-linked blunting ledger only when all three conditions hold:

1. the returning population is the Burgers-sign population emitted by the tensile crack-tip source on that system;
2. the effective transport stress at the crack/free-surface boundary is truly reversed relative to the positive-tension direction; and
3. positive left-boundary outflow is present.

Thus:

`physical return = emitted-sign population AND true reverse drive at x=0 AND left-boundary outflow`.

Raw left-boundary outflow without true reverse drive remains diagnostic only.

### 5.3 Irreversible and reversible ledgers

- Returned mobile line content may annihilate at the crack/free surface.
- Far-field escape remains a separate boundary fate.
- Retained/tangled content is **not** removed by the surface-return rule.
- Physical return cancels source-linked blunting only up to available uncancelled source slip.
- Preserve the cumulative source-slip ledger.
- Preserve separate cumulative returned-slip and nonnegative net-slip/blunting ledgers.
- Never allow returned slip to exceed cumulative source slip in net effect.

### 5.4 Cyclic transactionality

Preserve:

- waveform phase and same-cycle continuation after crack events;
- stochastic cleavage threshold and continuing RNG stream;
- event-length sampling and energy-gate semantics;
- atomic geometry commit and MPZ translation using the same admitted event distance;
- checkpoint/restart of signed state, cumulative ledgers, thresholds, RNG state, cycle phase, geometry signature, and high-cycle cache;
- separation of reversible active state from monotone cumulative diagnostic ledgers.

## 6. Integration architecture

Use a layered integration rather than merging monotonic lifecycle code into the fatigue engine.

Recommended structure:

1. `oneD_v2_four_class_screening_registry.csv`
   - versioned material rows only.
2. `FatigueMaterialRowAdapter`
   - exact field mapping from the V2 registry to existing fatigue barrier and MPZ configuration objects.
3. Existing fatigue cyclic engine
   - remains owner of signed cyclic transport, threshold lifecycle, event length, energy gate, geometry transaction, and cycle acceleration.
4. Optional sensitivity metadata
   - may inform screening order, but must not change physics.

Do not import the monotonic V2 backend lifecycle surrogate into the cyclic engine.

## 7. Required mapping audit

For every registry field record:

- V2 field name and units;
- fatigue target field/function;
- exact conversion, if any;
- whether it is material physics, shared geometry, inactive, or unsupported;
- default used by historical fatigue rows;
- before/after full-precision value;
- source hash.

Fail closed on an unmapped active coordinate. Keep backend-reduction parameters out of the material mapping.

## 8. Required regression ladder

### Gate A: registry and monotonic identity

- Exact full-precision registry parse for all four rows.
- Peak and DBTT historical rows remain bit-identical.
- New weak-T and ceramic rows reproduce the bounded monotonic PF transfer anchors within the recorded errors.
- No canonical registry overwrite.

### Gate B: no-reverse parity

At `R=0.1`, for any case with zero true reverse-drive fraction:

- baseline and reversible models must match in tip radius;
- cumulative cleavage hazard must match;
- event history and crack extension must match;
- raw left-boundary outflow must not cancel blunting.

### Gate C: signed transport and physical return

Use negative-R or an internally reversed transport state to prove:

- negative `K_transport` is retained;
- emitted-sign mobile velocity reverses;
- physical return becomes positive only under the three-condition rule;
- retained content is unchanged by the surface-return operation;
- returned-slip cancellation is bounded by available cumulative source slip.

### Gate D: cycle acceleration

Before production DMD/Poincare acceleration:

- explicit-cycle and accelerated trajectories must agree for signed active state;
- returned/far-field ledgers must agree;
- waveform phase must agree;
- hazard action and threshold crossing must agree;
- the accelerator must not positivity-project a signed reversible state;
- cumulative nonnegative ledgers may retain positivity-preserving treatment separately.

### Gate E: fatigue event transaction

- first passage uses the existing stochastic threshold;
- event-length proposal uses the existing threshold-correlated law;
- energy gate truncates the proposal without changing first passage;
- geometry and MPZ translation commit atomically;
- cycle continuation after the event is exact;
- checkpoint/restart reproduces event history and all reversible ledgers.

## 9. Initial fatigue validation matrix

Do not begin with a broad search. Use the final four shared rows in a bounded matrix:

- temperature: 300 K;
- fixed local `DeltaK`;
- `R = 0.1` baseline plus at least one negative-R reversibility diagnostic;
- frequency: 1000 Hz initially;
- identical seed set across rows and load points;
- explicit physical cycles for the first reversibility qualification;
- then the qualified high-cycle accelerator.

Suggested first cases:

- Peak and DBTT at one active `DeltaK` and one near-threshold `DeltaK`;
- weak-T and ceramic-like at one discriminating active load;
- one negative-R case for each of Peak and weak-T to exercise return physics.

Use event-resolved `da/dN`, moving-window developed rates, avalanche grouping, return/escape fractions, and signed-state histories. Do not fit a Paris law into the engine.

## 10. Sensitivity information to carry forward

The V2 972-case study found 84.8% overall PF/FEMCZM onset-sign agreement. Dominant directions agreed for cleavage and emission barrier coordinates, cleavage stress/shape coordinates, backstress, process-zone length, and blunting. Initial source density was the principal low-magnitude disagreement.

Use this sensitivity information to prioritize fatigue screening, but do not freeze parameters or alter cyclic state equations from the sensitivity table alone.

## 11. Domain and interpretation

The reduced screening rows were established over:

- 300–1200 K;
- unconditional 0–100 µm screening;
- PF rate factor 0.01–100;
- FEM/CZM rate factor 0.01–1.

The fatigue qualification has its own cyclic domain. Do not claim the monotonic screening domain automatically validates cyclic reversibility or VHCF acceleration.

These rows are response-class screening parameterizations, not uniquely identified material constants.

## 12. Branch and merge policy

Start from a branch that already contains the v9.14 signed reversible-transport and v4 physical-return corrections. Do not base the integration on an older fatigue branch and then reapply reversibility manually.

Recommended integration branch:

`codex/v9.14-oneD-v2-fatigue-integration`

Base:

`codex/v9.14-minimal-reversible-fatigue` or a verified descendant containing the same reversible commits and tests.

Keep the integration as a separate branch/PR until all gates above pass. Do not overwrite current production parameter options.

## 13. Final acceptance

The handoff is complete only when:

1. all four V2 material rows load at full precision;
2. no fatigue lifecycle or reversibility law has been replaced by monotonic V2 reduction logic;
3. no-reverse parity passes;
4. physical return passes under true reverse drive;
5. explicit-cycle and accelerated reversible-state parity passes;
6. event/energy/geometry transactionality passes;
7. the new rows are available through a versioned fatigue registry;
8. results and provenance are recorded without modifying the historical registry.
