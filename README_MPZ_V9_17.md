# MPZ v9.17: absolute hazard opening and opening-proportional source refresh

v9.17 addresses the decisive v9.16 failure: the trial interface opened on a
prescribed microsecond clock while the calibrated cleavage renewal time was of
order seconds, and fresh source sites were withheld until after the event.

## Constitutive protocol

For one nucleated physical crack quantum:

1. Enforce exactly one cleavage renewal (`--max-advances-per-step 1`).
2. Insert the initially intact adaptive-CZM trial segment.
3. Advance cohesive progress with the absolute calibrated cleavage hazard,

   `dq/dt = lambda_c(K, S, T)`,

   until the accumulated opening hazard reaches one.
4. Refresh the depleted source inventory incrementally with opened surface:

   `refresh_fraction = min(q * da / L_refresh, 1)`.

5. Evolve emission, transport, retention, recovery, shielding, FEM equilibrium,
   and cohesive opening over the same physical substeps.
6. Translate the moving MPZ frame and commit physical crack extension only after
   `q = 1`; the legacy end-of-advance source refresh is suppressed to avoid
   double refreshing.
7. If one target hazard increment requires longer than the external loading
   interval, resume external loading. The event restarts when the absolute
   hazard timescale becomes shorter than that interval. No normalized rate-ratio
   arrest threshold is used.

There is no `tau_event`, `lambda/lambda_nucleation` progress law, rate exponent,
or empirical arrest/resume ratio in the v9.17 protocol.

## Diagnostic gate

`audit_matched_stress_classes_v917.py` evaluates each selected manifest at its
own observed initiation stress before running FEM. It reports:

- absolute raw and multihit cleavage rates;
- one-renewal opening time;
- per-site and full-inventory emission times;
- opening/emission timescale ratio;
- source refresh fraction per 5 um advance;
- retained lifetime;
- emitted, mobile, retained, and shielding state after 100 and 4000 s.

This audit is intentionally allowed to expose a retention failure. v9.17 fixes
the event architecture; it does not silently refit the existing manifests.

## Initial gate

```bash
CONDA_ENV=arrhenius-fem-czm \
SEEDS="1" \
CLASSES="ceramic weakT DBTT" \
T_K=700 \
TARGET_EXT_UM=10 \
BULK_PLASTICITY_MODE=tip_only \
EVENT_TARGET_DQ=0.05 \
OUTROOT=runs/mpz_v9_17_hazard_clock_source_refresh_700K_10um_v1 \
AUDIT_OUT=runs/mpz_v9_17_matched_stress_audit_700K_v1 \
bash run_mpz_v9_17_hazard_clock_source_refresh_700K.sh
```

Inspect each case's:

- `absolute_hazard_event_relaxation_v917.json`
- `absolute_hazard_event_audit_v917.json`
- `v9_17_case_summary.csv/json`

The immediate mechanistic gate is that weakT and/or DBTT must show source refresh
and consequential event-window emission. A rising long-range R-curve is not
expected until the retention/recovery parameterization is separately corrected.

## Known remaining limitations

- The J/K driving field is still referenced by the inherited moving front rather
  than an explicit cohesive-damage front. Dual old-tip/cohesive-front integrals
  remain a required audit.
- Trial topology is inserted before physical commit. MPZ translation is deferred,
  but topology is not fully transactional.
- The inherited terminal topology guard may leave a zero-progress trial event
  beyond the requested committed analysis window; post-processing continues to
  treat committed extension as authoritative.
- Existing weakT/DBTT retention and recovery parameters may erase shielding over
  the inter-event loading time. This requires a retention-aware refit, not an
  event-controller patch.
