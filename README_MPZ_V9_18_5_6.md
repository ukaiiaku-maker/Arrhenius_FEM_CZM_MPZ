# MPZ v9.18.5.6 — explicit production-quality wrapper chain

## Why this revision exists

The v9.18.5.5 ceramic 700 K / 60 um run completed, but its audit files exposed
that the advertised quality wrapper did not execute:

- 12 physical topology events committed;
- `active_tip_resolution_veto_guard_v91854.json` contained zero accepted events;
- `resolution_audit_only_consecutive_veto_v91855.json` contained zero resolution
  warnings and no accepted-event rows.

The cause was nested monkeypatching.  v9.18.5 assigns an `_original` attribute to
whichever function occupies its quality-wrapper slot.  In v9.18.5.5 that
assignment rebound the outer wrapper directly to the raw adaptive-CZM backend,
bypassing the inherited strict-quality function.

## Correction

`mode_i_first_passage_v9_18_5_6.py` installs one self-contained wrapper directly
into the v9.18.5 slot and enters the runtime through the v9.18.5.3 corridor chain.
The wrapper directly performs:

- finite positive element-area validation;
- minimum affected-triangle quality;
- minimum child/parent area ratio;
- orphan-node detection;
- cohesive-endpoint bulk-support validation;
- consecutive geometry-veto fail-fast handling.

The active-tip `h_tip/da` ratio remains recorded as a resolution warning and is
not used as a topology-validity criterion.

The sweep runner refuses to report success unless the per-case audit contains
exactly `TARGET_EXT_UM / PHYSICAL_DA_UM` accepted quality-gate rows and zero
production-quality vetoes.

No barrier, hazard, cohesive-opening, MPZ, wake, shielding, loading, or material
law changes are included.

## Interpretation of the audited v9.18.5.5 ceramic run

Direct inspection of the archived result indicates that the completed geometry
was healthy despite the wrapper bypass:

- initial/final global minimum triangle quality: approximately 0.1595;
- orphan nodes: zero;
- minimum local update quality recorded by the CZM advance log: approximately
  0.288;
- minimum recorded child-area ratio: approximately 0.213.

The 50 um same-load sequence is a model-predicted unstable cascade rather than
duplicate event bookkeeping.  Ten separate 5 um events complete while the
cleavage rate approaches the correlation-window shelf `1/tau_c = 1e6 s^-1` and
the adaptive loading increment becomes negligible.  It should be interpreted as
unstable ceramic propagation, not as a rising material R-curve.
