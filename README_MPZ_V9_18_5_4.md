# MPZ v9.18.5.4 — active-tip resolution and bounded vetoes

The v9.18.5.3 ceramic 700 K gate reached two accepted advances and then entered
a deterministic geometry deadlock.  The rejected third event repeatedly reported

```
v9185_quality_veto:tip_h_over_da=1.097256>0.75
```

The run log contained thousands of identical vetoes because the v9.18.5.3 entry
path bypassed the separate bounded wrapper.

## Root causes

1. The event gate used `new_mesh.hbar_tip`.  That value is computed from the
   nearest two percent of all bulk elements and can broaden after topology
   updates.  It is not a direct resolution measure at the active cohesive
   endpoint.
2. v9.18.5.3 routed through v9.18.5.2 and v9.18.5.1 directly to v9.18.5, so
   `mode_i_first_passage_v9_18_5_bounded` was not active.

## v9.18.5.4 behavior

- Computes the mean length of unique, nonzero bulk edges in the geometric
  one-ring of the active cohesive endpoint.
- Uses that local value for the `h_tip / da` production gate.
- Stores the accepted local value back into `mesh.hbar_tip` for subsequent
  geometry tolerances and diagnostics.
- Integrates the identical-veto counter directly into the active strict-quality
  wrapper.  The default limit is 12.
- Retains triangle-quality, child-area, orphan-node, and cohesive-endpoint gates.
- Does not change barriers, hazards, cohesive opening, MPZ, wake, shielding, or
  material parameters.

## New audit

Each case writes:

```
active_tip_resolution_veto_guard_v91854.json
```

It records accepted-event one-ring metrics, all v9.18.5.4 vetoes, the legacy
stored `hbar_tip`, and any identical-veto abort.

## Acceptance gate

For the 700 K ceramic 60 µm diagnostic:

- the third event must not be rejected solely because the broad legacy
  `hbar_tip / da` is 1.097;
- accepted events must have `active_tip_h_over_da <= 0.75`;
- no geometry reason may repeat more than 12 times;
- the run must either reach the committed 60 µm target or fail explicitly with
  an auditable geometry reason.
