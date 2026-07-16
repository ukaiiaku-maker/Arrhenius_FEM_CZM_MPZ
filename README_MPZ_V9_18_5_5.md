# MPZ v9.18.5.5 — resolution audit-only and consecutive-veto guard

## Failure reproduced

The v9.18.5.4 700 K ceramic gate reached 15 um projected extension and then
rejected every subsequent exact-target event solely because
`active_tip_h_over_da = 1.311 > 0.75`.  The triangle-quality, child-area,
finite-area, orphan-node, and cohesive-endpoint checks did not report failures.
The run therefore restored the same pre-event state and retried the same
unconsumed renewal indefinitely.

The signature-based veto counter also failed because small runtime changes in
its signature prevented the count from reaching the configured limit.

## Correction

`h_tip / da` is retained as a reported resolution warning but is not used as a
topology-validity veto.  The adaptive CZM backend inserts the requested physical
endpoint exactly and the strict gate still enforces:

- finite positive triangle areas;
- minimum accepted triangle quality;
- minimum accepted child-area ratio;
- no orphan bulk nodes;
- supported cohesive endpoints.

These are the mechanical admissibility checks for the local topology
transaction.  The process-zone resolution remains separately reported.

Every consecutive rejected geometry transaction is now counted, independent of
small changes in p0, p1, or the formatted reason.  The counter resets only after
an accepted event and raises at `ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES`.

## Audit

Each case writes:

`resolution_audit_only_consecutive_veto_v91855.json`

The audit lists all accepted events whose measured `h_tip / da` exceeded the
requested warning threshold and records any consecutive-veto abort.

## Physics scope

No barrier, hazard, cohesive-opening, MPZ, wake, shielding, material, loading,
or crack-increment parameter is changed.
