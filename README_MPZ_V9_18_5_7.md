# MPZ v9.18.5.7 — subsegment-aware quality certification

The v9.18.5.6 ceramic 700 K, 60 µm run completed 12 physical cleavage
events, but the post-run runner rejected it because the quality audit contained
37 accepted records rather than 12.

That count is correct. The edge-aware adaptive-CZM backend realizes each
physical 5 µm event through one or more committed geometric subsegments. Its
recursive transaction chain records one accepted production-quality check per
committed subsegment/recursive level. In the audited run:

- 12 physical events were committed;
- those events contained 37 committed geometric subsegments;
- the v9.18.5.6 quality wrapper recorded 37 accepted transactions;
- every physical event's subsegments summed to exactly 5 µm;
- each event had a final v9.18.5.6 quality marker;
- there were no production-quality vetoes;
- minimum accepted triangle quality was 0.258546 (> 0.035);
- minimum child-area ratio was 0.212963 (> 0.08).

v9.18.5.7 changes no solver or constitutive physics. It keeps the
v9.18.5.6 explicit quality wrapper and replaces only the incorrect post-run
count assertion with a two-level certification:

1. **Physical-event level**
   - expected event IDs are present;
   - every event has the declared subsegment count and contiguous indices;
   - all subsegments are committed;
   - subsegment lengths sum to the physical `da`;
   - the final subsegment carries the event-level quality marker;
   - total committed length equals the requested target.

2. **Quality-transaction level**
   - accepted quality records equal committed subsegments;
   - all records pass triangle-quality and child-area floors;
   - no quality veto or consecutive-veto abort exists.

The certification is written to
`subsegment_aware_quality_certification_v91857.json`.
