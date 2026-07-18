# MPZ v9.10.4.7 analytical DBTT down-selection

This campaign restores a fidelity hierarchy before the expensive moving-interface calculation.

1. Sample the complete independent-shape parameter domain with a deterministic Sobol design.
2. Evaluate a zero-D first-event kinetic closure on 300--1100 K at 100 K increments.
3. Assign each viable candidate to its best adjacent 100 K DBTT bracket.
4. Retain the best candidates within every occupied bracket rather than allowing one transition temperature to dominate the shortlist.
5. Give each promoted candidate exactly four 1-D moving-interface temperatures spanning its own bracket: lower edge, one-third, two-thirds, and upper edge.
6. Rank the 1-D responses within each bracket and write the short-growth promotion manifest.

The analytical screen uses fixed zero direct cleavage temperature slopes by default. It retains finite source depletion, emission back stress, Peierls transport, Taylor retention/release, recovery, escape, blunting, cleavage shielding, and the cleavage hazard clock, but stops at the first event and does not translate the process zone or refresh sources.

Primary outputs:

- `analytical_all_candidates.csv`
- `analytical_promotion_manifest.csv`
- `analytical_bracket_summary.csv`
- `moving_interface_queue.csv`
- `dynamic_1d_all_candidates.csv`
- `dynamic_1d_temperature_detail.csv`
- `short_growth_promotion_manifest.csv`

Both stages checkpoint incrementally. The analytical stage checkpoints deterministic candidate batches. The 1-D stage checkpoints each completed candidate.
