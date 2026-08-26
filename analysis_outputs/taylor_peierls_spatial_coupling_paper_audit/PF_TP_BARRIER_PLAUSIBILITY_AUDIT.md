# PF Taylor/Peierls barrier plausibility audit

The audit evaluates `ExpFloorBarrier.values_eV/rate`, `TransportBarrier.as_surface`, the rho-dependent Taylor multiplicity, encounter law, fixed emission surface, fixed cleavage surface, and cleavage multihit transform. It does not classify raw H0 or activation entropy in isolation.

For every candidate and checkpoint/onset, the 1100-K table records source-effective stress, free-energy barrier, rate, barrier-floor flag, rate-saturation flag, and inactivity flag. The temperature artifact repeats characteristic observed stresses from 300 to 1200 K in 50-K increments. `pf_tp_source_relevant_barrier_rate_ratios.parquet` evaluates Peierls/emission and Taylor/emission barriers plus all requested rate ratios on the same candidate, temperature, state, system, and bin.

- Control: EFFECTIVE_MODEL_EXTREME; flags ["FLOOR_DOMINATED", "RATE_SATURATED"].
- Transport: EFFECTIVE_MODEL_EXTREME; flags ["FLOOR_DOMINATED", "RATE_SATURATED"].
- Retention: EFFECTIVE_MODEL_EXTREME; flags ["FLOOR_DOMINATED", "RATE_SATURATED"].
- Near-tip: EFFECTIVE_MODEL_EXTREME; flags ["FLOOR_DOMINATED", "RATE_SATURATED"].
- Wake: EFFECTIVE_MODEL_EXTREME; flags ["FLOOR_DOMINATED", "RATE_SATURATED"].
- Mobile: EFFECTIVE_MODEL_EXTREME; flags ["FLOOR_DOMINATED", "RATE_SATURATED"].
- Backstress: EFFECTIVE_MODEL_EXTREME; flags ["FLOOR_DOMINATED", "RATE_SATURATED"].
- Balanced: EFFECTIVE_MODEL_EXTREME; flags ["FLOOR_DOMINATED", "RATE_SATURATED"].

`OUTSIDE_OBSERVED_STRESS_DOMAIN` is not assigned because classification uses only saved source-effective stress states. Diagnostic landscape curves extend only from zero to each candidate's observed maximum and do not extrapolate the shielding atlas. Large or negative activation entropy alone was not used as a rejection rule.

All eight loaded trajectories spend most sampled states on the Peierls/Taylor floor or rate asymptote. The physically interpretable shortlist is therefore empty and no four-row minimal future set is certified. This is a fail-closed handoff decision, not evidence that the eight full-set rows are unusable as diagnostic mechanism extremes.
