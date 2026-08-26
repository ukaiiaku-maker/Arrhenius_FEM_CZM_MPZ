# Fracture-to-future-cyclic Taylor/Peierls handoff

`Taylor_Peierls_Microstructure_Option_Bank_v1` contains eight fracture-derived material options. Their fatigue response is unknown: every row has `fatigue_evaluated=false` and `fatigue_validation_status=NOT_EVALUATED`.

The rows were retained because Taylor/Peierls kinetics produce strongly different 2-D mobile/retained and wake states while giving nearly identical monotonic fracture topology and reload-separated softening. This makes them controlled options for a future test of whether cyclic loading is more sensitive to transport/retention partitioning than monotonic fracture.

- Complete mechanism set: all eight candidates.
- Physically interpretable shortlist (0): none passed the source-derived gate.
- Minimal future test set (0): not formed because fewer than four rows passed.

Future hypothesis — `HYPOTHESIS_NOT_EVALUATED`: under monotonic growth, rapid advance and current future-tip coupling suppress sensitivity to mobile/retained partitioning. Under repeated subcritical loading at a slowly advancing or stationary tip, the same transport/retention differences may accumulate into different blunting, shielding, return, or hazard histories. The future result may show strong, weak, or no divergence.

No PF lifecycle setting, mesh/wake control, reduced event rule, observer option, numerical tolerance, cache key, or avalanche-grouping rule is stored as a material coordinate. No fatigue repository or run was touched.
