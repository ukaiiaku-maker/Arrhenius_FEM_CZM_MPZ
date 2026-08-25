# One-dimensional V2 parameter-sensitivity validation

The deterministic matrix contains 972 cases (Peak/DBTT, 600/1000/1200 K, both providers, ±10% and ±25%). Every case reaches 100 µm. Shared-material first-onset sign agreement averages 84.8% across groups; all dominant coordinates agree in sign.

| parameter_owner | sensitivity_group | sign_agreement | onset_span | envelope_span |
| --- | --- | --- | --- | --- |
| SHARED_MATERIAL_ROW | emission_stress_free_barrier | 1 | 0.331 | 0.333 |
| SHARED_MATERIAL_ROW | cleavage_stress_free_barrier | 1 | 0.288 | 0.291 |
| CONTROLLED_COMMON_PHYSICS | backstress_coefficient | 1 | 0.287 | 0.351 |
| CONTROLLED_COMMON_PHYSICS | process_zone_length | 1 | 0.239 | 0.28 |
| SHARED_MATERIAL_ROW | cleavage_activation_stress_shape | 1 | 0.136 | 0.143 |
| SHARED_MATERIAL_ROW | blunting_coefficient | 1 | 0.134 | 0.15 |
| BACKEND_REDUCTION | lifecycle_hazard_progress | 1 | 0.0514 | 0.703 |
| SHARED_MATERIAL_ROW | emission_activation_stress_shape | 0.667 | 0.0166 | 0.0518 |
| SHARED_MATERIAL_ROW | initial_source_density | 0 | 0.0113 | 0.0173 |
| CONTROLLED_COMMON_PHYSICS | source_zone_length | 1 | 0.000196 | 0.113 |
| SHARED_MATERIAL_ROW | peierls_shape | 0.833 | 3.26e-17 | 3.26e-17 |
| SHARED_MATERIAL_ROW | peierls_barrier_entropy | 0.833 | 3.25e-17 | 3.25e-17 |

Dominant onset ranking is emission barrier, cleavage barrier, backstress, process-zone length, cleavage stress/shape, and blunting. Initial source density is the only clear low-magnitude cross-provider sign disagreement. Taylor/Peierls coordinates are nearly inactive for first onset in this local window but can affect later state.

PF loading-rate direction is checked against existing direct PF rate-extreme runs at 1000 K: Peak decreases from slow to fast in both models; DBTT increases in both. Maximum absolute reduced/direct onset error across 0.01×, 1×, and 100× is 22.1%. FEM/CZM is qualified only through 1×; 100× reaches the opening bound and is excluded.
