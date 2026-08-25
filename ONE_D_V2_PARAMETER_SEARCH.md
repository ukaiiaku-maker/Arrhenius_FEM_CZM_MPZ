# One-dimensional V2 parameter search

## Result

The provider-robust DBTT search screened 187 shared rows (20 archived candidates, 39 controlled parameter morphs, and 128 Sobol-local perturbations). Sixteen survivors were rerun at six temperatures with three hazard seeds. **Zero rows passed the full shared DBTT contract.**

The non-dominated conflict is structural. Rows closest to the target avalanche counts exceed the qualified radius map and/or lose state realism; complete bounded rows sacrifice the DBTT transition or remain 2–3 avalanches away from the paired PF/FEM targets. Therefore no candidate is promoted by a hidden weighted compromise.

| candidate_id | status_objective | class_topology_objective | onset_envelope_objective | provider_robustness_objective | state_realism_objective | pareto_nondominated_validation | selection_status |
|---|---|---|---|---|---|---|---|
| oneD_v2_dbtt_morph_014 | 0 | 2 | 2.38335 | 0.00499482 | 0 | True | REJECTED_FULL_DBTT_CONTRACT |
| oneD_v2_dbtt_morph_013 | 0 | 2.33333 | 2.35572 | 0.00928686 | 0 | True | REJECTED_FULL_DBTT_CONTRACT |
| oneD_v2_dbtt_local_a1_033 | 0 | 2.33333 | 2.46887 | 0.115317 | 0 | False | REJECTED_FULL_DBTT_CONTRACT |
| oneD_v2_dbtt_local_a1_040 | 0 | 2.33333 | 20.9596 | 1.50879 | 0.0354732 | False | REJECTED_FULL_DBTT_CONTRACT |
| oneD_v2_dbtt_local_a1_005 | 0 | 2.33333 | 2.73075 | 0.0493093 | 0 | False | REJECTED_FULL_DBTT_CONTRACT |
| oneD_v2_dbtt_morph_010 | 0 | 2.66667 | 2.37023 | 0.0100019 | 0 | False | REJECTED_FULL_DBTT_CONTRACT |
| oneD_v2_dbtt_morph_012 | 0 | 2.66667 | 2.38629 | 0.00989021 | 0 | False | REJECTED_FULL_DBTT_CONTRACT |
| oneD_v2_dbtt_morph_011 | 0 | 2.66667 | 2.39136 | 0.00951012 | 0 | False | REJECTED_FULL_DBTT_CONTRACT |

The machine population is `analysis_outputs/oneD_v2_predictive_model/oneD_v2_search_population.parquet`; independent objectives and active parameters are in `oneD_v2_pareto_candidates.csv`.
