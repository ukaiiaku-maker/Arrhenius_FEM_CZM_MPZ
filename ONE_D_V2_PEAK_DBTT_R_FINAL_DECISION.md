# Final Peak/DBTT R decision

**Peak: RETAIN_CONTROL; NO_CREDIBLE_R_ENRICHED_ROW.** The retained row is `v913_zeroD_sobol_0242980`.

**DBTT: RETAIN_CONTROL; NO_CREDIBLE_R_ENRICHED_ROW.** The retained row is `v913_zeroD_sobol_0202500`.

The reduced Peak signal is FEMCZM-only, while direct PF turns the DBTT finalist's apparent upper-temperature toughening into a negative reload-separated onset increment. Consequently neither finalist is a shared-provider replacement or qualified optional R-curve row. The tested variants remain fracture-side diagnostic options in the separately versioned bank:

| target response class | candidate id | parameter sha256 | direct PF validation status |
|---|---|---|---|
| DBTT | oneD_v2_dbtt_R_3e168d9381eafa7b | 3e168d9381eafa7b2c41022459483033f17ec54892536368fc6a8feb98f57e84 | NOT_RUN |
| DBTT | oneD_v2_dbtt_R_41a3f2096b3c4c54 | 41a3f2096b3c4c542d57f9e11a9e93839ec2b2fb01474fae92a5dc82c749846a | NOT_RUN |
| DBTT | oneD_v2_dbtt_R_7b3b5b45354e64ef | 7b3b5b45354e64efd68fd402352db131e1fda37744ea66f0eea86e0a251b4490 | NOT_RUN |
| DBTT | oneD_v2_dbtt_R_9f5160f509e713e2 | 9f5160f509e713e256ea747e2549d05063d6a269258810b4348dcb4bc7fa9d59 | FAILED_NO_POSITIVE_RELOAD_SEPARATED_RESISTANCE |
| DBTT | oneD_v2_dbtt_R_c2dabb42101cf812 | c2dabb42101cf8128031f0f082da14fb74123ba518e9982f23bb277d4b140a5d | NOT_RUN |
| DBTT | oneD_v2_dbtt_R_c4859a34963f15af | c4859a34963f15af8f0a768e9948a98385802b266845d1c9166f86669cb137c0 | NOT_RUN |
| Peak | oneD_v2_peak_R_41f8789bcbc1f097 | 41f8789bcbc1f097eae880694e11ac6711dfbb6884b86b2b78981bfdb22943fb | FAILED_NO_POSITIVE_RELOAD_SEPARATED_RESISTANCE |
| Peak | oneD_v2_peak_R_7834115ae79cd559 | 7834115ae79cd5599062495367866062658f85c02160a079d50b6ce39104d6f6 | NOT_RUN |
| Peak | oneD_v2_peak_R_a2e923a7587543db | a2e923a7587543dbcc24cd621145bbf3a6f2ba2c749a6bbd9c454558726210d2 | NOT_RUN |

This is a decisive screening outcome, not a failed search. The retained DBTT control already contains modest precursor topology—two PF physical avalanches and three FEM/CZM physical avalanches at 1000 K—whereas Peak has one avalanche in each provider. The result is specifically that no material-vector change produced a more positive, provider-robust, direct-PF-validated reinitiation envelope.

The 1-D model performed its intended role: it identified apparent R-propensity, separated reload-separated onset states from within-avalanche drive, rejected one-provider signals, and selected bounded direct-PF transfers. Direct PF then showed that the apparent positive reduced signal did not survive the sharp-wake geometry. Within the explored 29-coordinate space and current architecture, further local Peak/DBTT retuning is not warranted; a future rising-resistance study would require a new physically motivated persistent wake, path-memory, shielding, hardening, or renewal mechanism.

The focused R-curve registry contains the two retained controls plus the nine explicitly screened variants, labeled `DIAGNOSTIC_OPTION` and `NOT_PROMOTED_DIAGNOSTIC_ONLY`. None changes a production material row.

No FEM/CZM material row is changed; no new FEM/CZM simulation was run. Weak-T and ceramic-like selected rows remain unchanged. Fatigue was not evaluated.
