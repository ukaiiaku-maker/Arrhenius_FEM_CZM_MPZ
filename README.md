# Arrhenius FEM/CZM MPZ v9.11.1 response-mechanism options

The reduced spatial calibration and 2-D FEM/CZM validation now support named weak-temperature, primary DBTT, peak-type, broad-shielding, intrinsic-control, and moderate-shielding response options. The repository preserves these alternatives because a DBTT-like macroscopic trend does not identify one unique microscopic mechanism.

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.6 uncapped, detailed-balance Peierls–Taylor constitutive audit
3. v9.7 independent-entropy calibration audit
4. v9.8/v9.8.1 joint response optimization audit
5. v9.9 barrier-scaling and spatial-mismatch audit
6. [v9.10 unified Peierls-transport/Taylor-retention search](README_MPZ_V9_10_UNIFIED_SEARCH.md)
7. v9.10.1 one-shape-for-all-four audit
8. [v9.10.2 independent four-barrier EXP-floor search](README_MPZ_V9_10_2_INDEPENDENT_SHAPES.md)
9. v9.10.2 independent-shape spatial MPZ promotion
10. [v9.10.3 target-aware DBTT search](README_MPZ_V9_10_3_DBTT_TARGETED.md)
11. 2-D FEM/CZM validation of candidates whose response persists during crack growth
12. [v9.11.1 DBTT and peak response-mechanism options](README_MPZ_V9_11_1_RESPONSE_OPTIONS.md)

Version 9.10 removed the independent trap barrier: Peierls kinetics set mobile transport, forest encounters convert mobile lines to retained lines, and Taylor completion releases retained lines. Therefore `H_P < H_T` retains the intended physical interpretation.

Version 9.10.2 released the complete independent shape space:

```text
(alpha_c, n_c)
(alpha_e, n_e)
(alpha_P, n_P)
(alpha_T, n_T)
```

Every restart begins from a fresh full 29-dimensional Sobol population over common broad bounds. No previous shortlist is used as an initial population. The production defaults generate approximately 35,328 global evaluations before local refinement.

The v9.10.2 search produced ceramic and weak-temperature candidates that persisted through 500 micrometres of reduced spatial crack growth. Its DBTT objective, however, admitted an essentially zero 300 K toughness shelf because it scored only relative transition magnitudes. Version 9.10.3 retained the same 29-dimensional physics and search domain but introduced temperature-resolved DBTT targets and a finite low-temperature shelf gate. It imposed no mechanism-dominance requirement.

The v9.11 2-D results confirm that this mechanism-neutral formulation is necessary. A DBTT-like response may arise from intrinsic first-passage kinetics with negligible retained shielding, from a sharp onset of retained MPZ shielding, or from a broad mixed response. A nonmonotonic competition can instead produce a peak-type toughness curve. The primary production DBTT is `DBTT_restart04_candidate03`; `DBTT_restart05_candidate61` is preserved as the peak class; `DBTT_restart01_candidate68` is the broad-shielding alternate; and `DBTT_restart00_candidate103` is the intrinsic low-shielding control.

The microscopic attempt frequencies remain fixed at `1e12 s^-1` and `1e11 s^-1`. Independent activation entropies provide effective-prefactor variation. Strict `H_P < H_T` reference ordering and full same-stress `G_P(sigma,T) <= G_T(sigma,T)` surface ordering remain enforced.

The active installed package remains v0.9.6. Version 9.11 integrates the independent-shape moving process zone into the 2-D adaptive FEM/CZM path. No constitutive density, stress, rate, mobile-population, jump-length, Taylor-order, or shielding cap is introduced.

Implementation records:

- [README_MPZ_V9_11_1_RESPONSE_OPTIONS.md](README_MPZ_V9_11_1_RESPONSE_OPTIONS.md)
- [README_MPZ_V9_10_3_DBTT_TARGETED.md](README_MPZ_V9_10_3_DBTT_TARGETED.md)
- [README_MPZ_V9_10_2_INDEPENDENT_SHAPES.md](README_MPZ_V9_10_2_INDEPENDENT_SHAPES.md)
- [README_MPZ_V9_10_1_SHARED_SHAPE.md](README_MPZ_V9_10_1_SHARED_SHAPE.md)
- [README_MPZ_V9_10_UNIFIED_SEARCH.md](README_MPZ_V9_10_UNIFIED_SEARCH.md)
- [README_MPZ_V9_9_BARRIER_CONTINUATION.md](README_MPZ_V9_9_BARRIER_CONTINUATION.md)
- [README_MPZ_V9_8_JOINT_OPTIMIZATION.md](README_MPZ_V9_8_JOINT_OPTIMIZATION.md)

Use `run_mpz_v9_11_response_option.sh` for named option campaigns. It materializes an isolated parameter root and prints explicit campaign start, heartbeat, final-case, and completion records.

Do not interpret the raw median `K_J` over 200–500 micrometres minus `K_init` as a material R-curve without retained-state diagnostics and an intrinsic/geometry control comparison.
