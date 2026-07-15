# Arrhenius FEM/CZM MPZ v9.10.3 DBTT calibration

The reduced spatial calibration now has validated 500-micrometre ceramic and weak-temperature/FCC-like candidates. The remaining active calibration gate is the DBTT class.

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

Version 9.10 removed the independent trap barrier: Peierls kinetics set mobile transport, forest encounters convert mobile lines to retained lines, and Taylor completion releases retained lines. Therefore `H_P < H_T` retains the intended physical interpretation.

Version 9.10.2 released the complete independent shape space:

```text
(alpha_c, n_c)
(alpha_e, n_e)
(alpha_P, n_P)
(alpha_T, n_T)
```

Every restart begins from a fresh full 29-dimensional Sobol population over common broad bounds. No previous shortlist is used as an initial population. The production defaults generate approximately 35,328 global evaluations before local refinement.

The v9.10.2 search produced ceramic and weak-temperature candidates that persisted through 500 micrometres of reduced spatial crack growth. Its DBTT objective, however, admitted an essentially zero 300 K toughness shelf because it scored only relative transition magnitudes. Version 9.10.3 retains the same 29-dimensional physics and search domain but uses the temperature-resolved DBTT design targets and a finite low-temperature shelf gate. It imposes no mechanism-dominance requirement.

The microscopic attempt frequencies remain fixed at `1e12 s^-1` and `1e11 s^-1`. Independent activation entropies provide effective-prefactor variation. Strict `H_P < H_T` reference ordering and full same-stress `G_P(sigma,T) <= G_T(sigma,T)` surface ordering remain enforced.

The active installed package remains v0.9.6 with the v9.5 spatial MPZ state. The v9.10.2 independent-shape state is used only inside the search and promotion workflows until spatial and 2-D validation is complete. No constitutive density, stress, rate, mobile-population, jump-length, Taylor-order, or shielding cap is introduced.

Implementation records:

- [README_MPZ_V9_10_3_DBTT_TARGETED.md](README_MPZ_V9_10_3_DBTT_TARGETED.md)
- [README_MPZ_V9_10_2_INDEPENDENT_SHAPES.md](README_MPZ_V9_10_2_INDEPENDENT_SHAPES.md)
- [README_MPZ_V9_10_1_SHARED_SHAPE.md](README_MPZ_V9_10_1_SHARED_SHAPE.md)
- [README_MPZ_V9_10_UNIFIED_SEARCH.md](README_MPZ_V9_10_UNIFIED_SEARCH.md)
- [README_MPZ_V9_9_BARRIER_CONTINUATION.md](README_MPZ_V9_9_BARRIER_CONTINUATION.md)
- [README_MPZ_V9_8_JOINT_OPTIMIZATION.md](README_MPZ_V9_8_JOINT_OPTIMIZATION.md)

Run `bash verify_mpz_v9_4.sh` after pulling `main`.

Do not promote the v9.10.2 DBTT candidates directly into 2-D production. Use the v9.10.3 target-aware DBTT search, require persistence in the unchanged independent-shape spatial MPZ, and only then proceed to 2-D validation.
