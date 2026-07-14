# Arrhenius FEM/CZM MPZ v9.10.2 independent-shape search

The active calibration gate now uses one unified mobile/retained state model while allowing the four EXP-floor barriers to explore independent stress-shape coordinates:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.6 uncapped, detailed-balance Peierls–Taylor constitutive audit
3. v9.7 independent-entropy calibration audit
4. v9.8/v9.8.1 joint response optimization audit
5. v9.9 barrier-scaling and spatial-mismatch audit
6. [v9.10 unified Peierls-transport/Taylor-retention search](README_MPZ_V9_10_UNIFIED_SEARCH.md)
7. v9.10.1 one-shape-for-all-four audit
8. [v9.10.2 independent four-barrier EXP-floor search](README_MPZ_V9_10_2_INDEPENDENT_SHAPES.md)
9. v9.10.2 independent-shape spatial MPZ promotion
10. 2-D FEM/CZM validation of candidates whose response persists during crack growth

Version 9.10 removed the independent trap barrier: Peierls kinetics set mobile transport, forest encounters convert mobile lines to retained lines, and Taylor completion releases retained lines. Therefore `H_P < H_T` retains the intended physical interpretation.

The original v9.10 shape hypothesis was

```text
alpha_e = alpha_P = alpha_T
n_e     = n_P     = n_T
```

with an independent cleavage shape. The v9.10.1 audit imposed one shape on all four barriers. Its smoke search produced stable ceramic behavior and spatially persistent DBTT candidates, but no accepted weak-temperature/FCC-like candidate. Version 9.10.2 therefore releases the full shape space:

```text
(alpha_c, n_c)
(alpha_e, n_e)
(alpha_P, n_P)
(alpha_T, n_T)
```

Every v9.10.2 restart begins from a fresh full 29-dimensional Sobol population over common broad bounds. It does not use any previous shortlist as its initial population. The production defaults generate approximately 35,328 global evaluations per class before local refinement, matching the scale of the earlier approximately 35,000-point searches while covering the enlarged phase space.

The microscopic attempt frequencies remain fixed at `1e12 s^-1` and `1e11 s^-1`. Independent activation entropies provide effective-prefactor variation. Strict `H_P < H_T` reference ordering and full same-stress `G_P(sigma,T) <= G_T(sigma,T)` surface ordering remain enforced.

The active installed package remains v0.9.6 with the v9.5 spatial MPZ state. The v9.10.2 state is used only inside the search and promotion workflows until spatial and 2-D validation is complete. No constitutive density, stress, rate, mobile-population, jump-length, Taylor-order, or shielding cap is introduced.

Implementation records:

- [README_MPZ_V9_10_2_INDEPENDENT_SHAPES.md](README_MPZ_V9_10_2_INDEPENDENT_SHAPES.md)
- [README_MPZ_V9_10_1_SHARED_SHAPE.md](README_MPZ_V9_10_1_SHARED_SHAPE.md)
- [README_MPZ_V9_10_UNIFIED_SEARCH.md](README_MPZ_V9_10_UNIFIED_SEARCH.md)
- [README_MPZ_V9_9_BARRIER_CONTINUATION.md](README_MPZ_V9_9_BARRIER_CONTINUATION.md)
- [README_MPZ_V9_8_JOINT_OPTIMIZATION.md](README_MPZ_V9_8_JOINT_OPTIMIZATION.md)

Run `bash verify_mpz_v9_4.sh` after pulling `main`.

Do not promote the v9.10 or v9.10.1 smoke candidates directly into 2-D production. Use the v9.10.2 broad search, require persistence in the independent-shape unified spatial MPZ, and only then proceed to 2-D validation.
