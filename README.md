# Arrhenius FEM/CZM MPZ v9.10.1 shared-shape search

The active calibration gate now uses one mobile/retained state model and one common EXP-floor stress shape in both the analytical search and spatial crack-growth validation:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.6 uncapped, detailed-balance Peierls–Taylor constitutive audit
3. v9.7 independent-entropy calibration audit
4. v9.8/v9.8.1 joint response optimization audit
5. v9.9 barrier-scaling and spatial-mismatch audit
6. [v9.10 unified Peierls-transport/Taylor-retention global search](README_MPZ_V9_10_UNIFIED_SEARCH.md)
7. [v9.10.1 shared EXP-floor shape global search](README_MPZ_V9_10_1_SHARED_SHAPE.md)
8. v9.10.1 shared-shape spatial MPZ promotion
9. 2-D FEM/CZM validation of candidates whose response persists during crack growth

Version 9.9 demonstrated that its analytical retention proxy and spatial MPZ were not equivalent. Version 9.10 removed the independent trap barrier: Peierls kinetics set the mobile transport velocity, forest encounters convert mobile lines to retained lines, and Taylor completion releases retained lines. Therefore `H_P < H_T` has the intended physical interpretation.

The first v9.10 smoke search still used an asymmetric EXP-floor parameterization: cleavage and emission had independently optimized `alpha` and `n`, while Peierls and Taylor inherited the emission shape. Version 9.10.1 now optimizes one common pair

```text
alpha_c = alpha_e = alpha_P = alpha_T
n_c     = n_e     = n_P     = n_T
```

while retaining independent cleavage, emission, Peierls, and Taylor barrier heights and thermal coordinates. This tests whether one free-energy shape with different barrier heights is sufficient before adding four independent shape pairs.

Every v9.10.1 optimization begins from a fresh full 23-dimensional Sobol population over common broad bounds. It does not use the v9.8.1, v9.9, or v9.10 shortlist as an initial population. The microscopic attempt frequencies remain fixed at `1e12 s^-1` and `1e11 s^-1`.

The active installed package remains v0.9.6 with the v9.5 spatial MPZ state. The v9.10 unified state is used only inside the search and promotion workflows until spatial and 2-D validation is complete. No constitutive density, stress, rate, mobile-population, jump-length, Taylor-order, or shielding cap is introduced.

Implementation records:

- [README_MPZ_V9_10_1_SHARED_SHAPE.md](README_MPZ_V9_10_1_SHARED_SHAPE.md)
- [README_MPZ_V9_10_UNIFIED_SEARCH.md](README_MPZ_V9_10_UNIFIED_SEARCH.md)
- [README_MPZ_V9_9_BARRIER_CONTINUATION.md](README_MPZ_V9_9_BARRIER_CONTINUATION.md)
- [README_MPZ_V9_8_JOINT_OPTIMIZATION.md](README_MPZ_V9_8_JOINT_OPTIMIZATION.md)
- [README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md](README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md)

Run `bash verify_mpz_v9_4.sh` after pulling `main`.

Do not promote the v9.10 smoke candidates into 2-D production. Use the v9.10.1 shared-shape broad search, require persistence in the unified spatial MPZ, and only then proceed to 2-D validation. If a well-converged shared-shape search still fails, release four independent `(alpha,n)` pairs in a new full-space search rather than refining a prior shortlist.
