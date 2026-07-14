# Arrhenius FEM/CZM MPZ v9.10 unified broad search

The active calibration gate now uses one mobile/retained state model in both the analytical search and spatial crack-growth validation:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.6 uncapped, detailed-balance Peierls–Taylor constitutive audit
3. v9.7 independent-entropy calibration audit
4. v9.8/v9.8.1 joint response optimization audit
5. v9.9 barrier-scaling and spatial-mismatch audit
6. [v9.10 unified Peierls-transport/Taylor-retention global search](README_MPZ_V9_10_UNIFIED_SEARCH.md)
7. v9.10 unified spatial MPZ promotion
8. 2-D FEM/CZM validation of candidates whose response persists during crack growth

Version 9.9 demonstrated that its analytical retention proxy and spatial MPZ were not equivalent. The proxy could create shield-effective retention by suppressing Peierls motion, whereas the spatial state emitted mobile dislocations and required a separate legacy trap barrier to create retained lines. This caused the analytical DBTT R-curve to disappear during spatial crack growth and prevented the weak-temperature/FCC-like search from combining fast Peierls motion with substantial retained shielding.

Version 9.10 removes the independent trap barrier. Peierls kinetics set the mobile transport velocity, forest encounters convert mobile lines to retained lines, and Taylor completion releases retained lines. Therefore `H_P < H_T` has the intended physical interpretation: FCC-like candidates can combine a low Peierls transport barrier with a larger Taylor obstacle barrier, while ceramic-like candidates may retain slow Peierls mobility.

Every v9.10 optimization begins from a full 25-dimensional Sobol population over common broad bounds. It does not use the v9.8.1 or v9.9 shortlist as an initial population or class-specific down-selection. Cleavage, emission, Peierls, Taylor, finite source inventory, encounter efficiency, recovery, development length, and blunting are optimized together. The microscopic attempt frequencies remain fixed at `1e12 s^-1` and `1e11 s^-1`; independent activation entropies provide the effective-prefactor variation without reopening the prior eleven-decade prefactor degeneracy.

The active installed package remains v0.9.6 with the v9.5 spatial MPZ state. The v9.10 unified state is used only inside the new search and promotion workflows until spatial and 2-D validation is complete. No constitutive density, stress, rate, mobile-population, jump-length, Taylor-order, or shielding cap is introduced.

Implementation records:

- [README_MPZ_V9_10_UNIFIED_SEARCH.md](README_MPZ_V9_10_UNIFIED_SEARCH.md)
- [README_MPZ_V9_9_BARRIER_CONTINUATION.md](README_MPZ_V9_9_BARRIER_CONTINUATION.md)
- [README_MPZ_V9_8_JOINT_OPTIMIZATION.md](README_MPZ_V9_8_JOINT_OPTIMIZATION.md)
- [README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md](README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md)
- [CHANGELOG_MPZ_V9_6.md](CHANGELOG_MPZ_V9_6.md)
- [CHANGELOG_MPZ_V9_5.md](CHANGELOG_MPZ_V9_5.md)

Run `bash verify_mpz_v9_4.sh` after pulling `main`.

Do not resume the v9.9 continuation as the active calibration workflow and do not promote its DBTT or weak-temperature candidates into 2-D production. Their outputs remain useful evidence for the model mismatch. Use the v9.10 broad global search, then require persistence in the v9.10 spatial MPZ before 2-D validation.
