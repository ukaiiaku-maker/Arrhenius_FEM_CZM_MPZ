# Arrhenius FEM/CZM MPZ v9.9 continuation and promotion

The active workflow now uses the successful v9.8.1 basins as warm starts, reduces their absolute Peierls–Taylor barriers, and then tests the resulting candidates during finite crack growth:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.6 uncapped, detailed-balance Peierls–Taylor constitutive audit
3. v9.7 independent-entropy calibration audit
4. [v9.8/v9.8.1 joint cleavage–emission–Peierls–Taylor optimization](README_MPZ_V9_8_JOINT_OPTIMIZATION.md)
5. [v9.9 common-barrier scaling and spatial MPZ promotion](README_MPZ_V9_9_BARRIER_CONTINUATION.md)
6. 2-D FEM/CZM validation of the leading spatial candidates

Version 9.8.1 produced numerically strong ceramic and weak-temperature basins and a useful DBTT trend, but the selected absolute Peierls and Taylor barriers were large. Version 9.9 therefore scales `H_P` and `H_T` together, preserving their ratio and `H_P <= H_T`, while locally reoptimizing activation entropy, effective prefactors, Taylor correlation, finite source inventory, recovery, development length, and blunting.

The only class-specific mobility requirement is that the weak-temperature/FCC-like candidate must have a Peierls traverse number of at least one over the loading time. Slow Peierls motion is allowed for ceramic-like candidates, and no mechanism-dominance rule is imposed on DBTT candidates. DBTT acceptance is based on the low-temperature shelf, high-temperature branch, and persistence of the transition during crack extension rather than exact agreement with the provisional target points.

The active installed package remains v0.9.6 with the v9.5 spatial MPZ state. The v9.9 independent-entropy spatial adapter is used only inside the promotion workflow and is not activated globally. The prior first-passage `N_sat` and shielding coefficients remain historical benchmark coordinates, not production caps or optimizer parameters.

Implementation and change records:

- [README_MPZ_V9_9_BARRIER_CONTINUATION.md](README_MPZ_V9_9_BARRIER_CONTINUATION.md)
- [README_MPZ_V9_8_JOINT_OPTIMIZATION.md](README_MPZ_V9_8_JOINT_OPTIMIZATION.md)
- [README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md](README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md)
- [CHANGELOG_MPZ_V9_6.md](CHANGELOG_MPZ_V9_6.md)
- [CHANGELOG_MPZ_V9_5.md](CHANGELOG_MPZ_V9_5.md)
- [README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md](README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md)

Run `bash verify_mpz_v9_4.sh` after pulling `main`. The retained verifier filename now checks the v9.9 continuation helpers and independent-entropy promotion state in addition to the earlier production and calibration models.

Do not promote the raw v9.8.1 parameter values directly into 2-D production. Use the v9.9 continuation manifest, pass the selected candidates through the spatial MPZ gate, and then validate the retained ceramic, weak-temperature, and DBTT candidates in the 2-D FEM/CZM solver.
