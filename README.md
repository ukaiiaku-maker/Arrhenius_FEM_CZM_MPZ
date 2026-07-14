# Arrhenius FEM/CZM MPZ v9.6 + v9.7 calibration

The active workflow remains on the uncapped v9.6 production closure, but the next analytical gate is now the v9.7 Peierls–Taylor magnitude and entropy calibration:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.6 uncapped, detailed-balance Peierls–Taylor constitutive audit
3. [v9.7 reference-strength and independent-activation-entropy calibration](README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md)
4. broad intrinsic × calibrated Peierls–Taylor DBTT-capacity map
5. spatial moving-process-zone development and FEM/CZM validation only after these analytical gates pass

Version 9.6 removes the exploratory Peierls–Taylor caps and algebraic saturation functions that generated an artificial broad flow-stress plateau and late stress upturn. The production closure retains EXP-floor barriers, natural temperature and loading-rate dependence, exact zero-stress detailed balance, natural forest spacing, and an unbounded gamma waiting-time Taylor completion.

The v9.6 audit also showed that the historical `0.005/0.02` energy ratios do not provide a GPa-scale reference flow stress for the canonical emission surfaces. Version 9.7 therefore calibrates the reference-strength magnitude first and then varies Peierls and Taylor activation entropies independently. The entropy-calibration model is not activated in the MPZ or bulk FEM until its output is reviewed.

The prior first-passage `N_sat` and shielding coefficients remain analytical benchmark coordinates for the developed DBTT state, not production caps. A common Peierls–Taylor closure is not imposed before developed DBTT capacity is assessed.

Implementation and change records:

- [README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md](README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md)
- [CHANGELOG_MPZ_V9_6.md](CHANGELOG_MPZ_V9_6.md)
- [CHANGELOG_MPZ_V9_5.md](CHANGELOG_MPZ_V9_5.md)
- [README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md](README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md)
- [CHANGELOG_MPZ_V9_4.md](CHANGELOG_MPZ_V9_4.md)
- [README_MPZ_V9_3_PEIERLS_TAYLOR_SEARCH.md](README_MPZ_V9_3_PEIERLS_TAYLOR_SEARCH.md)
- [IMPLEMENTATION_STATUS_V9_3.md](IMPLEMENTATION_STATUS_V9_3.md)

Run `bash verify_mpz_v9_4.sh` after pulling `main`. The retained verifier filename checks package version 0.9.6, the active uncapped PT model, the spatial v9.5 MPZ state, exact detailed balance, absence of constitutive caps, and the v9.7 entropy-decoupling calibration model.

Do not resume the v9.5 continuation, developed-state optimization, or v9.6 broad DBTT map until the v9.7 magnitude/entropy calibration has been reviewed.

The full anisotropic, mixed-mode, multifront, branching, coalescence, fatigue, dwell, checkpoint, snapshot, and adaptive FEM/CZM capabilities remain present.
