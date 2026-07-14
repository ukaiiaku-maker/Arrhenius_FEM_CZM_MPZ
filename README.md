# Arrhenius FEM/CZM MPZ v9.6

The active workflow now returns to a broad analytical search before another developed-state or 2-D fracture calculation:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.6 uncapped, detailed-balance, emission-derived Peierls–Taylor constitutive audit
3. v9.6 broad intrinsic × Peierls–Taylor DBTT-capacity map
4. spatial moving-process-zone development and FEM/CZM validation only after the first two v9.6 gates pass

Version 9.6 removes the exploratory Peierls–Taylor caps and algebraic saturation functions that generated an artificial broad flow-stress plateau and late stress upturn. The production closure retains the emission-derived EXP-floor barriers, natural temperature and loading-rate dependence, exact zero-stress detailed balance, natural forest spacing, and an unbounded gamma waiting-time Taylor completion.

The prior first-passage `N_sat` and shielding coefficients are retained only as analytical benchmark coordinates for the developed DBTT state. They are not used as production caps. The broad map evaluates the complete refined atlas plus the exact historical four-class references and does not require a common Peierls–Taylor closure before developed DBTT capacity is assessed.

Implementation and change records:

- [CHANGELOG_MPZ_V9_6.md](CHANGELOG_MPZ_V9_6.md)
- [CHANGELOG_MPZ_V9_5.md](CHANGELOG_MPZ_V9_5.md)
- [README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md](README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md)
- [CHANGELOG_MPZ_V9_4.md](CHANGELOG_MPZ_V9_4.md)
- [README_MPZ_V9_3_PEIERLS_TAYLOR_SEARCH.md](README_MPZ_V9_3_PEIERLS_TAYLOR_SEARCH.md)
- [IMPLEMENTATION_STATUS_V9_3.md](IMPLEMENTATION_STATUS_V9_3.md)

Run `bash verify_mpz_v9_4.sh` after pulling `main`; the retained verifier filename now checks package version 0.9.6, the uncapped PT model, the spatial v9.5 MPZ state, exact detailed balance, and absence of active constitutive caps.

Do not resume the v9.5 continuation or developed-state optimization until the v9.6 uncapped constitutive audit and broad DBTT-capacity map have been reviewed.

The full anisotropic, mixed-mode, multifront, branching, coalescence, fatigue, dwell, checkpoint, snapshot, and adaptive FEM/CZM capabilities remain present.
