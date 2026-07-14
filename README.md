# Arrhenius FEM/CZM MPZ v9.5

The active workflow has three reduced front-end stages before 2-D transient fracture calculations:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.4 signed detailed-balance emission-derived Peierls–Taylor search, launched through `run_mpz_v9_3_peierls_taylor_search.sh` for command-line compatibility
3. v9.5 spatial local-density developed-state and branch-continuation audit

Version 9.5 retains the common emission-derived EXP-floor Peierls–Taylor chain and exact zero-stress detailed balance. It corrects the moving process zone so forest density and effective stress are evaluated locally in each spatial bin rather than from one global retained count and one uniform tip stress.

The v9.4 developed-state run is preserved as a diagnostic: all candidates stayed close to the virgin low-density branch because the global density conversion never activated the cooperative Taylor regime. Version 9.5 tests both virgin and physically seeded developed branches before another parameter optimization.

Implementation and change records:

- [CHANGELOG_MPZ_V9_5.md](CHANGELOG_MPZ_V9_5.md)
- [README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md](README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md)
- [CHANGELOG_MPZ_V9_4.md](CHANGELOG_MPZ_V9_4.md)
- [README_MPZ_V9_3_PEIERLS_TAYLOR_SEARCH.md](README_MPZ_V9_3_PEIERLS_TAYLOR_SEARCH.md)
- [IMPLEMENTATION_STATUS_V9_3.md](IMPLEMENTATION_STATUS_V9_3.md)
- [README_MPZ_V9_1_THREE_CLASS_TUNING.md](README_MPZ_V9_1_THREE_CLASS_TUNING.md)
- [README_MPZ_V9_0.md](README_MPZ_V9_0.md)

Run `bash verify_mpz_v9_4.sh` after pulling the current `main` branch; the retained filename now verifies package version 0.9.5. The v9.4 developed-state shortlist must not be treated as a calibrated parameter set.

The full anisotropic, mixed-mode, multifront, branching, coalescence, fatigue, dwell, checkpoint, snapshot, and adaptive FEM/CZM capabilities remain present.
