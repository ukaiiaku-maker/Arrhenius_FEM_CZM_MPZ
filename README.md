# Arrhenius FEM/CZM MPZ v9.4

The active workflow has three reduced front-end stages before 2-D transient fracture calculations:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.4 signed detailed-balance emission-derived Peierls–Taylor search, launched through `run_mpz_v9_3_peierls_taylor_search.sh` for command-line compatibility
3. [v9.4 common developed moving-process-zone search](README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md)

Version 9.4 retains the common emission-derived EXP-floor Peierls–Taylor chain for the moving process zone and bulk FEM plasticity, while restoring forward-minus-reverse detailed balance. Net plastic flow is exactly zero at zero effective stress. Taylor completion remains a correlated multi-hit renewal so the independent-Poisson high-density stress turnover is not hidden with a fitted dislocation-density cap.

The developed-state stage holds the strict common Peierls–Taylor closure and intrinsic cleavage/emission surfaces fixed. It searches only active moving-process-zone state parameters and applies one shared closure to ceramic, weak-temperature, and DBTT-precursor rows.

Implementation and change records:

- [README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md](README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md)
- [CHANGELOG_MPZ_V9_4.md](CHANGELOG_MPZ_V9_4.md)
- [README_MPZ_V9_3_PEIERLS_TAYLOR_SEARCH.md](README_MPZ_V9_3_PEIERLS_TAYLOR_SEARCH.md)
- [IMPLEMENTATION_STATUS_V9_3.md](IMPLEMENTATION_STATUS_V9_3.md)
- [README_MPZ_V9_1_THREE_CLASS_TUNING.md](README_MPZ_V9_1_THREE_CLASS_TUNING.md)
- [README_MPZ_V9_0.md](README_MPZ_V9_0.md)

Run `bash verify_mpz_v9_4.sh` after pulling the current `main` branch. The completed v9.3 v2 search is retained only as an audit of the missing detailed-balance condition; its shortlist must not be used for developed-state MPZ fitting.

The full anisotropic, mixed-mode, multifront, branching, coalescence, fatigue, dwell, checkpoint, snapshot, and adaptive FEM/CZM capabilities remain present.
