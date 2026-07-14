# Arrhenius FEM/CZM MPZ v9.3

The active workflow now has two analytical front-end stages before transient fracture calculations:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. [v9.3 emission-derived Peierls–Taylor monotonicity search](README_MPZ_V9_3_PEIERLS_TAYLOR_SEARCH.md)

Version 9.3 makes the moving process zone and bulk FEM plasticity use the same emission-derived EXP-floor Peierls–Taylor chain. Taylor completion uses a correlated multi-hit renewal so the independent-Poisson high-density stress turnover is not hidden with a fitted dislocation-density cap.

Implementation and change records:

- [IMPLEMENTATION_STATUS_V9_3.md](IMPLEMENTATION_STATUS_V9_3.md)
- [CHANGELOG_MPZ_V9_3.md](CHANGELOG_MPZ_V9_3.md)
- [README_MPZ_V9_1_THREE_CLASS_TUNING.md](README_MPZ_V9_1_THREE_CLASS_TUNING.md)
- [README_MPZ_V9_0.md](README_MPZ_V9_0.md)

The full anisotropic, mixed-mode, multifront, branching, coalescence, fatigue, dwell, checkpoint, snapshot, and adaptive FEM/CZM capabilities remain present.
