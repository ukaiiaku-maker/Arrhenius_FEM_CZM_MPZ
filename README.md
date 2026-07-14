# Arrhenius FEM/CZM MPZ v9.8 joint optimization

The active analytical workflow now optimizes fracture and plasticity together before another developed-state or 2-D calculation:

1. [v9.2 analytical virgin-tip first-passage atlas](README_MPZ_V9_2_ANALYTIC_ATLAS.md)
2. v9.6 uncapped, detailed-balance Peierls–Taylor constitutive audit
3. v9.7 independent-entropy calibration audit
4. [v9.8 joint cleavage–emission–Peierls–Taylor response optimization](README_MPZ_V9_8_JOINT_OPTIMIZATION.md)
5. spatial moving-process-zone promotion and FEM/CZM validation of the leading diverse candidates

Version 9.6 removed the exploratory Peierls–Taylor caps and algebraic saturation functions that generated an artificial broad flow-stress plateau and late stress upturn. Version 9.7 demonstrated that Peierls and Taylor parameters are not identifiable from one standalone bulk flow-stress target. Version 9.8 therefore optimizes cleavage, emission, Peierls, Taylor, and reduced state parameters together against initial toughness, developed resistance, R-curve increment, and temperature dependence.

The v9.8 global search uses atlas-seeded differential evolution followed by bounded Nelder–Mead refinement. It imposes no mechanism-dominance rules for ceramic-like, weak-temperature, or DBTT responses. It requires only positive barriers, exact detailed balance, uncapped kinetics, and `G_P(sigma,T) <= G_T(sigma,T)` on a common resolved-stress grid.

The prior first-passage `N_sat` and shielding coefficients remain historical benchmark coordinates, not production caps or direct optimizer parameters. The complete anisotropic, mixed-mode, multifront, branching, coalescence, fatigue, dwell, checkpoint, snapshot, and adaptive FEM/CZM capabilities remain present and unchanged.

Implementation and change records:

- [README_MPZ_V9_8_JOINT_OPTIMIZATION.md](README_MPZ_V9_8_JOINT_OPTIMIZATION.md)
- [README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md](README_MPZ_V9_7_PT_ENTROPY_CALIBRATION.md)
- [CHANGELOG_MPZ_V9_6.md](CHANGELOG_MPZ_V9_6.md)
- [CHANGELOG_MPZ_V9_5.md](CHANGELOG_MPZ_V9_5.md)
- [README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md](README_MPZ_V9_4_DEVELOPED_STATE_SEARCH.md)
- [CHANGELOG_MPZ_V9_4.md](CHANGELOG_MPZ_V9_4.md)

Run `bash verify_mpz_v9_4.sh` after pulling `main`. The retained verifier filename now checks the v9.8 optimizer helpers in addition to the active v9.6 production PT model, the v9.5 spatial MPZ state, and the v9.7 entropy-decoupling model.

Do not resume the earlier v9.5 continuation, v9.6 broad DBTT map, or standalone v9.7 calibration as a production workflow. Their outputs remain useful audits and seed information for v9.8.
