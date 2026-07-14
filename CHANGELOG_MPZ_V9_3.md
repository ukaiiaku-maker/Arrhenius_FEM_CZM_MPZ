# MPZ v9.3 changelog

- Added a shared emission-derived EXP-floor Peierls–Taylor kinetics module.
- Made production bulk FEM plasticity use sequential Peierls and Taylor rates instead of an independent rational barrier and additive flow stress.
- Added correlated multi-hit Taylor completion to prevent the independent-Poisson high-density flow-stress downturn without a constitutive total-density cap.
- Preserved Kocks–Mecking storage/recovery and local thermodynamic admissibility.
- Made the moving process zone use the same emission-derived Peierls rate for glide and correlated Taylor rate for release of retained content.
- Retained spatial capture, recovery, transport, shielding, blunting, branching, coalescence, fatigue, dwell, restart, and FEM/CZM architecture.
- Restored Peierls/emission and Taylor/emission scale centers of 0.005 and 0.02 throughout monotonic and cyclic paths.
- Added a post-emission parameter search that rejects unresolved or high-density-softening bulk closures before transient fracture calculations.
- Relegated the former additive Peierls–Taylor flow-stress model to an explicitly named legacy ablation.
