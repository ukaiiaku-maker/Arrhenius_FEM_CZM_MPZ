# One-dimensional V2 reload-separated resistance definition

## Qualified screening quantity

The reduced-model screening quantity is **R-curve propensity**, or **reload-separated resistance development**. A resistance candidate is the pre-event state of the first physical avalanche or the pre-event state after a qualified reload separating two physical avalanches.

For onset candidates `m`, the native metric is:

`DeltaK_reinit = max_m(K_onset[m]) - K_onset[1]`.

For FEM/CZM structural interpretation the separately qualified metric is:

`DeltaG_reinit = max_m(G_onset[m]) - G_onset[1]`.

`N_reinit` counts later reload-separated onset states. Target completion is right-censoring, not physical arrest.

## Excluded interpretation

Every numerical event inside one physical avalanche remains in the archive as a **model-native driving trajectory**. Its event-wise native J/KJ evolution, final KJ, and within-avalanche slope are not resistance points and are not fitted as an R-curve.

Quantitative R-curve claims are reserved for bounded direct 2-D PF cases after V2 event-transaction versus physical-avalanche grouping.
