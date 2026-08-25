# One-dimensional V2 final model architecture

## Outcome

The V2 implementation has exact controlled common-kernel parity and qualified candidate-independent PF and FEM/CZM mechanics/source-drive maps. The natural predictive lane is **partially qualified**, not universally predictive: PF Peak closes, while DBTT and absolute FEM/CZM onset do not.

## Runtime composition

1. `SharedBarrierHazardCore` owns the common barrier/rate equations.
2. Backend state adapters/factories construct production-source PF unified-MPZ and FEM/CZM moving-tip states.
3. `ProviderMechanicsMap` supplies native PF or native/qualified FEM/CZM coefficients with fail-closed extension bounds.
4. `SourceDriveMap` supplies exact tensor-probe normalized drive factors on extension/radius grids. The final grid spans 0–1000 µm extension and 1–100 µm radius; no clipping or extrapolation is permitted.
5. `run_zero_d_predictive` is a versioned V9.13 zero-D state/lifecycle surrogate. It is not claimed to reproduce the full backend-specific production event transaction policy.

The controlled composition result remains `CONTROLLED_SOURCE_COMPOSITION_QUALIFIED`. Natural prediction is assessed separately.

No production equation, production trajectory, or canonical parameter registry was modified. New 2-D PF runs: **0**. New 2-D FEM/CZM runs: **0**.
