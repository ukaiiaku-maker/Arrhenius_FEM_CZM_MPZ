# Changelog — v9.0

- Added `moving_pz` front-state model without removing the legacy scalar model.
- Added finite-site source kinetics and removed per-step emission capping from the new model.
- Added a conservative moving one-dimensional process-zone grid for every crack front.
- Added direct spatial `K_sh` evaluation and slip-based blunting.
- Disabled legacy fitted backstress, scalar saturation, stored-energy barrier lowering, and default stress capping in `moving_pz` mode.
- Preserved anisotropy, mixed-mode directional factors, branching, multifront evolution, coalescence, adaptive CZM insertion, cyclic mechanics, and existing diagnostics.
- Made branch-state splitting and geometry-veto rollback conserve the complete process-zone state.
- Routed fatigue through the same front material parameters used for monotonic and dwell loading.
- Added cap/ablation audit, repeated-growth calibration, fatigue, dwell, and small FEM/CZM validation runners.
- Added MPZ diagnostics to monotonic and fatigue histories.
- Fixed mixed-mode v8 summary writing for right-censored runs with unavailable probe values.
- Restored `arrhenius_fracture/sn_intact_fem.py`, which was referenced but absent from the supplied archive, by factoring the existing intact-FEM algorithms into a shared module; this preserves the stateful local-peridynamics initiation workflow.
- Added active-test discovery isolation from the frozen legacy tree and verified 139 tests plus 16 subtests.
