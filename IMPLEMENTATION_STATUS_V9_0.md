# v9.0 implementation status

## Implemented and smoke-tested

- Frozen v8 source snapshot.
- Runtime switch between `legacy_scalar` and `moving_pz`.
- Finite recoverable source-site renewal.
- Moving slip-system-resolved mobile, retained, and slip fields.
- Activated glide, trapping, detrapping, recovery, and optional annihilation.
- Conservative advection and crack-advance translation.
- Direct elastic `K_sh` and slip-based blunting.
- Complete state rollback when the crack backend vetoes a topology event.
- Conservative parent/daughter process-zone split at branch birth.
- Use of the same barrier and process-zone parameters in monotonic, cyclic, and dwell loading.
- Integration with the full anisotropic multifront sharp-front/CZM solver.
- Integration with the production mixed-mode v8 wrapper.
- Restored the shared intact-FEM helper layer required by the existing stateful local-peridynamics initiation model.
- Extended step and fatigue diagnostics.
- Legacy cap/ablation audit.
- Reduced repeated-growth fitting objective.
- Low/transition/high-temperature FEM/CZM validation runner.
- Unit tests, 1-D smoke, fatigue smoke, dwell smoke, anisotropic 2-D CZM smoke, branch-enabled anisotropic entry smoke, mixed-mode right-censored smoke, and stateful-PD core tests.
- Final active-package regression result: 139 tests and 16 subtests passed.

## Not yet completed

- Production optimization of the four material classes.
- Independent physical calibration of source density, shielding geometry factors, and transport/recovery barriers.
- Full convergence studies over `dK`, process-zone bin count, crack increment, FEM refinement, time step, and cycle-block size.
- Long-extension FEM/CZM validation of fitted classes.
- Full branch-enabled production validation.
- Full fatigue and dwell maps using fitted rows.
- Publication-level reruns and manuscript replacement figures.

## Scientific constraints

The legacy first-passage targets are retained as continuity references, not accepted as sufficient physical truth. A successful fit must also satisfy repeated-growth and state constraints. The `mpz_four_class_initial_guesses.csv` file is deliberately marked `INITIAL_GUESS_NOT_CALIBRATED`.
