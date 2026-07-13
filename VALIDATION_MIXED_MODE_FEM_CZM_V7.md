# Validation record: mixed-mode FEM/CZM v7

## Completed checks

- Python compilation of the mechanics, calibrator, campaign runner, plotter, and verifier.
- Bash syntax validation of the verification and campaign scripts.
- Thirteen unit/regression tests pass.
- The reported production-backend response matrix is included in the regression tests.
- The +30° basis solution retains its negative opening coefficient and round-trips through the full-circle coordinate.
- The phase reliability flag is independent of crystallographic directional-metric reliability.
- A synthetic backend reproduces the reported v6 negative-branch failure and the direct exact-backend root solver recovers the target.
- Barrier fingerprints, censoring states, and event-phase mismatch labels remain present.

## Not claimed

The full adaptive-CZM calibration cannot be reproduced in this execution environment because the exact active project/backend version is on the user's workstation. The three-angle v7 preflight is therefore the required integration test.

## Release criterion

Do not launch the nine-angle or crack-extension campaign unless all requested preflight calibration rows have `phase_converged=True` and the selected physical results have `event_phase_control_converged=True`.
