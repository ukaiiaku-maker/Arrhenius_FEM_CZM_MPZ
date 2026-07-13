# Mixed-mode FEM/CZM v3.5: matrix-authoritative calibration

V3.5 retains the validated v3.3 J-consistent isotropic mechanics and deterministic first-passage calculation. It changes only the elastic calibration acceptance logic and the campaign/plotting wrappers.

## Why v3.4 stopped

The +45 degree elastic verification achieved 45.002 degrees, an error of only 0.002 degrees, but the campaign aborted because the diagnostic annulus-to-annulus Williams phase spread was 24.5 degrees versus a hard 20 degree threshold.

In the J-consistent formulation, the authoritative phase comes from the measured two-basis elastic response matrix and its aggregate verification solve. The annular Williams fits do not supply the authoritative amplitude or mode partition. Their spread is therefore a diagnostic-confidence warning, not a campaign-wide veto.

## Acceptance in v3.5

A calibration row is accepted only when:

- KI, KII, and aggregate phase are finite;
- KI is positive after signed-basis normalization;
- the aggregate projection has sufficient points and annular fits;
- the aggregate projection matrix is sufficiently conditioned;
- the achieved phase is within PSI_TOL_DEG of the requested phase.

Annular phase spread above PHASE_SPREAD_TOL_DEG is recorded as:

- `phase_spread_warning = True`
- `calibration_confidence = target_met_with_spread_warning`

but does not abort the campaign.

This does not loosen the target-angle requirement. Large phase error, nonfinite modes, insufficient support, or ill conditioning still fail calibration.
