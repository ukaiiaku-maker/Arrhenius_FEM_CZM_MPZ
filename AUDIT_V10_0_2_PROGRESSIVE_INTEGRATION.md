# v10.0.2 progressive integration audit

## Status

The v10.0.2 branch is frozen at `922d13a`. Its isolated constitutive and lifecycle tests pass, but the full 2-D integration is not valid and must not be used for physics conclusions or penalty convergence.

## Evidence from the failed 700 K weakT run

1. The run reached 5 um physical crack extension, then raised `progressive run returned without activating the dedicated trial loop`.
2. `N_em` reached about 94 even though the weakT PF manifest has only about 2.4388 source sites per system and two systems, for a total finite source capacity of about 4.8776. The campaign finite-source state therefore was not controlling the executed emission path.
3. The legacy summary reported `advances=0` while `a_tip` moved from 0.500 to 0.505 mm.
4. The FEM startup printed `crystal theta=0.0 deg`, while the v9.11 directional summary reported `crystal_theta_deg=45.0`.
5. The v9.11 summary hard-coded `front_state_model_detail=moving_pz_v911_independent_shapes_2d_profile` even for a v10 campaign request.
6. `B_final` was reported as null although the accepted-step log printed a residual B of about 0.028 after the event.

## Root cause

`build_progressive_run_2d_v1002` compiles a transformed copy of `sharp_front.run_2d` with

```python
namespace = dict(original_run_2d.__globals__)
```

The v10.0.2 entry point performs that transformation before `mixed_mode_first_passage_v9_11.main` installs its live engine factory, J-integral wrapper, mixed boundary solver, process-zone profile wrapper, and bulk-plasticity wrapper. The transformed function therefore retains stale function references. Later monkeypatches of the `sharp_front`, `fem`, `j_integral`, and `plasticity` modules do not update the copied namespace.

Consequently, the full run can execute legacy v9.11 state and mechanics paths while the outer v10 wrapper assumes that the v10.0.2 lifecycle is active.

## Additional compatibility defects

- `CalibratedTipEngineMixin.predict_clock_increment` and `.step` dispatch `step_drives` only for `state_model == "moving_pz"`; `kinetic_campaign_czm` falls through to the legacy scalar pathway.
- Several `sharp_front.run_2d` diagnostics/export gates recognize only `moving_pz`, excluding the campaign state.
- The non-deflecting summary initializes `n_primary_adv=0` and updates it only inside the `if deflect:` block.
- The direct Mode-I wrapper and FEM parser can carry different crystal orientation defaults unless `--crystal-theta-deg` is explicit.
- Existing tests validate transforms and isolated controllers but do not verify live binding capture across the complete wrapper stack.

## v10.0.3 repair requirements

1. Transform `run_2d` lazily, after v9.11 has installed the live engine/J/mechanics/plasticity wrappers.
2. Fail before the FEM loop unless the transformed namespace captures the active v10 engine factory and all required live wrappers.
3. Use a v10-specific mixed-mode dispatch mixin that recognizes `kinetic_campaign_czm` and calls `predict_clock_increment_drives` / `step_drives` directly.
4. Preserve anisotropic FEM/J while prescribing the straight single-front Mode-I checkpoint path.
5. Treat `kinetic_campaign_czm` as a moving process-zone state for diagnostics and state export.
6. Correct non-deflecting advance accounting and cross-check summary extension against committed CZM geometry.
7. Require one explicit crystal orientation in both FEM and directional contexts.
8. Add an end-to-end binding test that would fail if a transformed function captures the stale base engine factory.
9. Do not authorize penalty convergence, longer growth, or temperature sweeps until the repaired one-segment gate passes all runtime audits.
