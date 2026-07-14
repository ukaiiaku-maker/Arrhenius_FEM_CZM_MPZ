# MPZ v9.3 implementation status

## Active production physics

- Crack opening: direct EXP-floor free-energy barrier and renewal clock.
- Crack-tip emission: direct EXP-floor free-energy barrier and finite source inventory.
- Peierls motion: EXP-floor barrier scaled from the active emission surface.
- Taylor obstacle escape: EXP-floor barrier scaled from the active emission surface and completed through a density-dependent correlated multi-hit renewal.
- Peierls and Taylor: sequential kinetic bottlenecks.
- Bulk FEM state: Kocks–Mecking storage/recovery with thermodynamic time-cone or Onsager admissibility.
- Moving process zone: separate mobile and retained fields, spatial capture, transport, release, recovery, shielding, blunting, and moving-frame translation.

## Explicitly non-production paths

- `legacy_additive_flow_stress`: former additive Peierls plus Taylor stress model.
- fixed independent MPZ glide/detrap barriers: retained only when emission-derived PT is explicitly disabled.
- total dislocation-density cap as a constitutive parameter: prohibited; the remote numerical ceiling is an overflow guard only.

## Validation state

- Unit and moving-process-zone regression tests pass locally.
- The v9.2 intrinsic first-passage atlas remains valid because it precedes post-emission transport.
- v9.3 Peierls–Taylor rows are analytical constitutive-screen candidates until steady-state and transient MPZ validation is completed.
