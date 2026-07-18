# v10.0.5.1 reduced 2-D slip-trace reporting

## Scope

v10.0.5.1 is a reporting-only point release layered on the v10.0.5 parallel crack-opening/emission coupling. It does not change the FEM solve, cohesive lifecycle, Arrhenius barriers, source inventory, Peierls/Taylor transport, back stress, shielding, blunting, or material manifests.

The two plastic channels in the current plane-strain model are reported as **reduced 2-D slip-trace channels**. They are not described as a complete three-dimensional BCC slip-system model.

## Certification meaning

Implementation certification requires that:

- the v10.0.5 tensor-coupling audit passed;
- per-channel drive diagnostics exist and are finite;
- the directional factor is applied only inside the emission hazard;
- the fit-derived shielding cap remains inactive;
- normalization does not modify any v10.0.5 source output.

Implementation certification does **not** require a nonzero emitted population. A zero-emission result is valid when the current material parameters and loading produce negligible emission hazard.

## New outputs

For a completed v10.0.5 directory, the normalizer writes:

- `slip_trace_reporting_v10_0_5_1.json` — reporting audit and source-file hashes;
- `mode_i_v10_0_5_1_results.json` — normalized result metadata;
- `slip_trace_channels_v10_0_5_1.csv` — long-form per-channel drive, stress, hazard, and emission increments.

Legacy `slip_system_*` fields in the v10.0.5 source outputs are preserved for compatibility. The normalized files use `slip_trace_channel_*` terminology.

## Existing output normalization

```bash
CONDA_ENV=arrhenius-fem-czm \
OUTROOT=runs/v10_0_5_parallel_weakT_700K_theta45_5um_v1 \
bash run_v10_0_5_1_normalize_existing_output.sh
```

This command launches no FEM solve and does not recompute physics.
