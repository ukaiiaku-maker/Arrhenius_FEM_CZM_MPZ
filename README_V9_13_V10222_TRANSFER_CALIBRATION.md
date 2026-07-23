# v9.13 transfer calibration to the v10.2.22 top-five 2-D campaign

This update calibrates the shared 1-D transfer geometry to the
`v10_2_22_top5_dbtt_physical_width_50um_theta45_crn3621_v1` archive. The five
candidate parameterizations are immutable: no cleavage, emission, Peierls,
Taylor, source-density, recovery, or blunting value is optimized or replaced.

## Shared transfer quantities

The archive fixes the following 1-D/2-D transfer assumptions:

| Quantity | Transferred value or law |
| --- | ---: |
| Activation-to-line conversion | `(0.9124087591240877, 0.9124087591240877)` |
| Encounter efficiency | `9.160246308716648` |
| Taylor local-stress multiplier limit | `20` |
| Channel-resolved spatial transport scale | `0` |
| Moving-tip convection increment | `0.1 µm`, or `moving_tip_cfl = 0.16` |
| Minimum front width | one Burgers vector, `2.74e-10 m` |
| Maximum front width | `50 µm` |
| Direct analytical 1-D shielding | disabled for this transfer |
| Emission projection | two-channel, piecewise constant in crack extension |

The zero spatial-transport scale does **not** disable Peierls kinetics or
forest encounters. It separates encounter storage from channel-resolved
advection/escape. This is required by the archived `validated_scalar`
transport mode: the reported anisotropic transport velocities are zero and
cumulative escape is negligible. Retaining full 1-D Peierls advection flushes
the fast-Peierls 0403 candidate out of the MPZ and destroys its 600--800 K
peak.

The extension-resolved emission projection is extracted only from the primary
candidate's plastic-free 300 K trajectory. It is held fixed for every other
candidate and temperature. A projection is held constant during one accepted
crack event and recomputed after the event, matching the 2-D mechanics/update
order.

The analytical 1-D shielding projection is disabled because it is not the
mechanically measured 2-D kernel. In the supplied campaign the measured signed
shielding is below approximately `0.05 MPa sqrt(m)`, while the observed
resistance changes are tens of `MPa sqrt(m)`. Blunting, backstress, source
width, and cleavage kinetics therefore remain the active transfer mechanisms.

## Calibration and holdout validation

Run:

```bash
PYTHONPATH=. python -u scripts/calibrate_v913_to_v10222_top5.py \
  --archive /path/to/v10_2_22_top5_dbtt_physical_width_50um_theta45_crn3621_v1.zip \
  --base-physics mpz_v9_13_persistent_sites_common_physics.json \
  --base-registry candidates/v9_13_persistent_sites_top5_registry.csv \
  --out runs/v9_13_v10222_transfer_calibration_v1
```

The command writes:

* `v9_13_v10222_transfer_common_physics.json`;
* `v9_13_v10222_fixed_candidate_registry.csv`;
* `v9_13_v10222_replay_validation.csv`; and
* `v9_13_v10222_transfer_manifest.json`.

The source archive contains 50 cases. One case—the primary candidate at
300 K—is used only to extract the plastic-free geometry projection. The other
49 histories are holdouts. For the 40 holdouts with maximum backstress above
`0.05 GPa`, the current replay results are:

| Observable | Median absolute relative error | 90th-percentile error | Worst active error |
| --- | ---: | ---: | ---: |
| Maximum backstress | 0.0756% | 1.372% | 1.774% |
| Minimum front width | 0.0941% | 1.355% | 1.751% |
| Maximum tip radius | 0.0764% | 1.210% | 1.445% |

The replay uses each archived case's accepted scalar `K`, time increment, and
crack advance. It validates the 1-D constitutive transfer and moving-frame
state update, but it is not itself an autonomous prediction of the 2-D
mechanical loading path.

That remaining step is implemented in
[`README_V9_13_AUTONOMOUS_RCURVE_CALIBRATION.md`](README_V9_13_AUTONOMOUS_RCURVE_CALIBRATION.md).
The autonomous driver closes the loop with the common stochastic event clock,
the crack-geometry loading map, and a shared hazard-coupled moving-tip
operator. Its acceptance tables compare the predicted and archived
event-resolved \(K(\Delta a)\) curves directly. The prescribed ten-segment
protocol remains a fast constitutive screen; it is not used as the autonomous
R-curve objective.
