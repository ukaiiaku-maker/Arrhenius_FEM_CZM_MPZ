# One-dimensional V2 PF mechanics-map qualification V3

Status: **PF_PRODUCTION_DISCRETE_MAP_QUALIFIED**.

The map applies the actual `SharpWakeBackend.advance` operation sequentially on
the production 36×72 graded plane-strain mesh. Each standardized straight event
is the source-defined 5 µm segment; nominal 2 µm therefore resolves to the first
complete 5 µm event and no segment is truncated. Binary killed-element damage is
accumulated, never independently reconstructed at each node.

The production cluster-domain J convention (`ell=80 µm`, exclusion radius twice
the kill radius, signed/effective aggregation) is retained. This is a
**model-native discrete map, not continuum G**. Maximum opening-scaling error is
`0.000e+00`. Bounded piecewise-linear interpolation is valid only over
0–1000 µm; extrapolation fails closed. Maximum LOO relative error is
`1.129e-01` and is retained as discrete interpolation
uncertainty rather than smoothed away.
