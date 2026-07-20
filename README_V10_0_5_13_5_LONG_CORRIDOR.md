# v10.0.5.13.5 long-corridor startup repair

## Failure reproduced

The 20 um DBTT smoke corridor started successfully, but the 100 um production
matrix failed before the first mechanics solve. All inherited multi-cloud
candidate layouts had minimum initial triangle quality below the unchanged
production floor of 0.035.

The detailed candidate audit identified near-duplicate nodes from overlapping
330 um radial-ring clouds. The old rounded-coordinate binning could leave
points 0.056987 um apart when they fell on opposite sides of a bin boundary.
Delaunay then formed a sliver triangle with edge lengths approximately
0.056987, 3.863520, and 3.876599 um.

## Repair

v10.0.5.13.5 replaces only the multi-cloud point deduplication operation:

- use Euclidean-radius neighbor detection with deterministic union-find;
- preserve exact corridor-center nodes;
- preserve exact specimen-boundary nodes;
- retain the same radial-ring placement, 330 um physical radius, Delaunay
  triangulation, 0.035 quality floor, and all FEM/CZM/MPZ physics.

The launcher also rejects a full campaign whose `OUTROOT` still appears to be a
20 um smoke-test directory.

## Real 100 um regression

The CI regression constructs the actual production startup mesh using:

- 100 um committed crack extension plus 10 um guard;
- 2.5 um requested tip spacing;
- 1.15 grading ratio;
- 330 um physical refinement radius;
- 5 um physical crack increment;
- 0.035 minimum initial triangle quality.

The validated selected mesh has:

- 3 corridor centers;
- 2,654 nodes;
- 5,186 triangles;
- minimum initial triangle quality = 0.0521124475;
- no orphan nodes;
- no nonpositive elements;
- maximum sampled h_tip/da = 1.00120884, recorded as a resolution warning only.

The focused long-corridor test and all inherited barrier, tip-only, adaptive-CZM,
refinement-support, same-mesh energy, and mechanics contracts pass.
