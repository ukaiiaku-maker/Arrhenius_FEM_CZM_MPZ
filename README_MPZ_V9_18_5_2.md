# MPZ v9.18.5.2 — compact corridor mesh

Branch: `v9.18.5.2-compact-corridor-mesh`

## Failure addressed

The v9.18.5 and v9.18.5.1 ceramic diagnostics exited before the first physical
FEM event.  The multi-center corridor combines several radial point clouds
before Delaunay triangulation.  Delaunay/Qhull may omit duplicate or
numerically redundant input points from the returned simplices.  If those
unused points remain in `TriMesh.nodes`, their displacement DOFs have zero
stiffness and they appear as isolated one-node connected components.

## Correction

v9.18.5.2:

1. generates the same requested pre-refined Mode-I corridor;
2. finds the nodes actually referenced by bulk triangles;
3. removes unused input nodes;
4. remaps triangle connectivity;
5. rebuilds element areas, gradients, B matrices, and local/global mesh scales;
6. verifies zero orphan nodes and finite positive areas;
7. requires the initial triangle-quality floor before the first FEM assembly.

No constitutive, barrier, hazard, cohesive, MPZ, wake, source-refresh, or
shielding parameter is changed.

## New audit

Each case writes:

`compact_corridor_mesh_v91852.json`

The file records input/compacted node counts, removed node IDs, orphan count,
initial minimum triangle quality, and initial tip resolution.  It is also
written when startup raises an exception.

The sweep shell prints the tail of every per-case matrix log automatically when
a subprocess fails, so future startup failures are no longer hidden behind the
campaign summary.

## Required diagnostic

Repeat only the 700 K ceramic 60 um gate.  Acceptance requires:

- nonzero FEM/event output;
- `orphan_node_count_after_compaction = 0`;
- initial and event triangle-quality gates pass;
- no singular matrix or non-finite displacement;
- 60 um committed extension;
- immediate post-target exit;
- all v9.16–v9.18 wake/commit gates pass.
