# MPZ v9.18.4: mechanically valid crack-front regularization

Branch: `v9.18.4-mechanically-valid-front-regularization`

## Result of the v9.18.3 700 K archive

v9.18.3 fixed the original repeated element-2126 geometry veto. Ceramic, weakT,
and DBTT all crossed the previous 45 um stall and inserted the 50 um event.
Immediately after that event, however, all three runs reported:

```text
MatrixRankWarning: Matrix is exactly singular
```

The non-finite FEM displacement then propagated into moving-process-zone
advection and terminated as:

```text
ValueError: cannot convert float NaN to integer
```

Thus the edge-aware element geometry was valid in the narrow sense of positive
triangle area, but the accepted crack-front topology was not mechanically
well-posed.

## Cause

At a target lying exactly on an existing edge, v9.18.3 split the edge and left
the newly opened endpoint exactly on that mesh edge. The bulk triangles had
positive area, but the failed cohesive interface has zero tensile/shear
stiffness. The resulting split could therefore leave an orphan crack-surface
node or a disconnected bulk component with a rigid-body mode. The next linear
solve was singular.

## Correction

Only when the requested endpoint lies numerically on an edge, v9.18.4 replaces
the edge-front endpoint with a same-length endpoint rotated by a small bounded
angle into an adjacent triangle. Candidate angles are, by default:

```text
0.01 0.025 0.05 0.1 0.2 0.35 degrees
```

Both signs are evaluated. Candidates retain the exact physical advance length,
are ranked by minimum child-element quality, and are attempted transactionally.
The requested and regularized coordinates, angular perturbation, and length
error are stored in the CZM advance metadata.

Before a geometry event is accepted, v9.18.4 requires:

- finite positive bulk elements;
- no orphan bulk node;
- every cohesive endpoint to have bulk support;
- every disconnected bulk component to possess both x and y Dirichlet
  anchoring under the actual model boundary conditions.

A failed candidate rolls back the mesh, displacement, cohesive network, front
bookkeeping, and physical renewal together. The next candidate can then be
tried without consuming the Arrhenius event.

The mechanics solve also fails immediately on non-finite displacement or
reaction, preventing a singular solve from entering MPZ transport.

## Physics unchanged

v9.18.4 does not change:

- cleavage or emission barriers;
- absolute-hazard clocks;
- source refresh;
- cohesive opening progress;
- physical crack-advance length;
- Peierls--Taylor state evolution;
- persistent-wake transfer;
- wake shielding.

The only approximation is a recorded sub-degree local geometry regularization
when an exact endpoint coincides with an edge.

## Required diagnostic

Do not reuse the failed v9.18.3 output. First run ceramic to 60 um in a new
output root. This crosses both the old 45 um geometry veto and the v9.18.3 50 um
singular-topology failure.

```bash
TEMPS="700" \
CLASSES="ceramic" \
TARGET_EXT_UM=60 \
STEPS=15000 \
OUTROOT_BASE=runs/mpz_v9_18_4_ceramic_700K_60um_v1 \
bash run_mpz_v9_18_4_persistent_plastic_wake_sweep.sh
```

Acceptance requires committed extension of 60 um, no singular-matrix warning,
no NaN transport exception, and a complete v9.18 wake audit.
