# MPZ v9.18.3: edge-aware exact-ray geometry recovery

Branch: `v9.18.3-edge-aware-geometry-recovery`

## Reproduced failure

The v9.18.2 ceramic 700 K, 100 um run advanced to 45 um and then repeated:

```text
GEOMETRY VETO front 0:
local_hrefine_error:degenerate elements after topology update: [2126]
-- renewal retained in B=1.000
```

The same front, target, element, and reason repeated indefinitely. Increasing
`STEPS` cannot resolve this deterministic geometry state.

## Cause

After repeated r-adaptive steering, the exact 5 um target lay on an existing
triangle edge. The incident-triangle h-refinement fallback accepted a
barycentric coordinate near zero as an interior point and split the parent into
three children. One child necessarily had zero area, so `rebuild_tri_mesh`
rejected element 2126.

This is a topology-classification error, not a kinetic, hazard, or material
parameter problem.

## Correction

v9.18.3 classifies the target before constructing topology:

1. target at an existing vertex: reuse the vertex;
2. target on an existing edge: split the edge at the exact-ray target and
   inherit integration-point histories from each parent;
3. target strictly inside a triangle: perform the three-child split only after
   an explicit child-area and quality precheck.

The complete Arrhenius fracture event remains transactional. If geometry still
fails, the mesh, cohesive network, front state, and physical renewal are rolled
back together.

A repeated-identical-veto guard aborts after 12 identical backend failures by
default. This prevents a pathological geometry state from consuming tens of
thousands of steps. The guard does not consume the physical event or alter the
cleavage action.

## Physics unchanged

The branch does not change:

- cleavage or emission barriers;
- absolute hazard integration;
- opening relaxation;
- source refresh;
- Peierls--Taylor retention;
- persistent wake transport;
- wake shielding;
- physical crack quantum.

## First diagnostic

Stop the stuck v9.18.2 process and run only ceramic through 60 um, which crosses
the previous 45 um failure location:

```bash
TEMPS="700" \
CLASSES="ceramic" \
TARGET_EXT_UM=60 \
STEPS=15000 \
OUTROOT_BASE=runs/mpz_v9_18_3_geometry_recovery_ceramic_700K_60um_v1 \
bash run_mpz_v9_18_3_persistent_plastic_wake_sweep.sh
```

Acceptance criteria:

- extension passes 45 um;
- no repeated element-2126 veto loop;
- `geometry_recovery_v9183.json` exists;
- committed-event and persistent-wake audits pass;
- output is marked complete by the v9.18.2 committed-completion handshake.

## Follow-up 100 um comparison

After the ceramic diagnostic passes, run the three classes with a fresh output
root. Do not reuse the incomplete v9.18.2 directory.
