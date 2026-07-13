# Arrhenius FEM/CZM migration branch

## Purpose

This branch preserves the complete existing hazard/fracture/fatigue model and replaces only the crack-geometry realization behind a backend interface.

The retained production physics includes:

- `FractureBarrier`, including direct EXP-floor free-energy surfaces;
- `FrontEngine` emission, back-stress, blunting, stored-energy coupling, cooperative cleavage renewal, wake retention, and recovery;
- Peierls/Taylor sequential transport and escape kinetics;
- fatigue cycle quadrature and adaptive cycle blocks;
- spatial process-zone state (`emit`, `store`, `mobile`, `escape` fields);
- anisotropic directional competition;
- branch-specific first-passage clocks;
- multi-front inventory, unresolved branch clusters, local/cluster J handoff, and stagnant-branch bookkeeping.

The migration is therefore **not** a reduced cohesive-zone model that must later regain the mature physics.

## Backend boundary

The new interface is in:

- `arrhenius_fracture/crack_backend.py`

Available backends:

### `sharp_wake`

Compatibility backend reproducing the legacy geometry update:

1. FrontEngine completes an Arrhenius renewal.
2. The explicit crack path advances.
3. Bulk stiffness is removed around the new path segment.

This backend is retained for regression comparison.

### `edge_split_czm`

First working discrete-CZM migration backend:

1. The existing directional hazard machinery selects the physical propagation direction.
2. The event stepper limits topology updates to one completed renewal per geometry solve.
3. The backend selects the best existing mesh edge in the selected direction.
4. Endpoint node stars are split into the two crack sides.
5. Coincident crack-surface nodes are created.
6. A zero-thickness cohesive interface record is inserted.
7. The interface broken fraction is committed by the Arrhenius renewal event.
8. Bulk element **count and ordering do not change**, so Gauss-point plasticity, density, and process-zone histories remain aligned without field projection.

### `adaptive_czm`

Angle-faithful local r-adaptive CZM backend:

1. The directional hazard machinery supplies the physical propagation direction.
2. The requested physical increment `da_phys` defines the exact target point on that ray.
3. A geometric one-ring neighbor of the active tip is selected for relocation.
4. The candidate relocation is accepted only if every incident triangle preserves orientation, exceeds the minimum area-ratio constraint, and exceeds the minimum triangle-quality constraint.
5. The nodal displacement at the relocated point is interpolated from the pre-update FEM field.
6. The existing topology split then creates the two crack faces and inserts the zero-thickness Arrhenius cohesive interface.

Element connectivity, element count, and element ordering remain unchanged during the local r-adaptation, so Gauss-point plasticity, density, and process-zone arrays retain their indices. This removes the mesh-angle steering seen in `edge_split_czm` without repeatedly projecting the full process-zone history to a new mesh.

For very long or highly tortuous growth, accumulated element-quality degradation remains possible; the backend records local area ratio, triangle quality, and node motion for each event. Full patch retriangulation remains a possible later fallback, not a prerequisite for the oriented temperature sweep.

## Pure Arrhenius cohesive formulation

The cohesive representation does not introduce an independent failure criterion.

The existing EXP-floor barrier remains authoritative:

\[
G_c^*(\sigma,T)
=
G_{\rm floor}(T)
+
\left[G_0(T)-G_{\rm floor}(T)\right]
\exp\left[-a\left(\sigma/\sigma_c(T)\right)^n\right].
\]

The existing raw and cooperative rates remain:

\[
\lambda_c^{\rm raw}
=\nu_c\exp[-G_c^*/k_BT],
\]

\[
\lambda_c^{\rm renew}
=
P(m,\lambda_c^{\rm raw}\tau_c)/\tau_c.
\]

A completed renewal is the only fracture event. The CZM backend only changes how the created crack surface is represented mechanically.

The initial production setting is discrete link failure:

- intact before the event;
- one completed renewal;
- newly inserted interface has `damage=1` by default.

This is the closest spatial translation of the current sharp-front renewal model and the discrete Arrhenius CZM proof-of-concept. The data structure also supports partially active interfaces (`damage < 1`) for a later progressive hazard-driven separation law, but no empirical traction or opening failure threshold is added.

## Why topology-only edge splitting is the first backend

The critical migration risk is history transfer in the high-gradient process zone. Repeated remeshing can diffuse or perturb:

- plastic strain;
- dislocation density;
- retained/mobile/escaped process-zone populations;
- residual stresses derived from those fields.

Because hazards are exponentially sensitive to their driving state, small transfer errors can produce large kinetic errors.

The edge-split backend therefore changes only nodal topology. Bulk elements and their integration-point indices remain fixed. This gives a robust full-physics migration baseline before local patch retriangulation is introduced.

## Event stepping

For CZM backends, monotonic adaptive event stepping is enabled automatically. `FrontEngine` now honors `max_advances_per_step`:

- at most one topology event is consumed per accepted CZM geometry solve;
- additional completed renewal count remains in the cleavage clock `B`;
- no hazard or event count is discarded.

Fatigue retains the existing adaptive cycle-block controller. The same global physical cycle count is used by all active fronts.

## Key new files

- `arrhenius_fracture/cohesive.py` — cohesive element/network state and interface traction assembly.
- `arrhenius_fracture/crack_backend.py` — backend interface, legacy compatibility backend, topology-preserving CZM backend.
- `arrhenius_fracture/material_state.py` — first-class state container for the later persistent state-carrier migration.
- `arrhenius_fracture/local_patch.py` — patch-selection boundary for the future local-remeshing backend.
- `tests/test_czm_backend.py` — topology and cohesive assembly tests.
- `run_fem_czm_fullphysics_smoke.sh` — complete fatigue/plasticity/CZM smoke campaign using the existing v8 comparison driver.

## Running

Legacy regression:

```bash
python -m arrhenius_fracture.sharp_front \
  --mode 2d \
  --crack-backend sharp_wake \
  ...existing arguments...
```

Discrete Arrhenius CZM:

```bash
python -m arrhenius_fracture.sharp_front \
  --mode 2d \
  --crack-backend edge_split_czm \
  --cleave-barrier-kind exp_floor \
  ...existing arguments...
```

Full-physics smoke:

```bash
bash run_fem_czm_fullphysics_smoke.sh
```

Existing production campaign with the CZM backend:

```bash
CRACK_BACKEND=adaptive_czm \
CASE_FILTER="plastic_shielded_case64_M1" \
KLIST_OVERRIDE="6.0 7.0 8.0" \
bash run_v8_material_response_production_2d.sh
```

## Current limitations

1. `edge_split_czm` still follows the best existing edge and is retained as a regression backend. Use `adaptive_czm` for anisotropic path studies.
2. `adaptive_czm` uses exact-angle local r-adaptation, not full patch retriangulation. Long highly tortuous runs should monitor the per-event `min_area_ratio` and `min_triangle_quality` diagnostics for cumulative mesh degradation.
3. The initial cohesive event is abrupt (`--czm-event-damage 1`). Progressive hazard-driven cohesive opening is structurally supported but not yet coupled to a pre-inserted active tip segment.
4. Branch bookkeeping is retained, but daughter-tip topology and long multi-branch growth still require dedicated regression testing before large production campaigns.
5. The initial notch remains represented by the legacy damaged notch band; newly propagated segments use topological CZM surfaces.

These limitations are explicit so that the first comparisons test the architecture cleanly rather than mixing unverified remeshing, field projection, and new constitutive assumptions at once.
