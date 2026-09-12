# One-dimensional V3 stateful-voiding crosswalk

## Mission boundary and provenance

This reduction starts from one-dimensional V2 commit
`5045feed8c87440676c580fe7f74989be4c4a995` on
`codex/oneD-v2-terminal-predictive-program`. The two-dimensional oracle is the
read-only V5 implementation on `ukaiiaku-maker/PF-fracture-fatigue` PR 63. Its
runtime hardening commit is
`b58997bdb18cf4e9a32c251c073115d8b405bb27`, and its successful bounded
attestation commit is `c7583ecd0a259f28ce92780833d21a358a920f45`.

No two-dimensional branch is merged into this repository. V2 fracture state,
barriers, hazard evaluation, mechanics ownership, renewal, translation,
backend veto policy, and material rows remain unchanged. V3 composes the V2
fracture state with a separate void state and reserves a third independent
slot for later fatigue state:

```text
OneDStateV3 = (ReducedState V2, OneDVoidState | None, fatigue_state | None)
```

The reduction changes the mechanics provider. It does not reduce lifecycle
physics or state ownership to a scalar damage variable.

## Ownership classes

| Class | V5 source | V3 owner | Rule |
| --- | --- | --- | --- |
| Physical | `VoidPhase`, site/cavity coordinate, radius, 2-D area, inventory, lineage | `OneDVoidState` | Phase names are identical. Area remains `pi*R_void^2` per unit thickness. |
| Stochastic | birth, stabilization, healing, ligament, and downstream actions/thresholds; RNG state | `HazardClock` fields and `rng_state` | Checkpoint directly. An intermediate birth hit requires both a renewed threshold and the post-draw RNG state. |
| Numerical representation | subgrid versus explicit resolved cavity | `MechanicsRepresentation` | Kept separate from physical quantities. Promotion cannot change radius, inventory, clocks, ledgers, or RNG. |
| Fracture process | sharp-front emission/cleavage state and `r_tip` | unchanged V2 `ReducedState` | V3 wraps the same object at the no-void boundary. `R_void` never enters the `r_tip` law. |
| Observables | local/marginal kinetic G, cavity tensor, reaction, compliance, energy, topology | topology-specific mechanics response and oracle rows | Observables are read-only outputs and never own thresholds or lifecycle transitions. |
| Length accounting | fractured ligament, ordinary crack, free cavity span, active-front travel | `LengthLedgers` | Newly fractured material, pre-existing free span, and front advance remain separate. |

## Phase crosswalk

| V5 phase | V3 phase | Mechanics representation | Active sharp tip | Owned clocks and transitions |
| --- | --- | --- | --- | --- |
| `AVAILABLE_SITE` | exact | `INACTIVE` | existing root only | Birth clock; repeated hits renew threshold and RNG explicitly. |
| `EMBRYO` | exact | `INACTIVE` | existing root only | Stabilization and healing clocks compete by first passage. |
| `HEALED_SITE` | exact | `INACTIVE` | existing root only | Terminal site state for this single-site prototype. |
| `STABLE_SUBGRID_VOID` | exact | `SUBGRID` | existing root only | Radius and inventory may evolve; no represented-geometry event is counted. |
| `RESOLVED_VOID` | exact | `SURROGATE_ACTIVE` | existing root only | Preconnection response map supplies candidate G and cavity tensor. |
| `CONNECTED_VOID` | exact | `CONNECTED_CAVITY` | none | Far cavity surface owns the downstream clock. No sharp-tip G or `r_tip` exists. |
| `DOWNSTREAM_FRONT_ACTIVE` | exact | `CHILD_ACTIVE` | downstream child | Child uses the unchanged V2 fracture state and renewal policy. |
| `MERGED_OR_CONSUMED` | exact | `CONSUMED` | topology dependent | Terminal inventory/lineage state; outside the initial trajectory. |

The initial implementation is aligned, circular, centerline, monotonic, and
limited to one active void. Lateral offset, branching, multiple sites/voids,
fatigue stepping, and experimental calibration are excluded.

## Field crosswalk

| V5 field | V3 field | Conversion |
| --- | --- | --- |
| `VoidSite.site_id` / `Cavity2D.parent_site_id` | `site_id` | Exact string. |
| `center_m[0]` | `x_void_m` | Exact laboratory coordinate. |
| `center_m[1]` | absent dynamic coordinate | Must equal zero for the aligned domain; retained in oracle rows for rejection/audit. |
| `Cavity2D.radius_m` | `radius_m` | Exact; distinct from `ReducedState.tip_radius_m`. |
| `Cavity2D.area_m2` | `area_m2` | Exact `pi*R_void^2` convention. |
| `VoidSite.hits`, `required_hits`, `candidate_weight` | `hit_count`, `required_hit_count`, `candidate_weight` | Exact. |
| site birth/stabilization/healing clocks | same named `HazardClock` fields | Exact action and threshold. |
| ligament directional hazard | `ligament` | Separate from birth/stabilization/healing. |
| connected far-surface hazard | `downstream` | Separate owner; never assigned to a dormant sharp tip. |
| available/consumed defect inventory | same named V3 fields | Area-conserving debit/return. |
| cavity/event lineage | `lineage` | Append-only event identity. |
| production RNG state | `rng_state` | JSON checkpoint payload, retained independently of V2 fracture RNG state. |
| V5 length-ledger family | `LengthLedgers` | Reduced to three independent physical accounts plus connected free-surface extent. |
| active branch endpoint | `active_tip_position_m` | `None` throughout the connected dormant interval. |
| complete accepted-state fingerprint | `source_v5_state_sha256` / oracle identity | Provenance only; never a surrogate state coordinate. |

## Geometry and length rules

For the aligned geometry:

```text
x_near = x_void - R_void
x_far  = x_void + R_void
L      = x_near - x_tip
```

Ligament rupture adds `L` to newly fractured material. Connection records the
void diameter as connected free-surface extent while leaving the active tip
undefined. Downstream child creation adds `2*R_void` to the free-span ledger
and the finite child length to newly fractured material. Child continuation
adds only its committed material advance. The void diameter is never counted
as fractured material.

## Mechanics-map contract

Three independent maps are required:

1. `PRECONNECTION`: additive candidate `Delta G_void`, the four cavity tensor
   components, compliance correction, and potential-energy correction.
2. `CONNECTED_CAVITY`: cavity tensor, compliance, and energy; sharp-tip G is
   structurally forbidden.
3. `DOWNSTREAM_CHILD`: additive candidate `Delta G_void`, compliance, and
   energy; cavity tensor is optional when the child law does not consume it.

The first implementation uses complete Cartesian grids in `L/R_void` and
`R_void/l_star`, with piecewise-linear tensor-product interpolation. Queries
outside the qualified rectangle raise `MechanicsMapDomainError`. The no-void
path bypasses the map and returns the exact V2 values. Load-factorized map
coordinates are permitted only after the oracle data establish stress-linear
and energy/G-quadratic scaling.

Mechanics-map rows contain mechanics observables rather than Arrhenius rates.
Material parameters must remain frozen while these maps are identified.

## Retained-evidence audit result

The exporter reads only retained V5 checkpoint JSON and static-case JSON. Each
normalized row records the source file digest, source implementation identity,
accepted-state identity when present, topology identity, geometry, phase,
mechanics observables, and length ledgers.

The retained evidence supplies the complete lifecycle state sequence, energy,
static reaction/compliance, and length ledgers. It does not supply a complete
aligned geometry-family table containing both candidate `G_kinetic` and the
four fixed-arc cavity tensor components. The exporter therefore classifies the
initial map dataset as `BLOCKED_ORACLE_COORDINATES_INCOMPLETE`. No values are
inferred from stress extrema, no Arrhenius rate is fitted, and no broad new
two-dimensional campaign is launched.

## Initial terminal classifications

```text
ONE_D_V3_VOIDING_STATE_AND_EVENTS = PASS_BOUNDED_ALIGNED_MONOTONIC_PROTOTYPE
ONE_D_V3_VOID_MECHANICS_MAP = CONTRACT_PASS_FIT_BLOCKED_ORACLE_COORDINATES_INCOMPLETE
ONE_D_V3_MONOTONIC_2D_TRANSFER = BLOCKED_UNTIL_ALIGNED_MAP_FAMILIES_EXIST
FATIGUE_IMPLEMENTATION = NOT_STARTED_BY_CONTRACT
R_TIP_EQUALS_R_VOID = FORBIDDEN
SCALAR_VOID_DAMAGE_VARIABLE = ABSENT
```
