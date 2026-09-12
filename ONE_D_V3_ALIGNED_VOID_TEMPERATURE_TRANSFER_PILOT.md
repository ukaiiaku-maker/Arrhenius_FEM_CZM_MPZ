# OneD V3 aligned-void temperature-transfer pilot

## Bounded preflight result

The four-family material amendment is installed as a fail-closed pilot
contract. The paired mechanics and trajectory campaign did not start because
the required Phase M2 material mapping gate fails before any new 2-D solve:

```text
MATERIAL_MAPPING_GATE = BLOCKED_UNMAPPED_ACTIVE_FIELD
PAIRED_MATERIAL_TEMPERATURE_CELLS_RUN = 0/12
NEW_2D_MECHANICS_SOLVES = 0
NEW_2D_TRAJECTORIES = 0
NEW_1D_TRAJECTORIES = 0
```

The pinned, read-only V5 one-void driver creates its `FrontEngine`, barriers,
and initial material state internally. It does not accept a versioned fracture
row. Barrier and manifest classes elsewhere in the 2-D repository provide
possible structural destinations, but the one-void runtime does not bind those
destinations to a selected row. Running it now would silently use defaults and
would not be a material-row transfer.

Peak and DBTT identities occur in historical 2-D registries, but those
registries are not consumed by this V5 runtime. The exact focused weak-T and
ceramic-like rows occur only in the OneD V2 registry at this checkpoint. No
registry was copied into or merged with the 2-D repository.

## M0 — versioned material bundles

Each paired case owns five independent row identities:

```text
MaterialBundle = (
    fracture_material_row_id,
    void_kinetics_row_id,
    elastic_row_id,
    site_population_row_id,
    specimen_loading_row_id,
)
```

The exact fracture rows are unchanged:

| Response family | Exact OneD V2 fracture row |
| --- | --- |
| Peak | `v913_zeroD_sobol_0242980` |
| DBTT | `v913_zeroD_sobol_0202500` |
| weak-T | `oneD_v2_focused_weak_T_0016` |
| ceramic-like | `oneD_v2_focused_ceramic_like_0018` |

All four bundles use the same
`v5.voiding-config.reference-one-void/1` void-kinetics row. Fracture and void
parameter ownership remain separate. No family-specific void barrier was fit.

## M1 — common reference void kinetics

The exact V5 `arrhenius_fracture.voiding_v5.arrhenius_rates` function was
evaluated at a declared frozen hydrostatic in-plane tensor of 3 GPa on the
25 K grid from 300 through 1200 K. This was a rate-only calculation and did
not assemble or solve a 2-D mechanics state.

Because the same void row is used for every family, this scout is intentionally
material-independent. It cannot identify four material-specific feature
temperatures without the fracture rows being exactly bound to both trajectory
runtimes. `T_low = 300 K` and `T_high = 1200 K` are recorded; all four feature
anchors remain `NOT_FROZEN_M2_MATERIAL_BINDING_GATE_BLOCKED`. No paired result
was inspected when making this decision.

## M2 — exact field mapping

The machine-readable audit contains 204 records: every one of the 51 source
columns for each of the four exact rows. Each record contains the source field,
units, original full-precision CSV text, OneD target, potential 2-D structural
target, conversion, activity, classification, and reason.

Thirty-eight active fields per material lack an exact binding in the pinned V5
one-void runtime. These include the cleavage and emission surfaces, transport
parameters, initial source density and dimensions, process-zone length,
correlation controls, and `c_blunt`. A class member such as
`FrontEngine.f.c_blunt` is reported as a structural target but remains
unsupported until the selected material row is actually passed to every V5
engine construction and checkpoint/renewal path.

Seven zero-valued disabled pathways and six identity/provenance columns are
classified inactive. Unsupported active fields are never treated as defaults,
and the audit raises `BLOCKED_UNMAPPED_ACTIVE_FIELD` when execution is
requested.

## M3–M7 — gated execution ledger

The 12 planned material/temperature cells are present in the bounded ledger.
Every no-void 1-D baseline, no-void 2-D companion, void 1-D run, void 2-D run,
and void-induced increment is marked `NOT_RUN_M2_GATE_BLOCKED`. Therefore no
raw observable, increment, mechanics-map fit, fixed-state comparison, event
comparison, or trajectory comparison is claimed.

Peak, DBTT, and weak-T remain the development sentinels. Ceramic-like remains
a complete held-out material family and was not used to tune a map, tolerance,
interpolation order, or temperature.

The unexecuted downstream contract remains frozen: 18 prospective mechanics
states comprise three topology regimes, three geometry levels, and two load
levels. Each applied 2-D load requires a matched no-void companion, with
`K0_1D = sqrt(Eprime * max(G0_2D, 0))`; raw boundary opening is never used as
OneD K. Map criteria remain 5% for held-out G and cavity tensors and 2% for
compliance and energy, with complete-family holdouts and fail-closed bounds.
Raw observables and `Delta O_void = O_void - O_no_void` are separate outputs.

Every future paired run must preserve exact equality of thresholds, RNG state,
site state, initial void state, fracture state, `r_tip`/emission/shielding
state, load history, temperature history, and stop condition. Complete transfer
can pass only when all four material families pass both fixed-state and
trajectory comparisons.

The required family decisions are:

```text
MATERIAL_ROW_TRANSFER_PEAK = BLOCKED_UNMAPPED_ACTIVE_FIELD
MATERIAL_ROW_TRANSFER_DBTT = BLOCKED_UNMAPPED_ACTIVE_FIELD
MATERIAL_ROW_TRANSFER_WEAK_T = BLOCKED_UNMAPPED_ACTIVE_FIELD
MATERIAL_ROW_TRANSFER_CERAMIC_LIKE = BLOCKED_UNMAPPED_ACTIVE_FIELD

ONE_D_V3_PEAK_TEMPERATURE_TRANSFER = NOT_RUN_M2_GATE_BLOCKED
ONE_D_V3_DBTT_TEMPERATURE_TRANSFER = NOT_RUN_M2_GATE_BLOCKED
ONE_D_V3_WEAK_T_TEMPERATURE_TRANSFER = NOT_RUN_M2_GATE_BLOCKED
ONE_D_V3_CERAMIC_LIKE_TEMPERATURE_TRANSFER = NOT_RUN_M2_GATE_BLOCKED

ONE_D_V3_MECHANICS_MAP_FIT = NOT_RUN_M2_GATE_BLOCKED
ONE_D_V3_FIXED_STATE_2D_TRANSFER = NOT_RUN_M2_GATE_BLOCKED
ONE_D_V3_LOW_T_TRAJECTORY_TRANSFER = NOT_RUN_M2_GATE_BLOCKED
ONE_D_V3_TRANSITION_T_TRAJECTORY_TRANSFER = NOT_RUN_M2_GATE_BLOCKED
ONE_D_V3_HIGH_T_TRAJECTORY_TRANSFER = NOT_RUN_M2_GATE_BLOCKED
ONE_D_V3_COMPLETE_ALIGNED_MONOTONIC_TRANSFER = BLOCKED
FATIGUE_IMPLEMENTATION = NOT_STARTED_BY_CONTRACT
```

## Preserved scientific boundaries

No value was inferred for `G_kinetic`, fixed-arc cavity `sigma_nn`, `sigma_tt`,
`tau_nt`, or cavity mean stress. The existing `r_tip` law is unchanged and
`r_tip != R_void`. Fractured material, free span, and front coordinate remain
independent ledgers. Temperature remains a kinetics input. No fatigue,
lateral-offset, branching, multiple-void, calibration, experimental fitting,
or broad parameter campaign was started.

The next bounded prerequisite is a separately qualified material-bundle
injection contract in the V5 runtime. It must bind the selected row through
initial engine construction, process-owner checkpoint restore, renewal, and
child continuation before this preflight can permit any oracle or trajectory
solve.

## Machine-readable records

- `analysis_outputs/oneD_v3_aligned_temperature_transfer/material_bundles.json`
- `analysis_outputs/oneD_v3_aligned_temperature_transfer/material_field_mapping_audit.json`
- `analysis_outputs/oneD_v3_aligned_temperature_transfer/common_void_rate_scout.json`
- `analysis_outputs/oneD_v3_aligned_temperature_transfer/paired_case_ledger.csv`
- `analysis_outputs/oneD_v3_aligned_temperature_transfer/decision.json`

## V4 source-readiness update

The source-conforming V4 checkpoint is pinned to 2-D head
`9d2cb7cb63fc74f02e5f0b54b5a2d1a394fe7d05`. The retained V3 operator and its
manufactured/Kirsch PASS remain unchanged. V3 central tensor convergence is
`NOT_RUN` because its owned polygon vertex failed nominal-circle registration.

V4 centers a polygon facet on the source ray and adds exact degree-two near and
far source nodes at the nominal radius. The central DBTT evaluation therefore
passes geometry registration. Across 32, 64, and 128 angular segments it also
passes the final tensor-change, traction, fixed-window, and patch-conditioning
predicates. It fails the unchanged normal-resolution, tangential-resolution,
and global mesh-quality predicates. The terminal source decision is:

```text
DBTT_SOURCE_READINESS = BLOCKED_WITH_EXACT_V4_FAILURE_CLASS
V4_FAILURE_CLASS = NORMAL_DIRECTION_RESOLUTION + TANGENTIAL_DIRECTION_RESOLUTION + MESH_QUALITY
ORACLE_STATES_ACCEPTED = 0_OF_18
PAIRED_TRAJECTORIES_RUN = 0_OF_12
ONE_D_V3_MECHANICS_MAP_FIT = BLOCKED_CENTRAL_DBTT_V4_RESOLUTION_AND_QUALITY
ONE_D_V3_COMPLETE_ALIGNED_MONOTONIC_TRANSFER = BLOCKED_ZERO_OF_18_ORACLE_STATES
FATIGUE_IMPLEMENTATION = NOT_STARTED_BY_CONTRACT
```

No mechanics field was inferred, no map was fit, and no paired trajectory was
run. A finite activation-zone observable was not derived in this mission.
