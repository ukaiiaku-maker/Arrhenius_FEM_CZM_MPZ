# PF Reference Artifact Regeneration Contract

This document is for the **Codex PF workspace**. It specifies exactly what
this FEM/CZM parity workspace needs regenerated, with every parameter
captured from this workspace's own audit trail before the live artifacts
became unavailable. It does not ask the PF workspace to change its own
conventions — only to reproduce (or knowingly re-issue) two specific
artifacts and hand back verifiable checksums.

## Executable request (concise — read this section first)

```text
PF repository:            https://github.com/ukaiiaku-maker/PF-fracture-fatigue.git
Branch (best candidate):  v10.4.1-bulk-detailed-balance
Exact commit (best candidate):
                           59da598c7b19b70068cf8f308b09d1362db05f6f
                           ("Keep selective reuse flag at campaign level
                           only", 2026-07-31 21:36:56 -0700)
Already checked out at:   /Volumes/Data/Data/Nanopillar_calculation/
                           PF-fracture-fatigue_v10_4_1_bulk_detailed_balance
                           (a git worktree of the same repository, already
                           at this exact commit — confirmed identical to
                           origin/v10.4.1-bulk-detailed-balance's tip, no
                           further commits on that branch)

Confidence caveat:        NOT independently verified as the literal
                           generating commit -- PF run outputs do not
                           record git provenance (checked several sibling
                           campaign manifest/lock JSON files in this
                           repository; none contain a commit hash). This
                           is the best available candidate, identified by
                           searching commit messages for "v10.4.1",
                           "selective_reuse"/"selective reuse", "base3621",
                           and "rate1x" and cross-checking commit dates
                           against the missing run family's naming
                           convention. If output columns or first-passage
                           values do not match the reference values below
                           once regenerated, try nearby commits on the
                           same branch (see "Development note" below).

Option ID:                v913_paper_peak01_0242980_persistent_sites
                           (candidate_id: v913_zeroD_sobol_0242980)
Temperature:               1000 K
Theta:                     0 deg
Hazard seed:                8666
Base site seed:             3621
Loading/timestep:           dU=2.0e-7 m, dt=8.4 s, n_stagger=2,
                             tip_h_fine=1.0e-6 m, tip_ratio=1.20,
                             da_phys=5.0e-6 m, adaptive_event_target=0.15,
                             maximum_fronts=1, loading rate multiplier
                             rate1x, target_projected_extension_um=1000,
                             bulk_plasticity_mode=full_field
Kernel construction:        see "Artifact 2" configuration block below
                             (configuration_fingerprint
                             1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5)
Expected output columns:    step, dt_cur_s, J_effective_direct_J_per_m2,
                             J_signed_direct_J_per_m2, KJ_Pa_sqrtm, B,
                             N_em, crack_extension_m (plus any other
                             columns the PF steps table normally emits --
                             these 8 are the ones this workspace's reader
                             requires and validates)
First-passage reference:    row index 176 (0-based; PF step=177),
                             physical time 1486.8 s (=177*8.4),
                             KJ=53.9019 MPa√m, J≈6531 J/m² -- this is
                             the first row with crack_extension_m>0 in
                             the original (now unavailable) artifact;
                             use it to sanity-check a regenerated file
                             before handing it back
Destination for artifacts:  hand back the two files (or their new paths)
                             for this workspace to copy into
                             Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude/
                             reference_inputs/pf_v10_4_1_theta0_peak1000K_seed8666/
                             (this workspace performs the copy + checksum
                             step itself -- do not write into this FEM/CZM
                             workspace directly)
Checksums to verify:        steps_1000K.csv sha256 =
                             666d839ae8ef4267ae4ac7ab52a4220de676525bfb70917744e25a119c291b0c
                             family.json sha256 =
                             a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a
                             (bit-exact match is the goal but not
                             guaranteed -- see the kernel section's
                             caveat about the same configuration_fingerprint
                             already having produced two different
                             family_sha256 values in this repository)
```

### Development note: nearby commits if `59da598` does not reproduce

```text
706ea8f Bump package metadata to v10.4.1 detailed balance   (2026-07-31 14:38:54, earlier same day)
2bc95e0 Restore frozen v10.4.1 launcher builder
4f2a846 Validate native and previously reused v10.4.1 cases correctly
22710b6 Add v10.4.1 completed-case materializer for v10.4.2
bbfbcb1 Add audited reuse of completed v10.4.1 fracture cases
9a5e8d3 Include selective reuse tests in v10.4.1 validation
6ab901c Permit only audited selective reuse in v10.4.1 scheduler
5f5fb09 Add selective reuse materialization CLI
345fe70 Add selective reuse audit CLI
018bb4a Add audited v10.4.0 to v10.4.1 case reuse machinery
```
All found via `git log --all --oneline --grep="v10.4.1\|bulk_PT\|selective_reuse" -i`
in this repository; none independently confirmed as more or less likely
than `59da598` without actually running them.

## Why this is needed

This FEM/CZM workspace
(`Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude`, branch
`claude/v10.0.5.18.4.0-j-controlled-loading`) is pinned by CLAUDE.md and
FEM_CZM_HANDOFF.md to two exact artifacts in the PF reference repository.
Both became unavailable at their documented paths during this session's
work (2026-08-04), due to what appears to be unrelated, concurrent
regeneration activity in that repository's live `runs/` directory (a
"v10.4.x bulk plasticity campaign" with file timestamps in the same
window). See `CLAUDE_PROGRESS.md`'s "EXTERNAL BLOCKER" section in this
workspace for the full incident record.

A systematic, time-boxed search of this whole `Nanopillar_calculation`
tree (40 candidate `steps_1000K.csv` files, 2 candidate `family.json`
files, hashed exhaustively) found no exact-hash match for either
artifact. No local filesystem/APFS snapshots exist for the volume
(`tmutil`/`diskutil apfs listSnapshots` both empty for
`/Volumes/Data`). The search is closed; do not re-run it — regenerate
instead.

## Artifact 1: `steps_1000K.csv` (PF prefracture/growth trajectory)

**Required output**: the PF step-table CSV for the Peak parameterization,
1000 K, theta=0°, hazard seed 8666, run family
`v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1`.

**Expected SHA-256** (the exact artifact this workspace audited before it
disappeared): `666d839ae8ef4267ae4ac7ab52a4220de676525bfb70917744e25a119c291b0c`

**Exact run parameters** (captured from this workspace's own
`campaign_configuration.txt`/`local_preflight.json` records and from the
FEM production entry's forwarded PF-matched controls, before the external
change):

```text
parameter_option:              v913_paper_peak01_0242980_persistent_sites
candidate_id:                  v913_zeroD_sobol_0242980
temperature_K:                 1000
theta_deg:                     0
hazard_seed:                   8666
base_site_seed:                3621   (from the run-family name "...base3621_v1")
target_projected_extension_um: 1000
dU_m:                          2.0e-7
dt_s:                           8.4
n_stagger:                      2
tip_h_fine_m:                   1.0e-6
tip_ratio:                      1.20
da_phys_m:                      5.0e-6
adaptive_event_target:          0.15
maximum_fronts:                 1
bulk_plasticity_mode:           full_field
bulk_model:                     v10.4.1_bulk_peierls_taylor_detailed_balance_exact_port
crack_backend (PF side):        sharp_wake
loading rate multiplier:        rate1x  (the run family has rate0p01x and
                                 rate100x siblings elsewhere in this repo;
                                 the required artifact is specifically the
                                 rate1x case)
```

**Selected-row barrier/kinetic parameters** (the exact "shared_contract"
values this workspace's own production run captured from the persistent
site registry for this candidate — reproduce with the same row, do not
re-derive):

```text
candidate_fingerprint_sha256:  ce47ec844b3d7a9812aa0b077351fc6d022c292b39b84bd4056871395f7f4f52
barrier_fingerprint_sha256:    8d8639702ffdafa621649d714ad55948ef40733cc8b1f46a1280aafe1d4d2876
parameter_source:               PF_v10.2.22_exact_persistent_site_registry

Tref_K:                         481.33
cleave_G00_eV:                  4.011803912930191
cleave_gT_eV_per_K:              0.0081468212408944
cleave_sigc0_GPa:                7.349638969637454
cleave_sT_GPa_per_K:            -0.0012099820682778
cleave_exp_a:                    1.3227705595549195
cleave_exp_n:                    1.2558047282509506
cleave_floor_frac:               0.0184513317135146
emit_G00_eV:                     2.0446315124630927
emit_gT_eV_per_K:                 0.000941995373927
emit_sigc0_GPa:                  7.205215461552143
emit_sT_GPa_per_K:                0.0013325903080403
emit_exp_a:                       0.0858719444833695
emit_exp_n:                       1.1423492641188204
emit_floor_frac:                  0.0058420097608974
peierls_H0_eV:                    6.368386985268444
peierls_activation_entropy_kB:   12.2265643812716
peierls_exp_a:                    1.9097672308795155
peierls_exp_n:                    1.1372957159765065
peierls_nu0_s:                    1.0e12
taylor_H0_eV:                     1.6986502355337143
taylor_activation_entropy_kB:    38.45378890633583
taylor_exp_a:                     0.4561179580539465
taylor_exp_n:                     1.2506696581840515
taylor_nu0_s:                     1.0e11
source_sites_per_system:          141.0590567476921
rho_source0_m2:                   1.4102520844479966e+16
encounter_efficiency:             9.160246308716648
c_blunt:                          2.687256144359708
taylor_corr_rho_c_m2:             58780073970578.01
taylor_corr_scale:                0.8638954413677038
L_pz_um_recommended:              50.0
n_bins_recommended:               80
rho_forest_floor_m2:               5.0e12
peierls_stress_fraction:          0.5773502691896258
taylor_stress_fraction:           0.5773502691896258
```

**Required output location** (produced by the PF workspace, path is its
own choice, but report it back exactly):
`.../runs/v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/steps_1000K.csv`

## Artifact 2: kernel `family.json`

**Required output**: the signed shielding-kernel family used by the above
PF case.

**Expected SHA-256** (the exact artifact this workspace audited):
`a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a`

**Exact generating configuration** (recovered from this workspace's
read-only inspection of the PF repo's own
`runs/v10_2_28_kernel_cache/registry.json`, entry keyed by
`configuration_fingerprint = 1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5`):

```json
{
  "active_station_policy_id": "v10.2.14_exact_first_last_active_stations",
  "atlas_anchor_spacing_m": 0.0002,
  "boundary_condition_id": "v10.2.27_mode_I_displacement",
  "branching_mode": "single_front",
  "crystal_C11_Pa": 523000000000.0,
  "crystal_C12_Pa": 203000000000.0,
  "crystal_C44_Pa": 160000000000.0,
  "da_phys_m": 5e-06,
  "elasticity_policy": "cubic_plane_strain_explicit_constants",
  "extra": {
    "burgers_m": 2.74e-10,
    "geometry_anchor_source": "direct_prescription_not_stochastic_capture",
    "hazard_seed_required": false,
    "kinetic_packet_length_m": 2.5e-10,
    "material_parameter_option_required": false,
    "prescribed_crack_path_policy": "forward_100_cleavage_trace",
    "prior_kernel_family_required": false
  },
  "front_direction_convention": "v10.2.27_front_direction_fix",
  "initial_crack_length_m": 0.0005,
  "interaction_length_m": 2e-06,
  "kernel_provider_id": "v10.2.28_direct_prescribed_geometry_fem_v1",
  "maximum_fronts": 1,
  "measurement_mesh_policy_id": "v10.2.28_direct_prescribed_geometry_endpoint_mesh",
  "measurement_tip_h_fine_m": 7.8125e-09,
  "measurement_tip_ratio": 1.2,
  "mechanics_backend": "v10.2.28_sharp_front_fem_direct_kernel",
  "mesh_nx": 36,
  "mesh_ny": 72,
  "mesh_policy_id": "v10.2.28_production_evolution_plus_direct_measurement_mesh",
  "minimum_elements_per_process_zone": 3.0,
  "nominal_crack_angle_deg": 0.0,
  "normalization_policy": "v10.2.28_unchanged_code_defined_activation_line_conversion",
  "notch_half_thickness_m": 8e-05,
  "process_zone_bins": 80,
  "process_zone_length_m": 5e-05,
  "process_zone_policy_id": "dynamic_tip_radius_physical_front_width",
  "profile_id": "v10_2_28_direct_prescribed_geometry_single_front",
  "schema": "v10.2.27_mechanical_kernel_configuration_v6",
  "signed_channel_convention": "v10.2.27_signed_mobile_retained_channels",
  "specimen_geometry_id": "v10.2.27_rectangular_single_edge_notch",
  "specimen_length_x_m": 0.002,
  "specimen_length_y_m": 0.004,
  "temperature_K": null,
  "temperature_dependent_mechanics": false,
  "theta_deg": 0.0,
  "tip_h_fine_m": 1e-06,
  "tip_ratio": 1.2
}
```

**Resulting kernel artifact schema** (from the FEM side's audited
compatibility record of this exact kernel, for cross-checking the
regenerated file's shape/metadata even if the byte-level content differs):

```text
schema:                          v10.2.14_active_only_real_signed_2d_shielding_atlas
point_release:                    10.0.5.14.1
artifact_kind:                    crack_extension_kernel_family
state_ids:                        E0000000, E0000200, E0000400, E0000600,
                                   E0000800, E0001000, E0001030
crack_extension_levels_m:         [0.0, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.00103]
source_active_grid_points:        80
interpolation:                    inverse_distance, neighbors=4, power=2.0
spatial_projection:                piecewise_linear_with_endpoint_hold
wake_kernel_forced_zero:           true
activation_to_line_content_by_system: [0.9124087591240877, 0.9124087591240877]
constitutive_K_shield_cap:         false
counts_are_signed_burgers_lines:   true
normalization_is_mechanically_derived: true
```

**Important caveat — regeneration may not reproduce the exact original
hash.** While investigating, this workspace observed that the *same*
`configuration_fingerprint` (`1447653d...`) in the PF repo's own live
`registry.json` currently maps to a **different** `family_sha256`
(`d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3`) than
the one this workspace originally audited
(`a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a`) — for
the identical configuration key. The registry entry also carries
`"self_consistency_converged_iteration": null`, suggesting the kernel
comes from an iterative self-consistency solve that may not be bit-exact
across separate generation runs, or that the underlying kernel-generation
code changed between the original generation and now. **If bit-exact
reproduction of the original SHA-256 is not achievable, say so explicitly
rather than silently supplying a different-hash artifact** — this FEM/CZM
workspace's CLAUDE.md pins the specific hash above as the audited
contract, and changing it anywhere in this workspace requires an explicit,
human-approved decision, not a silent substitution.

## What to hand back

1. Whether an exact-hash match for either artifact was found elsewhere (in
   which case, just the path — no regeneration needed).
2. If regenerated: the two files, plus their SHA-256, plus whichever of
   the above configuration/parameter blocks were actually used (so this
   workspace can record any deviation).
3. If regeneration cannot reproduce the exact original SHA-256: an
   explicit statement of that fact, plus the SHA-256 of what *was*
   produced, so a human can decide whether to re-audit this workspace's
   contract against the new artifact or to keep searching/wait.
4. Do not modify the live PF `runs/` tree's ongoing "v10.4.x bulk
   plasticity campaign" activity to accommodate this request — treat this
   as a side, read-adjacent regeneration task, not a reason to interrupt
   or reorganize whatever else is running there. That activity lives in
   the `PF-fracture-fatigue_v10_2_21_persistent_sites_top1` worktree
   (currently on branch `v10.2.30-hazard-energy-gated-fatigue-events`).
   Regenerate from the **separate**
   `PF-fracture-fatigue_v10_4_1_bulk_detailed_balance` worktree instead —
   it is already checked out at the candidate commit `59da598` and does
   not need to touch the other worktree at all.

## What happens on this side once artifacts are available

Per this workspace's own `CLAUDE_PROGRESS.md` "Exact next command": the
two files get copied into a dedicated immutable local reference directory
inside the FEM/CZM workspace (outside the live PF `runs/` tree),
independently re-verified against the SHA-256 values above (or the
newly-issued ones, if a human approves treating regeneration as
authoritative), and only then does Gate 1 (controlled prefracture
mechanics) of `FEM_PF_PARITY_SCORECARD.md` proceed.
