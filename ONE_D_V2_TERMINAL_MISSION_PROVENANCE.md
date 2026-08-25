# One-dimensional V2 terminal-mission provenance

## Mission boundary

- Terminal branch: `codex/oneD-v2-terminal-predictive-program`
- Starting commit: `e2ed84c9bad5086a953505a8e24129c07a494cdc`
- Frozen diagnostic tag: `oneD-v2-predictive-model-diagnostic-baseline`
- Frozen-baseline classification: **V2 predictive-model diagnostic baseline, partially qualified, no new material row promoted.**
- This branch was created in a clean worktree. Untracked archives in the parent worktree were not copied, modified, or adopted as inputs.

The source commits below are the commits recorded by the qualified V2 audit artifacts. Ambient repository heads and uncommitted files are not silently substituted for these pinned sources.

## Source repositories and commits

| Role | Repository | Qualified source commit | Ambient head at mission start | Use policy |
|---|---|---|---|---|
| Primary reduced-model repository | `/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v9_13_dbtt_temperature_shelf` | `e2ed84c9bad5086a953505a8e24129c07a494cdc` | terminal worktree starts at the same commit | Writable only through the terminal branch |
| Authoritative PF source and archived trajectories | `/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1` | `ab6279331d050919f78e2e7f278bc332466e34f8` | `9e884fb0b0845da621d2612bdf1042e481b8df49` | Pinned source content and clean analysis-only worktrees for bounded validation; canonical results are read-only |
| Corrected FEM/CZM source and trajectories | `/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude` | `931bed66913afc970117bce900805ecc9b6225f8` | `17f8bf67c0138b453eadff4b509afa80ebaa2c04` | Read-only; no new 2-D FEM/CZM production run |

The PF and FEM/CZM source repositories had ambient state beyond the pinned audit commits at mission start. That state is not evidence for this mission unless a later artifact explicitly records a clean source worktree, exact commit, and content hashes.

## Mechanics and local-drive inputs

| Input | Semantics | Qualified interval | SHA-256 |
|---|---|---:|---|
| `analysis_outputs/oneD_v2_mechanics_maps_and_baselines/oneD_v2_pf_mechanics_map.csv` | PF production-discrete native J/KJ mechanics; finite-width stiffness-killed wake | 0–1000 µm | `b508d5ba6a7562e53f5280ffecf015d3c6f8c801d0e03abe8bc9d7ddef8f0fc1` |
| `analysis_outputs/oneD_v2_mechanics_maps_and_baselines/oneD_v2_fem_native_mechanics_map.csv` | FEM/CZM production-native domain J/KJ mechanics; traction-free codimension-one interface | 0–1000 µm | `21b13ca048482b560e2427dfd6beb359fcc406c421cafa488567c5fc27b4ffbd` |
| `analysis_outputs/oneD_v2_mechanics_maps_and_baselines/oneD_v2_fem_qualified_G_map.csv` | Qualified FEM/CZM structural G/K_G interpretation output | 0–1000 µm | `1e8a77a642d6c086582963f9eb576227ebb658b6a392910ce81f3c94186e0607` |
| `analysis_outputs/oneD_v2_predictive_model/oneD_v2_pf_source_drive_map.csv` | PF deterministic normalized source-drive probes | extension 0–1000 µm; radius 1–100 µm | `cdd5760511ff782f55333f7492ed3b1b01ebc2fea5eb239a254896217b589eb8` |
| `analysis_outputs/oneD_v2_predictive_model/oneD_v2_fem_source_drive_map.csv` | FEM/CZM deterministic normalized source-drive probes | extension 0–1000 µm; radius 1–100 µm | `59b9a04006bf92ea2dfc69ed6720698d68bc276f56e057d1d115bfdff3a3e247` |
| `mpz_v9_13_v10222_transfer_common_physics.json` | Shared common-physics configuration | fixed audited input | `bd356ce6d9d17c5aa7be1f6d975577d6484706daa7e39a19eb52f006b7ee6beb` |

All qualified maps use bounded interpolation and fail closed outside their recorded domain. PF native J is not relabelled as continuum G. FEM/CZM native J remains the kinetic input; qualified structural G remains a separate interpretation output.

## Lifecycle evidence

| Backend | Accepted path | Rejected trials / rollback | Post-event renewal | Terminal behavior | Evidence SHA-256 |
|---|---|---|---|---|---|
| PF | Qualified production adaptive predictor/retry and exactly-once commit | Rejected-trial state preservation qualified | Qualified | Late geometry veto terminates fail-closed; rollback-and-continue is unsupported by production | `745a1719777d5769687d3c6b23103e0cd67e6373a19fafed3038ef7fef7bc9f3` |
| FEM/CZM | Qualified normal accepted-event transaction | Joint mechanics, topology, process-zone, bulk-history, clock, threshold, and RNG restoration qualified | Qualified | Transactional rollback-and-continue qualified | `ff587965dd02b9eb3e8e22ccf614cfb6e51b2d096a6cbbb3e340fc6b80e7f5f6` |

The lifecycle evidence is backend-specific. No common event-lifecycle surrogate is assumed by the terminal program.

## Historical material registry

Registry:

`/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1/arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv`

SHA-256: `4ba723c80abcfdd7101cee0afaa9b4104eebf2a8e847892528128664829c7760`

| Response class | Historical candidate ID | Status at mission start |
|---|---|---|
| Peak-like | `v913_zeroD_sobol_0242980` | Reference row; not newly promoted |
| DBTT | `v913_zeroD_sobol_0202500` | Reference row; not newly promoted |
| weak-T / FCC-like | `v913_zeroD_sobol_0129902` | Reference row; not newly promoted |
| ceramic-like | `v913_zeroD_sobol_0077080` | Reference row; not newly promoted |

One shared material row per class is required across both mechanics providers. Backend-specific lifecycle and reduction parameters are not material parameters.

## Canonical comparison references

The authoritative PF 1000 K, rate-1×, theta-0 references are:

- Peak, seed 8666: `runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1/rate1x/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666`
- DBTT, seed 1008666: `runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1/rate1x/v913_paper_dbtt01_0202500_persistent_sites/T1000K_th0_seed1008666`

The committed derived transaction table is `analysis_outputs/oneD_v2_mechanics_maps_and_baselines/pf_2d_event_transactions_v2.csv`; it records one target-censored Peak physical avalanche and two DBTT physical avalanches. Existing corrected Peak and DBTT FEM/CZM trajectories, as represented by the committed V2 baseline/audit artifacts and pinned FEM/CZM source commit above, are the only FEM/CZM 2-D references permitted in this mission.

## Explicit exclusions

- `/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_codex_v10_2_30` is outside mission scope. Its active production state will not be stopped, resumed, checked out, reset, modified, or reused.
- Canonical PF result directories are not output targets.
- No new 2-D FEM/CZM simulation is authorized.
- The `e2ed84c` reports and machine outputs remain the diagnostic baseline; terminal results will be written as a superseding, separately versioned artifact set.
