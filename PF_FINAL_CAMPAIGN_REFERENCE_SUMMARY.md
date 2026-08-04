# PF Final Campaign (v10.2.30, four-class, theta=0, rate extremes) — Reference Summary

Principal PF physical-correspondence reference bank for this workspace,
per explicit instruction, superseding the earlier, now-unavailable
v10.4.1 selective_reuse Peak/1000K case as the primary reference (that
older case's own provenance record remains valid for what it audited at
the time — see `PF_REFERENCE_REGENERATION_CONTRACT.md`).

Campaign root (read-only; never modified by this workspace):
```
/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1/
  runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1/
```

Full per-case detail: `PF_FINAL_CAMPAIGN_REFERENCE_MANIFEST.json` (this
workspace root). Derived comparison tables:
`runs/pf_final_campaign_v10230/pf_{case_inventory,first_passage_vs_temperature,
rcurve_summary,event_statistics,material_class_ordering,rate_sensitivity_summary}.csv`
(regenerate with `scripts/inventory_pf_final_campaign_v10230.py` then
`scripts/characterize_pf_final_campaign_v10230.py`; both are read-only
against the PF repository and safe to re-run any time).

## Campaign structure

- 3 loading-rate directories: `rate0p01x`, `rate1x`, `rate100x` (plus a
  `rate_overlay_by_mechanism/` directory of pre-existing cross-rate plots
  and a `rate_overlay_source_data.csv`, not re-derived here).
- 4 material-class options per rate, exact production option IDs:
  - Peak: `v913_paper_peak01_0242980_persistent_sites`
  - DBTT: `v913_paper_dbtt01_0202500_persistent_sites`
  - weak-T: `v913_paper_weakT01_0129902_persistent_sites`
  - ceramic: `v913_paper_ceramic01_0077080_persistent_sites`
- 12 temperatures per class per rate: 300, 600, 800, 900, 950, 1000, 1050,
  1100, 1150, 1200, 1250, 1300 K. theta=0 throughout. Seeds are
  temperature-specific (e.g. Peak/1000K uses seed 8666, matching the
  earlier v10.4.1 reference case exactly), base site seed 3621.
- **3 × 4 × 12 = 144 total cases. All 144 are `complete_target_extension`
  (COMPLETE marker present, exit_code 0, `stage3_case_status.json`
  `status_complete: true`) with zero missing/malformed fields detected by
  the inventory script.** This is a fully completed, internally consistent
  campaign — do not assume otherwise from directory naming alone (per
  instruction, this was verified by inspecting files, not inferred).

## Kernel provenance caveat (read before any controlled-K_J diagnostic)

Every case's `command.sh` references the same kernel cache path:
```
.../PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/v10_2_28_kernel_cache/
  1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json
```
This campaign's case files are dated **2026-07-29**. This exact
kernel_cache path was independently observed **during this workspace's
own session on 2026-08-03/04** to change content (SHA-256 changed from
`a85b57ad9e...` to `d41b08f69a...`) due to unrelated, concurrent activity
in that same live PF repository. Re-hashing it now for this manifest
gives `d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3`
for every case — but **this is not guaranteed to be what the campaign
actually consumed at runtime**, since no case's own audit JSON records a
byte-level kernel fingerprint (only shape/policy metadata, e.g.
`v10_2_30_hazard_energy_gate_audit.json`'s
`engine.signed_burgers_shared_physics.kernel` block). Per the
reference-artifact policy: freeze whatever kernel bytes are used for a
specific controlled diagnostic into
`reference_inputs/pf_final_v10_2_30/<rate>/<option>/<temperature_seed>/`
at the time of use, record its actual hash, and do not assume it
reproduces exactly what generated the original PF steps table — the
steps table itself (the actual PF *result*) is unaffected by this
caveat, since it is a frozen output file, not something regenerated from
the kernel on read.

## Key characterization findings (rate1x, first-passage K_J vs. temperature)

| T (K) | Peak | DBTT | weak-T | ceramic |
|---|---|---|---|---|
| 300  | 23.7 | 23.4 | 19.3 | 13.6 |
| 600  | 22.5 | 23.3 | 20.7 | 13.3 |
| 800  | 41.2 | 22.7 | 20.6 | 11.9 |
| 900  | 56.7 | 22.4 | 20.6 | 12.6 |
| 950  | 56.0 | 21.7 | 19.5 | 12.3 |
| 1000 | 53.9 | 56.2 | 20.7 | 12.0 |
| 1050 | 53.5 | 55.3 | 20.8 | 11.4 |
| 1100 | 52.7 | 55.6 | 20.3 | 10.7 |
| 1150 | 51.6 | 54.7 | 20.4 | 10.4 |
| 1200 | 49.1 | 54.1 | 19.9 | 10.1 |
| 1250 | 50.4 | 54.7 | 20.0 | 9.8 |
| 1300 | 48.7 | 50.5 | 20.0 | 9.3 |

(K_J values in MPa√m, first-passage/initiation, i.e. `Kc_first_MPa_sqrt_m`
from each case's `stage3_case_status.json`.)

**The four classes behave exactly as their names suggest** — this is a
strong, clean qualitative signature to reproduce in FEM:

- **Peak**: a *broad, gradual* transition starting around 800 K and
  developing into a plasticity-toughened plateau (~50-57 MPa√m) by 900-950 K,
  slowly decaying at higher T. This is the paper's primary reference class.
- **DBTT**: a *sharp, narrow* transition — essentially flat at ~21-23
  MPa√m through 950 K, then jumping to ~50-56 MPa√m at exactly 1000 K and
  staying there through 1300 K. Far sharper than Peak's transition, as its
  name implies.
- **weak-T**: essentially *flat* across the entire 300-1300 K range
  (19.3-20.8 MPa√m) — no transition at all, confirming "weak temperature
  dependence" by construction.
- **ceramic**: *monotonically decreasing* with temperature (13.6 → 9.3
  MPa√m), no transition, lowest toughness of the four classes throughout
  — consistent with a brittle, ceramic-like response with no compensating
  thermally-activated plasticity.

**Ordering** (first-passage K_J, low to high) is `ceramic < weak-T <
{Peak, DBTT}` at essentially every temperature and rate, with Peak and
DBTT swapping relative order across their respective transitions (both
converge to a similar ~50-56 MPa√m plateau above ~1000 K, but arrive there
via different transition shapes). See
`pf_material_class_ordering.csv` for the full per-(rate, temperature)
ranking table.

## Rate sensitivity (Peak class, first-passage K_J, MPa√m)

| T (K) | rate0.01x | rate1x | rate100x | max/min ratio |
|---|---|---|---|---|
| 300  | 22.7 | 23.7 | 24.7 | 1.08 |
| 600  | 24.7 | 22.5 | 23.4 | 1.10 |
| 800  | 55.0 | 41.2 | 23.4 | 2.35 |
| 900  | 56.1 | 56.7 | 26.8 | 2.11 |
| 1000 | 66.9 | 53.9 | 39.0 | 1.72 |
| 1150 | 92.1 | 51.6 | 54.0 | 1.79 |
| 1300 | 89.0 | 48.7 | 51.1 | 1.83 |

**Rate sensitivity is weak at low temperature (300-600 K, ratio ~1.1,
cleavage-dominated, little time-dependent plasticity) and strong once
thermally-activated plasticity/emission become active (800 K and above,
ratio 1.7-2.35)** — slower loading (rate0.01x) consistently gives *higher*
first-passage toughness in that regime (more time for stress relaxation
via plasticity/emission before cleavage fires), a physically sensible
signature. See `pf_rate_sensitivity_summary.csv` for the full table
(all four classes).

## Recommended diagnostic anchor cases (Stage A)

Based on the characterization above:

| Class | Anchor temperatures | Rationale |
|---|---|---|
| Peak | 300 K, 1000 K | Below transition (cleavage-dominated, ~24 MPa√m) and on the plasticity plateau (~54 MPa√m). 1000 K/seed 8666 is also the pre-existing reference case this workspace's controller infrastructure was already built and tested against (`ad2c446` onward). |
| DBTT | 900 K, 1000 K, 1100 K | Below (~22), essentially at (transition occurs between 950 and 1000 K exactly), and above (~56) the sharp transition. |
| weak-T | 300 K, 1300 K | Confirms flatness holds at both extremes (~19-20 MPa√m expected at both). |
| ceramic | 300 K, 1300 K | Confirms the monotonic-decrease trend holds at both extremes (~13.6 → ~9.3 MPa√m expected). |

This gives 9 anchor cases (vs. the user's suggested ~7; DBTT gets a third
point to bracket its sharp transition precisely, per the governing
instructions' explicit "below, near, and above transition" requirement).
Use `rate1x` for all anchors initially, per instruction, unless campaign
analysis (above) establishes another rate as more appropriate — Peak's
rate sensitivity table shows rate1x sits in the middle of the
rate0.01x/rate100x spread at every anchor temperature, making it a
reasonable default baseline.

## What this characterization does NOT yet cover

- Full R-curve *shape* (K_J vs. crack extension) beyond the summary
  first-passage/final-value pair in `pf_rcurve_summary.csv` — the full
  per-case K_J(a) trajectory is available in each case's `steps_*K.csv`
  (already hashed and path-recorded in the manifest) but not yet
  extracted into a shared comparison table. Extract this for a specific
  anchor case if/when a native-loading FEM R-curve needs a point-by-point
  PF comparison, not speculatively for all 144 cases.
- B(t)/N_em(t)/MPZ *time series* (only initial/final scalars are in the
  manifest so far, from `stage3_case_status.json` and the steps table's
  last row) — same policy: extract per-case time series only when a
  specific controlled-K_J diagnostic needs them.
- Event-length *distributions* (only count/mean/min/max per case so far)
  — sufficient for Statistical-correspondence screening; extract full
  distributions if a specific comparison needs them.

Per instruction: stop extending generic characterization infrastructure
here. The next new extraction should be driven by a specific selected
comparison (an anchor case's controlled-K_J diagnostic or native-loading
FEM run), not built speculatively ahead of need.
