# v10.0.5.6 — KJ audit, stochastic diagnostics, and first-passage bracket

## Scope

v10.0.5.6 is a diagnostic and campaign-control point release on top of
v10.0.5.5. It does not alter the Arrhenius cleavage or emission surfaces,
finite source capacity, transport rates, shielding, blunting, cohesive law,
or bulk constitutive response.

It addresses three issues:

1. relabel the stochastic outputs using authoritative scheduler and MPZ fields;
2. audit the conversion from remote reaction force to nominal stress and from
   the FEM state to the domain-integral KJ;
3. find a remote stress range that brackets stochastic first passage at a
   user-selected cycle horizon.

No repeated run of the same stochastic seed is required.

## Remote-stress audit

The plane-strain FEM returns the top reaction force per unit out-of-plane
thickness. Gross nominal stress is therefore

```text
sigma_gross = abs(F_top) / specimen_width
```

For the default specimen:

```text
width  = 2.0 mm
height = 4.0 mm
a0     = 0.5 mm
a0/W   = 0.25
```

The existing stress calibration uses this gross-section definition correctly.
The principal defect was the inherited primary/cluster J contour: with
`L_pz=100 um`, the default cluster parameter produced an actual outer radius of
8 mm. The crack tip is only 0.5 mm from its nearest specimen boundary, so that
circular domain cannot close inside the body.

v10.0.5.6 therefore requires an explicit, audited actual cluster-J outer
radius. The contour sweep:

- rejects any contour that intersects a specimen boundary;
- uses a conservative 80% boundary-clearance limit;
- requires at least 12 active J-domain elements;
- selects a plateau in KJ/sigma across at least three radii;
- compares the plateau with a finite-width single-edge-crack reference;
- checks reaction-stress accuracy and KJ linearity using a second load amplitude.

The first-passage bracket is not authorized unless all checks pass.

## Corrected stochastic diagnostics

The no-rerun relabeler and new campaign runner write:

```text
stochastic_scheduler_mode
stochastic_event_rate_per_cycle
stochastic_expected_state_events
base_cycle_limiter
final_cycle_limiter
cycle_limiter_label
source_budget_total
source_consumed_final
source_remaining_final
mobile_count_final
retained_count_final
active_count_final
cumulative_emitted
stochastic_emission_channel_events
```

`stochastic_emission_channel_events` is the number of nonzero system/channel
realizations. It is not the emitted source count. The authoritative emitted
source count is `cumulative_emitted` or `source_consumed_final`.

## First-passage bracket

After a contour plateau is selected, the bracket runner uses its measured
KJ/sigma slope to convert requested KJmax values to remote stress ranges. It:

1. evaluates ascending KJ targets;
2. expands the range if all cases fail or all survive;
3. identifies the highest physical no-first-passage case and the lowest
   first-passage case;
4. performs configurable bisection refinements.

A case stopped by `MAX_BLOCKS` is rejected as numerical censoring and cannot
serve as either side of the physical bracket.

## New files

```text
arrhenius_fracture/kj_audit_v10056.py
run_v10_0_5_6_stochastic_delta_sigma.py
run_v10_0_5_6_stochastic_delta_sigma_audited.py
run_v10_0_5_6_kj_audit_bracket.py
run_v10_0_5_6_kj_audit_bracket_audited.py
run_v10_0_5_6_kj_audit_bracket.sh
relabel_v10_0_5_5_stochastic_outputs_v10056.py
tests/test_v10056_kj_audit_bracket.py
tests/test_v10056_audited_wrappers.py
```

## Local validation

```bash
python -m py_compile \
  arrhenius_fracture/kj_audit_v10056.py \
  run_v10_0_5_6_stochastic_delta_sigma.py \
  run_v10_0_5_6_stochastic_delta_sigma_audited.py \
  run_v10_0_5_6_kj_audit_bracket.py \
  run_v10_0_5_6_kj_audit_bracket_audited.py \
  relabel_v10_0_5_5_stochastic_outputs_v10056.py

bash -n run_v10_0_5_6_kj_audit_bracket.sh

pytest -q \
  tests/test_v10053_audited_source_preflight.py \
  tests/test_kinetic_fatigue_v10053.py \
  tests/test_v10054_vhcf_first_passage.py \
  tests/test_v10055_stochastic_vhcf.py \
  tests/test_v10056_kj_audit_bracket.py \
  tests/test_v10056_audited_wrappers.py
```

## Relabel the completed five seeds without rerunning

```bash
python relabel_v10_0_5_5_stochastic_outputs_v10056.py \
  runs/v10_0_5_5_DBTT_700K_350MPa_stochastic_smoke_seed1_v1 \
  runs/v10_0_5_5_DBTT_700K_350MPa_stochastic_smoke_seed2_v1 \
  runs/v10_0_5_5_DBTT_700K_350MPa_stochastic_smoke_seed3_v1 \
  runs/v10_0_5_5_DBTT_700K_350MPa_stochastic_smoke_seed4_v1 \
  runs/v10_0_5_5_DBTT_700K_350MPa_stochastic_smoke_seed5_v1
```

## Recommended staged run

First run only the contour and loading audit:

```bash
ACTION=audit \
MATERIAL=DBTT \
TEMPERATURE_K=700 \
CONTOUR_OUTER_UM="80 100 140 180 240 300 360 400" \
OUTROOT=runs/v10_0_5_6_DBTT_700K_KJ_contour_audit_v1 \
bash run_v10_0_5_6_kj_audit_bracket.sh
```

The audit must report `plateau_selected`, pass the KJ/LEFM ratio check, and pass
the second-amplitude stress/KJ linearity check.

Then run the first-passage bracket using the selected contour:

```bash
ACTION=bracket \
MATERIAL=DBTT \
TEMPERATURE_K=700 \
STOCHASTIC_SEED=1 \
CYCLES_MAX=1e7 \
TARGET_KJMAX="2 4 6 8 10 12 16 20 24" \
SELECTED_CONTOUR_JSON="$PWD/runs/v10_0_5_6_DBTT_700K_KJ_contour_audit_v1/selected_KJ_contour_v10_0_5_6.json" \
OUTROOT=runs/v10_0_5_6_DBTT_700K_first_passage_bracket_seed1_v1 \
bash run_v10_0_5_6_kj_audit_bracket.sh
```

The bracket is defined at the selected `CYCLES_MAX`; it is not claimed to be a
cycle-independent material threshold.
