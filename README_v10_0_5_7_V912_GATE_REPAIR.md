# v10.0.5.7 — v9.12 campaign bookkeeping and publication-gate repair

This point release is additive on top of `v10.0.5.6-kj-audit-first-passage-bracket`.
It does not change the Arrhenius barriers, moving-process-zone state, FEM,
cohesive insertion, shielding, blunting, transport, or crack-renewal physics.

## Confirmed bugs repaired

### 1. Root-level temperature-summary path

`run_mpz_v9_11_mode_i_rcurve_3T.py` writes:

```text
<run_root>/rcurve_temperature_summary.csv
```

The v9.12 campaign previously attempted to read:

```text
<run_root>/<class>/rcurve_temperature_summary.csv
```

Consequently `status`, `final_extension_um`, `target_completed`, `control_state`,
and `K_init_MPa_sqrt_m` were silently absent from v9.12 case and campaign tables.

The v10.0.5.7 runner reads the authoritative root-level file immediately after
each subprocess, selects the requested class/temperature row, copies the file
into the case directory, and rewrites `v9_12_case_summary.csv/json`.

### 2. Vacuous publication gate

The v9.12 audit passed whenever no geometry-dominated pair and no missing field
image were found. It did not require the subprocess to succeed, the target crack
extension to complete, or first passage to be observed. A one-material campaign
also passed the pairwise test vacuously.

The v10.0.5.7 gate requires, for every case:

- `subprocess_returncode == 0`;
- `status` is `complete` or `skipped_complete`;
- `target_completed == true`;
- `control_state == first_passage` with finite initiation toughness;
- the required full-field image exists.

At least one pairwise material comparison is required. The audit writes explicit
case-level failure reasons and lists failed, incomplete/censored, non-first-
passage, missing-summary, and missing-field cases separately.

### 3. Full-field entry composition

The publication entry point composes the bookkeeping repair with the validated
v9.12 command substitution that selects `run_mpz_v9_12_mode_i_rcurve.py` and the
mapped MPZ field renderer. The non-full-field repair module remains importable for
unit testing, but it is not the publication campaign entry point.

## Entry point

```bash
python run_mpz_v10_0_5_7_tip_only_material_rcurve_fullfield.py \
  --classes "ceramic weakT DBTT" \
  --T-K 700 \
  --seeds 1 \
  --target-extension-um 100 \
  --outroot runs/v10_0_5_7_material_transfer_700K_v1
```

The new audit outputs are:

```text
material_rcurve_audit_v10_0_5_7.json
material_rcurve_case_audit_v10_0_5_7.csv
material_rcurve_pairwise_audit_v10_0_5_7.csv
```

Each case also receives an immutable copy of the authoritative subprocess table:

```text
rcurve_temperature_summary_v9_11.csv
```

## Validation

```bash
python -m py_compile \
  arrhenius_fracture/material_rcurve_audit_v10057.py \
  run_mpz_v10_0_5_7_tip_only_material_rcurve.py \
  run_mpz_v10_0_5_7_tip_only_material_rcurve_fullfield.py

pytest -q \
  tests/test_material_rcurve_audit_v912.py \
  tests/test_mpz_v9_12_runner.py \
  tests/test_v912_paper_transfer_contract.py \
  tests/test_field_snapshots_v912.py \
  tests/test_v10057_v912_gate_and_summary.py
```

## Deliberately not changed

The following scientific audits remain separate and must not be hidden inside a
bookkeeping release:

- anisotropic `J -> K` conversion versus the current isotropic `sqrt(J E')` convention;
- full elastic energy versus tensile-filtered energy in the domain J integral;
- finite starter-notch reference versus sharp-edge-crack LEFM;
- cooperative-cleavage renewal-rate ceiling;
- tip-zoom visualization and emitted-ledger naming in field snapshots.

Those items require explicit benchmark evidence and dedicated point releases.
