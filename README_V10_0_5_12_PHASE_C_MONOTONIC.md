# v10.0.5.12 — four-option Phase-C monotonic FEM/CZM campaign

## Purpose

This release starts the controlled Phase-C monotonic campaign after completion of
0-D, 1-D, and initial 2-D parameter transfer. It adds no fitting and changes no
Arrhenius surface. It integrates the validated 330 micrometer physical
refinement policy into the ordinary progressive FEM/CZM path and exposes the
four exact v9.11.1 response options:

| Option | Candidate | MPZ grid |
|---|---|---|
| `ceramic_primary` | `ceramic_restart02_candidate00` | 100 um / 200 bins |
| `weakT_primary` | `weakT_restart00_candidate00` | 100 um / 200 bins |
| `dbtt_primary` | `DBTT_restart04_candidate03` | 50 um / 80 bins |
| `peak_primary` | `DBTT_restart05_candidate61` | 50 um / 80 bins |

The registry is packaged at
`arrhenius_fracture/data/mpz_v9_11_1/response_registry.json`. Every case writes
an exact one-row CSV, a SHA-256 fingerprint, and a production manifest before the
mechanics solve starts.

## Fixed scientific and mechanics contract

- deterministic expected finite-source emission;
- one straight active front; branching disabled;
- cubic anisotropic elasticity at 45 degrees;
- physical crack increment 5 um;
- 330 um physical refinement radius through initial corridor and remesh paths;
- selected cluster-J outer radius 240 um;
- local-J outer radius 100 um;
- requested tip spacing 2.5 um and grading ratio 1.15;
- reference loading increment `dU=2e-7 m`, `dt=8.4 s`;
- mobile shielding fraction zero; explicit unresolved shielding derives from
  retained state only.

## Validation

```bash
python -m py_compile \
  arrhenius_fracture/mpz_response_registry_v100512.py \
  arrhenius_fracture/mode_i_first_passage_v10_0_5_12_phase_c.py \
  run_v10_0_5_12_phase_c_monotonic.py

pytest -q tests/test_v100512_phase_c.py
```

A command-only matrix audit can be run without FEM:

```bash
python run_v10_0_5_12_phase_c_monotonic.py \
  --mode full \
  --dry-run \
  --outroot runs/v10_0_5_12_phase_c_full_dryrun_v1
```

This must write exactly 40 cases and apply the class-specific MPZ grids.

## Launch sequence

### 1. Current-stack smoke

```bash
MODE=smoke \
MAX_JOBS=1 \
OUTROOT=runs/v10_0_5_12_phase_c_smoke_DBTT_700K_50um_v1 \
bash run_v10_0_5_12_phase_c_monotonic.sh
```

The smoke is `dbtt_primary`, 700 K, 50 um. It must pass the v10.0.5.2
completion manifest, exact parameter fingerprint, and physical-refinement
runtime check.

### 2. Sixteen anchor cases

```bash
MODE=anchors \
MAX_JOBS=2 \
OUTROOT=runs/v10_0_5_12_phase_c_500um_theta45_v1 \
bash run_v10_0_5_12_phase_c_monotonic.sh
```

This runs all four options at 300, 700, 900, and 1200 K to 500 um. The full
launch deliberately uses the same output root so these 16 completed cases are
reused rather than recomputed.

### 3. Complete forty-case Phase-C matrix

```bash
MODE=full \
MAX_JOBS=2 \
OUTROOT=runs/v10_0_5_12_phase_c_500um_theta45_v1 \
bash run_v10_0_5_12_phase_c_monotonic.sh
```

The runner uses `--skip-existing` by default. Reuse is accepted only when the
case-level `.phase_c_complete` marker exists and the completion, parameter, and
refinement manifests still pass.

## Outputs

The campaign root contains:

- `phase_c_registry_snapshot.json`;
- `phase_c_provenance.json`;
- `phase_c_matrix.csv/json`;
- `phase_c_campaign.json`;
- incremental and final `phase_c_summary.csv/json`;
- `phase_c_completion.json`.

Each case contains the exact command, input option/fingerprint, selected
one-row manifest, v10.0.5.2 lifecycle output, the v10.0.5.12 production manifest,
raw step/event histories, R-curve output, snapshots, and a normalized Phase-C
case summary.

The summary separates first passage from propagation resistance and reports
`K_FP`, the 0-50 um median, the 200-300 um median, the late-window median,
`delta_K_R`, normalized R-curve area, completion, and manifest status.

## Promotion rule

Do not launch stochastic long R-curves or fatigue from a partially passing
matrix. Phase C passes only when every requested case is complete or explicitly
classified, exact fingerprints match, physical refinement is verified, and no
case-level lifecycle failure is hidden.
