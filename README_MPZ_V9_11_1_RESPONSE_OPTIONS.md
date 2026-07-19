# MPZ v9.11.1 response-mechanism options

## Purpose

The 2-D validation results show that a ductile-to-brittle transition temperature
(DBTT) response is not evidence for one unique microscopic mechanism. Within
the same Arrhenius first-passage/MPZ framework, a DBTT-like macroscopic trend
can arise from:

1. intrinsic cleavage first-passage kinetics with negligible retained shielding;
2. a sharp onset of a retained moving process zone and sustained crack-tip shielding;
3. a broad mixed response in which shielding is substantial but much of the raw
   `K_J(delta a)` curve remains controlled by intrinsic kinetics and specimen
   geometry; or
4. nonmonotonic competition that produces a peak-type toughness response.

The repository therefore preserves named response options rather than replacing
all alternatives with one nominal DBTT parameterization.

## Frozen roles

| Option | Candidate | Role |
|---|---|---|
| `dbtt_primary` | `DBTT_restart04_candidate03` | Main production DBTT; sharp high-temperature retained-zone onset |
| `peak_primary` | `DBTT_restart05_candidate61` | Preserved peak-type response |
| `dbtt_broad_shielding` | `DBTT_restart01_candidate68` | Broad sustained-shielding alternate |
| `dbtt_intrinsic_control` | `DBTT_restart00_candidate103` | Negligible-shielding intrinsic first-passage control |
| `dbtt_moderate_shielding_reference` | `DBTT_restart00_candidate04` | Intermediate moderate-shielding reference |
| `weakT_primary` | `weakT_restart00_candidate00` | Weak-temperature reference |

The main paper calculations should use `dbtt_primary`. A shorter mechanism
section can compare `dbtt_primary`, `dbtt_intrinsic_control`,
`dbtt_broad_shielding`, and `peak_primary` to demonstrate non-uniqueness of the
underlying physics.

## Interpretation requirements

The raw 2-D `K_J(delta a)` histories contain a strong common specimen/loading
trajectory. A large value of

```text
median K_J over 200-500 um - K_init
```

must not be labeled a material R-curve without additional controls. Report at
least:

- `K_init`;
- late-growth slope or endpoint relative to initiation;
- retained MPZ population sampled immediately before renewal;
- active `K_shield`;
- comparison with an intrinsic low-shielding control; and
- an initiation-matched or otherwise justified geometry/control subtraction.

This distinction is central to the paper argument: DBTT-like behavior can exist
with or without shielding, and an apparent rising `K_J(delta a)` curve can exist
with or without a plasticity-generated material R-curve.

## Files

- `mpz_v9_11_response_options.json`: scientific roles and recommended sweeps.
- `mpz_v9_11_dbtt_option_rows.csv`: exact DBTT parameter rows retained from the
  validated top-five package.
- `prepare_mpz_v9_11_response_option.py`: builds an isolated runnable parameter
  root without altering the canonical root.
- `run_mpz_v9_11_response_option.sh`: launches one named option with explicit
  start, heartbeat, final-case, and completion reporting.

## Examples

List the available keys:

```bash
python - <<'PY'
import json
from pathlib import Path
r = json.loads(Path('mpz_v9_11_response_options.json').read_text())
for key, value in r['options'].items():
    print(f"{key:36s} {value['candidate_id']:32s} {value['role']}")
PY
```

Run the primary DBTT transition refinement:

```bash
OPTION=dbtt_primary \
TEMPS="950 1000 1050 1100 1150" \
STEPS=25000 \
bash run_mpz_v9_11_response_option.sh
```

Run the peak-class refinement:

```bash
OPTION=peak_primary \
TEMPS="800 850 900 950 1000 1050 1100" \
STEPS=25000 \
bash run_mpz_v9_11_response_option.sh
```

Run the intrinsic control:

```bash
OPTION=dbtt_intrinsic_control \
TEMPS="700 900 1000 1100 1200" \
STEPS=25000 \
bash run_mpz_v9_11_response_option.sh
```

Every materialized root contains `response_option_selection.json`; every run
prints `CAMPAIGN_START`, periodic `HEARTBEAT` records, final case statuses, and
`CAMPAIGN_COMPLETE` or `CAMPAIGN_INCOMPLETE_OR_FAILED`.
