# v10.0.5.4 VHCF first-passage correction

This point release retains the audited v10.0.5.3 FEM/CZM/MPZ constitutive and
cohesive lifecycle. It corrects campaign-level integration and audit behavior
needed for very-high-cycle calculations.

## Corrections

- The physical cycle horizon may be as large as `1e14` cycles.
- `max_block_cycles` defaults to `inf`; the remaining physical horizon and the
  hazard/state-increment limits determine each cycle jump.
- Runtime counters verify that the live tensor-resolved MPZ engine supplies the
  one-cycle predictor. Any legacy scalar-predictor call fails the invocation.
- Exhaustion of `MAX_BLOCKS` is reported as `right_censored_max_blocks`, not as a
  physically complete fatigue calculation.
- Tip-only VHCF runs default to maximum-load FEM plus cycle-integrated local
  hazards. Phase-resolved cyclic FEM remains available with
  `RESOLVE_CYCLIC_MECHANICS=1` for cases with evolving bulk cyclic mechanics.
- `fatigue_block_diagnostics_v10_0_5_4.csv` distinguishes accepted block
  increments from per-cycle predictor quantities.

## Validation

```bash
python -m py_compile \
  arrhenius_fracture/mode_i_first_passage_v10_0_5_4_vhcf.py \
  run_v10_0_5_4_vhcf_delta_sigma.py \
  run_v10_0_5_4_vhcf_delta_sigma_compat.py

pytest -q \
  tests/test_v10053_audited_source_preflight.py \
  tests/test_kinetic_fatigue_v10053.py \
  tests/test_v10054_vhcf_first_passage.py
```

## One-point transition smoke test

This test is intended to pass through the finite source-depletion transient and
then demonstrate a large hazard-limited jump toward the physical cycle horizon.

```bash
MODE=smoke \
MATERIAL=DBTT \
TEMPERATURES="700" \
DELTA_SIGMA_MPA="350" \
R=0.1 \
FREQUENCY_HZ=1000 \
CYCLES_MAX=1e5 \
MAX_BLOCK_CYCLES=inf \
MAX_BLOCKS=500 \
TARGET_EXTENSION_UM=5 \
RESOLVE_CYCLIC_MECHANICS=0 \
OUTROOT=runs/v10_0_5_4_DBTT_700K_350MPa_vhcf_transition_smoke_v1 \
bash run_v10_0_5_4_vhcf_delta_sigma.sh
```

The case is accepted as physically complete only if it reaches the cycle
horizon or target extension. Check:

- `run_completion_v10_0_5_4_vhcf.json`
- `progressive_fatigue_v10_0_5_4_vhcf.json`
- `K_vs_delta_sigma.csv`
- `fatigue_block_diagnostics_v10_0_5_4.csv`

The VHCF audit must report at least one authoritative predictor call, zero
legacy predictor calls, and either `cycle_horizon` or `target_extension` as the
termination condition.
