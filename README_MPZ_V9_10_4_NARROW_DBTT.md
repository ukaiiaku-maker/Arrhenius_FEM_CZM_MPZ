# MPZ v9.10.4: mechanism-first narrow-DBTT campaign

## Purpose

The v9.10.3 DBTT target and acceptance gate rewarded a broad increase from the
low-temperature shelf to the high-temperature branch. It did not measure the
transition width, and the selected candidates used large direct cleavage
`gT`/`sT` terms. v9.10.4 replaces that calibration path with a staged,
mechanism-first campaign.

The transition temperature is not prescribed. Every candidate is evaluated at
300--1100 K in 100 K increments, and every admissible adjacent temperature pair
is tested as the possible DBTT interval.

## Reduced model changes

`arrhenius_fracture/reduced_campaign_front_v9104.py` is a prescribed-K reduced
counterpart of the finalized moving-tip model. It implements:

- unshielded opening stress;
- cleavage stress reduced only by active elastic shielding;
- emission stress reduced only by a local Taylor back stress;
- exact bounded finite-source emission;
- Peierls transport and Taylor trapping/release;
- retained recovery and mobile escape;
- slip-based blunting;
- continuous cleavage action and moving-tip translation;
- crack-advance-only source refresh;
- adaptive kinetic substeps and event localization.

It deliberately excludes the 2-D FEM solve, branching, mesh topology, wake
shielding, geometry/source feedback, and any time-based source recycling.

## Transition gate

For a freely selected adjacent split `[T_j,T_{j+1}]`, the first-passage gate
requires approximately:

- high/low median toughness ratio >= 2;
- robust min(high)/max(low) ratio >= 1.8;
- at least 75% of the shelf change in the selected adjacent jump;
- low-shelf span <= 15%;
- high-shelf span <= 20%;
- matched plasticity-off high/low ratio <= 1.25;
- at least 60% of the full shelf jump attributable to the full-minus-
  plasticity-off response.

The default phase-A search fixes the direct cleavage temperature slopes to zero.
An optional `narrow` mode permits small slopes with strong regularization only
after the fixed-zero feasibility result has been assessed.

## Staged campaign

### Stage 0: audit historical candidates

```bash
MANIFEST=/path/to/historical/DBTT/spatial_promotion_manifest.csv \
OUT=runs/mpz_v9_10_4_current_dbtt_audit_v1 \
bash scripts/run_mpz_v9_10_4_stage0_audit.sh
```

This evaluates the old and new reduced models on all nine temperatures.

### Stage 1: reduced first-passage search

```bash
OUT=runs/mpz_v9_10_4_narrow_dbtt_first_passage_fixed_zero_v1 \
CLEAVAGE_SLOPE_MODE=fixed_zero \
RESTARTS=4 \
DE_MAXITER=80 \
DE_POPSIZE=8 \
LOCAL_MAXITER=250 \
bash scripts/run_mpz_v9_10_4_first_passage_search.sh
```

Do not begin with `CLEAVAGE_SLOPE_MODE=narrow`. Use it only if the fixed-zero
search finds no feasible mechanism-driven candidate.

### Stage 2: short-growth refinement

```bash
INPUT_MANIFEST=runs/mpz_v9_10_4_narrow_dbtt_first_passage_fixed_zero_v1/short_growth_promotion_manifest.csv \
OUT=runs/mpz_v9_10_4_narrow_dbtt_short_growth_100um_v1 \
bash scripts/run_mpz_v9_10_4_short_growth.sh
```

This refines and validates candidates through 100 um of crack extension.
Initiation and plateau transitions must select the same temperature interval.

### Stage 3: long-growth refinement

```bash
INPUT_MANIFEST=runs/mpz_v9_10_4_narrow_dbtt_short_growth_100um_v1/narrow_dbtt_short_growth_promotion_manifest.csv \
OUT=runs/mpz_v9_10_4_narrow_dbtt_long_growth_500um_v1 \
bash scripts/run_mpz_v9_10_4_long_growth.sh
```

This validates the surviving candidates through 500 um.

### Stage 4: prepare the 2-D validation matrix

```bash
INPUT_MANIFEST=runs/mpz_v9_10_4_narrow_dbtt_long_growth_500um_v1/narrow_dbtt_long_growth_promotion_manifest.csv \
OUT=runs/mpz_v9_10_4_2d_validation_manifest_v1.csv \
bash scripts/run_mpz_v9_10_4_prepare_2d_validation.sh
```

For each selected candidate, the generated matrix includes one temperature below
the transition, the two 100 K bracket temperatures, and one temperature above.
The required 2-D ablations are `full`, `plasticity_off`, `backstress_off`,
`shielding_off`, and `blunting_off`.

## Promotion rule

No parameter set replaces the active DBTT manifest until it passes:

1. fixed-zero or strongly regularized first-passage gate;
2. 100 um short-growth gate;
3. 500 um long-growth gate;
4. fine 25--50 K transition-width check around the selected bracket;
5. 2-D full and matched-ablation validation.

If the fixed-zero and narrow-slope searches cannot produce the required factor-
of-two step, the correct conclusion is that the current mechanism set cannot
support the requested DBTT without a strong intrinsic cleavage-temperature
term. The optimizer must not recreate the historical large activation entropy
as an undocumented shortcut.
