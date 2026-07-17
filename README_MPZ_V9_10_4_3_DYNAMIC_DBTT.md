# MPZ v9.10.4.3: Dynamic Narrow-DBTT Campaign

## Temperature protocol

The DBTT temperature is not prescribed.

1. **Coarse first passage:** evaluate every candidate at 300, 400, ..., 1100 K.
2. **Bracket selection:** select the candidate's best adjacent 100 K transition interval.
3. **Refined first passage:** place four total points across that interval using `linspace(Tlow, Thigh, 4)`.
4. **Shelf retention:** keep two broad low-shelf and two broad high-shelf anchors from the coarse grid.
5. **Short growth:** use the same candidate-specific schedule to 100 micrometres.
6. **Long growth:** use the same candidate-specific schedule to 500 micrometres.
7. **2-D validation:** evaluate one shelf point below, all four refined transition points, and one shelf point above.

The coarse stage alone chooses the transition location. Refined stages hold that bracket fixed while optimizing the constitutive parameters.

## Numerical correction

v9.10.4.1 limited timesteps using microscopic rates whenever any nonzero source or dislocation population remained. Continuous source refresh can leave a trace population, so extremely large rates forced hundreds of thousands of steps even when the maximum possible state change was negligible.

v9.10.4.3 limits the **absolute count changed per step**, normalized by the finite source capacity. Exact finite-source depletion and exact exchange updates remain unchanged. A trace population smaller than the permitted count change may be exhausted in one step.

## Stage sequence

### A. Re-run the historical candidate diagnostic

Use the existing Stage 0 audit with the current one-row DBTT manifest. A three-temperature run is diagnostic only and reports an invalid transition reason rather than attempting shelf scoring.

### B. Coarse first-passage search

Use `scripts/run_mpz_v9_10_4_first_passage_search.sh` with the full default 300--1100 K grid.

### C. Build candidate-specific schedules

```bash
INPUT_MANIFEST=<coarse-output>/short_growth_promotion_manifest.csv \
OUT=<coarse-output>/dynamic_refinement_manifest.csv \
bash scripts/run_mpz_v9_10_4_3_prepare_dynamic_refinement.sh
```

### D. Refined first passage

```bash
INPUT_MANIFEST=<coarse-output>/dynamic_refinement_manifest.csv \
OUT=runs/mpz_v9_10_4_3_refined_first_passage_v1 \
bash scripts/run_mpz_v9_10_4_3_refined_first_passage.sh
```

### E. Short growth

```bash
INPUT_MANIFEST=runs/mpz_v9_10_4_3_refined_first_passage_v1/short_growth_promotion_manifest.csv \
STAGE=short \
OUT=runs/mpz_v9_10_4_3_short_growth_100um_v1 \
bash scripts/run_mpz_v9_10_4_3_growth.sh
```

### F. Long growth

```bash
INPUT_MANIFEST=runs/mpz_v9_10_4_3_short_growth_100um_v1/narrow_dbtt_v91043_short_growth_promotion_manifest.csv \
STAGE=long \
OUT=runs/mpz_v9_10_4_3_long_growth_500um_v1 \
bash scripts/run_mpz_v9_10_4_3_growth.sh
```

### G. 2-D validation manifest

Use `prepare_mpz_v9_10_4_2d_validation.py` on the long-growth promotion manifest.

## Acceptance logic

The candidate must retain:

- at least a factor-of-two high/low toughness ratio;
- robust separation between the shelves;
- a dominant narrow jump inside its selected coarse 100 K bracket;
- flat low- and high-temperature shelves;
- a comparatively flat plasticity-off response;
- initiation and plateau transitions in the same candidate-specific bracket;
- a non-collapsing high-temperature R-curve.
