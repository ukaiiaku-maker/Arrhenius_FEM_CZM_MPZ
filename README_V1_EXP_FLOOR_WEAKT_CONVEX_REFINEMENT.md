# V1 EXP-floor weak-T convex refinement

This calibration treats the prior weak-T class differently from the peak class.
The old weak-T interior jump is **not** used as a fitting target. Instead, the
old 300 and 1200 K toughness values anchor a smooth weakly temperature-dependent
convex reference

\[
K_{\rm ref}(T)=K_\infty + A\exp[-(T-T_0)/\tau],
\]

with `TAU_K=450 K` by default. The optimizer rewards:

- small error relative to this smooth weak-T reference;
- low total temperature variation;
- small adjacent-temperature jumps;
- nondecreasing slope with temperature (convexity);
- weakly decreasing rather than oscillatory behavior;
- smoothness and correct overall toughness scale.

Convexity/topology metrics are evaluated on the 100 K grid used by the PF
comparison, which avoids treating fine `dK` first-passage quantization as a
physical loss of convexity.

The search reuses two existing populations:

1. `runs/v1_exp_floor_four_class_tuning` — original 24,576-candidate broad screen;
2. `runs/v1_exp_floor_peak_expanded_search` — continuous global search over both
   emission and cleavage EXP-floor barriers.

It then performs three local refinement generations with continuous mutation of
both barrier channels and the state-coupling parameters.

The improved peak result is treated as locked. The harness includes a copy of
`locked_peak_expanded_recommendation.csv` and
`locked_peak_expanded_curve_dense.csv`; the runtime output also copies the peak
recommendation from `PEAK_RUN` when available.

## Smoke test

```bash
SMOKE=1 \
SOURCE_RUN=runs/v1_exp_floor_four_class_tuning \
PEAK_RUN=runs/v1_exp_floor_peak_expanded_search \
OUT=runs/v1_exp_floor_weakT_convex_smoke \
bash run_v1_exp_floor_weakT_convex_refinement.sh
```

## Production run

```bash
SOURCE_RUN=runs/v1_exp_floor_four_class_tuning \
PEAK_RUN=runs/v1_exp_floor_peak_expanded_search \
OUT=runs/v1_exp_floor_weakT_convex_refinement \
N_SEEDS=64 \
GEN1_PERTURB=64 \
GEN2_SEEDS=36 \
GEN2_PERTURB=72 \
GEN3_SEEDS=20 \
GEN3_PERTURB=64 \
FINALISTS=12 \
RESUME=1 \
bash run_v1_exp_floor_weakT_convex_refinement.sh
```

## Main outputs

- `recommended_weakT_convex.csv`
- `recommended_weakT_curve_dense.csv`
- `weakT_convex_fit.png`
- `saved_peak_recommendation.csv`
- `saved_peak_curve_dense.csv`
- `saved_regime_recommendations.csv`
- generation candidate/score tables and compressed `Kc` arrays

`TAU_K` is configurable. Larger values make the smooth reference closer to a
weak linear decrease; smaller values make the reference more strongly convex.
