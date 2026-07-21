# MPZ v9.12 physics-informed active-learning campaign

## Why the original Sobol space is too large

The original profile has 19 continuously varied coordinates:

- 4 emission-surface coordinates;
- 4 Peierls coordinates;
- 4 Taylor coordinates;
- 2 Taylor-correlation coordinates;
- source density and refresh length;
- 3 recovery coordinates.

A 4096-point Sobol design in 19 dimensions provides only
`4096^(1/19) = 1.55` effective levels per coordinate. Even a two-level full
factorial would require 524288 cases. The initial design is therefore a broad
exploration set, not a thorough characterization of the 19-D domain.

The Arrhenius competition is more naturally represented by

```text
log(lambda_i/lambda_j)
  = log(nu_i/nu_j) - [G_i-G_j]/(kB T)
```

and by residence-time groups such as

```text
Pi_store   = k_encounter * t_reference
Pi_release = lambda_Taylor * t_reference
Pi_recover = k_recovery * t_reference.
```

Raw barrier ratios may correlate with response, but barrier differences divided
by `kB*T` are the direct rate-ratio coordinates. Absolute barriers still set the
stress/time scale.

## State-focused profile

`mpz_v9_12_emergent_gnd_search_bounds_state_focused.json` reduces the first
campaign from 19 to 11 continuous coordinates by holding the EXP-floor shape
exponents, Taylor correlation scale, and recovery prefactor at their seed-row
values. Distinct shallow/base/steep barrier shapes should be treated as separate
seed families rather than diluted into one high-dimensional design.

The 11 varied coordinates retain:

- emission level and temperature slope;
- Peierls level and activation entropy;
- Taylor level and activation entropy;
- Taylor correlation density;
- physical source density;
- source-refresh length;
- recovery level and activation entropy.

## Numerical stiffness correction

The original adaptive timestep used the largest raw Arrhenius attempt/event rate
even though source depletion, mobile-retained exchange, and first-order recovery
were already updated by exact bounded exponentials. Near-prefactor rates could
therefore request more than one million unnecessary substeps.

`emergent_gnd_state_v912_stiff.py` now:

- advances those local stiff reactions with their existing exponential updates;
- limits constitutive feedback re-evaluation using a common physical interval;
- applies a CFL restriction only to explicit 1-D transport;
- leaves the homogeneous 0-D gate free of a nonexistent transport constraint.

Environment controls are:

```text
MPZ_V912_MAX_FEEDBACK_SUBSTEP_S   default 0.1 s
MPZ_V912_TRANSPORT_CFL            default 0.25
```

## Recommended multi-fidelity loop

1. Generate 2048-4096 initial state-focused Sobol candidates.
2. Run the dense-temperature 0-D screen.
3. Build physics-informed ML features.
4. Train the ensemble surrogate.
5. Generate a much larger unevaluated candidate pool.
6. Select 256-512 candidates using predicted utility, uncertainty, feasibility,
   and feature-space diversity.
7. Run those candidates in 0-D and retrain for several rounds.
8. Promote a diverse 256-512 candidate subset to 1-D.
9. Train a second surrogate for the 0-D-to-1-D correction.
10. Promote approximately 20-50 candidates to the exact 2-D implementation.

Nelder-Mead is not recommended for this discontinuous, multimodal, constrained
problem. Particle swarm or cross-entropy can be useful as proposal generators,
but the active-learning surrogate preserves unsuccessful and unresolved cases
as information and supports explicit diversity and multi-fidelity promotion.

## Install ML dependencies

```bash
python -m pip install -e '.[ml]'
```

## Build the initial state-focused design

```bash
python -u scripts/generate_mpz_v9_12_sobol_candidates.py \
  --base-registry mpz_v9_11_dbtt_option_rows.csv \
  --base-candidate-id DBTT_restart04_candidate03 \
  --bounds-json mpz_v9_12_emergent_gnd_search_bounds_state_focused.json \
  --n 4096 \
  --seed 912 \
  --prefix v912_state_focused_r0 \
  --out candidates/v9_12_state_focused_r0_4096.csv
```

## Build the training table after the 0-D run

```bash
python -u scripts/build_mpz_v9_12_ml_table.py \
  --candidate-registry candidates/v9_12_state_focused_r0_4096.csv \
  --ranking-csv runs/v9_12_state_focused_r0_0d/ranking.csv \
  --physics-json mpz_v9_12_emergent_gnd_common_physics.json \
  --bounds-json mpz_v9_12_emergent_gnd_search_bounds_state_focused.json \
  --temperatures 300 500 700 900 1100 1200 \
  --K-values 15 25 35 \
  --reference-time-s 8.4 \
  --out ml/v9_12_state_focused_r0_training.csv
```

The table includes raw coordinates and engineered features at each sampled
stress-temperature point:

- `(G_emit-G_cleave)/(kB*T)`;
- `(G_Peierls-G_emit)/(kB*T)`;
- `(G_Taylor-G_Peierls)/(kB*T)`;
- emission rate and Peierls velocity;
- `Pi_store`, `Pi_release`, and `Pi_recovery`;
- physical source inventory;
- refresh length divided by MPZ length;
- Taylor correlation-order increment at the forest floor.

## Train the surrogate

```bash
python -u scripts/train_mpz_v9_12_surrogate.py \
  --ml-table ml/v9_12_state_focused_r0_training.csv \
  --trees 500 \
  --seed 912 \
  --out-model ml/v9_12_state_focused_r0_surrogate.joblib \
  --out-dir ml/v9_12_state_focused_r0_fit
```

The trainer fits separate ensemble models for:

- numerical completion/feasibility;
- formal campaign pass probability when both classes are present;
- campaign score;
- `Delta_K_micro` amplitude;
- largest-jump localization;
- transition width;
- maximum shielding;
- maximum signed-GND content.

## Propose an active batch

First generate a large unevaluated Sobol pool and build its feature table without
a ranking CSV. Then select a diverse batch:

```bash
python -u scripts/propose_mpz_v9_12_active_batch.py \
  --model ml/v9_12_state_focused_r0_surrogate.joblib \
  --pool-ml-table ml/v9_12_state_focused_pool_100000_features.csv \
  --pool-registry candidates/v9_12_state_focused_pool_100000.csv \
  --batch-size 512 \
  --beta 1.5 \
  --diversity-weight 0.35 \
  --out candidates/v9_12_state_focused_r1_active_512.csv
```

The acquisition rewards predicted response quality and ensemble uncertainty,
multiplies by completion and pass probabilities, and then applies greedy
feature-space diversity so one narrow predicted optimum does not consume the
whole batch.

## Interpretation

Feature importance is a screening diagnostic, not proof of causality. Confirm
important coordinates by controlled 0-D and 1-D ablations. In particular, use
the surrogate to identify whether favorable response is organized primarily by
barrier differences, source inventory, refresh length, Taylor release,
recovery, or their interactions; then rerun selected slices with the remaining
coordinates fixed.
