# MPZ v9.9 barrier continuation and spatial promotion

## Purpose

The v9.8.1 joint optimization found good ceramic and weak-temperature response
basins and a useful DBTT trend, but it also selected large absolute Peierls and
Taylor barriers. Version 9.9 uses those basins as warm starts rather than
restarting a global 17-parameter search.

The workflow has two gates:

1. scale the Peierls and Taylor reference enthalpies together while locally
   reoptimizing entropy, effective prefactors, and state parameters;
2. promote the accepted candidates to the spatial moving-process-zone model and
   test whether the response persists during finite crack extension.

The active production package remains v0.9.6 with the v9.5 spatial state. The
v9.9 independent-entropy state is imported only by the promotion workflow and
is not activated globally.

## Common barrier scaling

For each v9.8.1 basin,

```text
H_P' = s H_P
H_T' = s H_T
```

with default scales

```text
1.0 0.8 0.6 0.4 0.3
```

This preserves both `H_T/H_P` and `H_P <= H_T`. The local continuation varies:

- Peierls and Taylor activation entropies;
- Peierls and Taylor effective prefactors;
- Taylor correlation density and scale;
- mobile fraction;
- finite source inventory;
- recovery rate;
- source-development length;
- blunting coefficient.

The prefactor warm start approximately preserves the zero-stress rate at 700 K
after barrier scaling, but it is not enforced. Powell refinement then optimizes
the coupled response.

## Class-specific acceptance

Only one mechanism-specific requirement is imposed:

- ceramic: slow Peierls motion is permitted;
- weakT: the Peierls traverse number must be at least one over the loading time;
- DBTT: no Peierls/Taylor dominance requirement is imposed.

The DBTT objective uses the low-temperature shelf, high-temperature branch,
R-curve increase, and overall transition trend. It does not require exact
agreement with the provisional point-by-point target curve.

## Spatial promotion

`arrhenius_fracture/moving_process_zone_v99.py` adapts the v9.5 local-density
state to the v9.7 independent-entropy Peierls–Taylor model. The historical
`pt_*_entropy_ratio` configuration slots are interpreted as `S*/k_B` only in
this promotion state. The global package state is unchanged.

The spatial promotion uses the existing sharp-front reduced R-curve driver and
records:

- initiation and plateau toughness;
- R-curve increment and slopes;
- retained and mobile populations;
- source depletion;
- shielding and blunting;
- final Peierls and Taylor rates;
- Peierls process-zone traverse number;
- complete event histories versus crack extension.

A DBTT candidate passes only when the high-minus-low plateau separation and
high-temperature R-curve remain after finite crack growth.

## Verification

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v9_1_full_git
git pull --ff-only origin main
conda activate arrhenius-fem-czm
python -m pip install -e '.[dev]'
bash verify_mpz_v9_4.sh
```

## Continuation smoke test

```bash
OUTROOT=runs/mpz_v9_9_barrier_continuation_smoke_v1 \
INPUT_ROOT=runs/mpz_v9_8_1_joint_response_optimization_v1 \
CLASSES='ceramic weakT DBTT' \
SCALES='1.0 0.6' \
CANDIDATES_PER_CLASS=1 \
LOCAL_MAXITER=30 \
DK=1.0 \
bash run_mpz_v9_9_barrier_continuation.sh
```

## Continuation production run

```bash
nohup env \
OUTROOT=runs/mpz_v9_9_barrier_continuation_v1 \
INPUT_ROOT=runs/mpz_v9_8_1_joint_response_optimization_v1 \
CLASSES='ceramic weakT DBTT' \
SCALES='1.0 0.8 0.6 0.4 0.3' \
CANDIDATES_PER_CLASS=3 \
LOCAL_MAXITER=400 \
DK=0.5 \
bash run_mpz_v9_9_barrier_continuation.sh \
> runs/mpz_v9_9_barrier_continuation_v1.nohup.log 2>&1 &
```

## Spatial smoke test

```bash
OUTROOT=runs/mpz_v9_9_spatial_promotion_smoke_v1 \
MANIFEST=runs/mpz_v9_9_barrier_continuation_v1/spatial_promotion_manifest.csv \
MAX_PER_CLASS=1 \
TARGET_EXTENSION_UM=50 \
DA_UM=5 \
DK=0.5 \
MPZ_LENGTH_UM=50 \
MPZ_N_BINS=80 \
bash run_mpz_v9_9_spatial_promotion.sh
```

## Spatial production run

```bash
nohup env \
OUTROOT=runs/mpz_v9_9_spatial_promotion_v1 \
MANIFEST=runs/mpz_v9_9_barrier_continuation_v1/spatial_promotion_manifest.csv \
MAX_PER_CLASS=2 \
TARGET_EXTENSION_UM=500 \
DA_UM=5 \
DK=0.25 \
KMAX=80 \
MPZ_LENGTH_UM=100 \
MPZ_N_BINS=200 \
bash run_mpz_v9_9_spatial_promotion.sh \
> runs/mpz_v9_9_spatial_promotion_v1.nohup.log 2>&1 &
```

## Outputs

Continuation:

- `barrier_continuation_all.csv`
- `barrier_continuation_accepted.csv`
- `barrier_continuation_temperature_detail.csv`
- `spatial_promotion_manifest.csv`
- `barrier_continuation_summary.json`
- `barrier_continuation_config.json`

Spatial promotion:

- `spatial_promotion_metrics.csv`
- `spatial_promotion_events.csv`
- `spatial_promotion_summary.csv`
- `spatial_promotion_accepted.csv`
- `spatial_promotion_report.json`
- `spatial_promotion_config.json`

Accepted spatial candidates remain reduced-model results. They require final
2-D FEM/CZM validation before being treated as production parameterizations.
