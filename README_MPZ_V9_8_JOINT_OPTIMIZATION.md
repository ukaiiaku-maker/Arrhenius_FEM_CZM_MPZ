# MPZ v9.8 joint response optimization

## Purpose

The v9.7 Peierls–Taylor calibration showed that the plasticity barriers cannot
be identified reliably from a single bulk flow-stress target.  Version 9.8
therefore optimizes cleavage, emission, Peierls, Taylor, and reduced
process-zone parameters together against the requested fracture-response
class.

The current implementation is analytical fidelity 0.  It locates admissible
parameter basins efficiently before the leading diverse candidates are
promoted to the spatial moving-process-zone and 2-D FEM/CZM solvers.

## Optimizer

Each response class is optimized independently with identical governing
physics and bounds:

1. select diverse first-passage atlas seeds, including the prior canonical row;
2. initialize a global population with half local perturbations and half Sobol
   coverage;
3. run SciPy differential evolution;
4. refine the best point in each basin with bounded adaptive Nelder–Mead;
5. retain a score-and-distance diverse shortlist rather than only one optimum.

No assumption is made about whether Peierls or Taylor must control a named
response class.

## Joint parameter vector

The initial v9.8 vector contains 17 variables:

- cleavage `G00`, `gT`, and critical stress;
- emission `G00`, `gT`, and critical stress;
- absolute Peierls reference barrier `H_P`;
- non-negative Taylor increment `Delta H_PT`, with `H_T=H_P+Delta H_PT`;
- independent Peierls and Taylor activation entropies;
- Taylor correlation density and correlation scale;
- mobile fraction;
- finite source-site inventory;
- recovery rate;
- source refresh/development length;
- blunting coefficient.

The EXP-floor shape exponents and floor fractions remain fixed within each
atlas-seeded basin.  This limits the global dimension while preserving several
distinct shape families.  They can be released in a later local refinement if
the retained basins require it.

## Physical constraints

Every objective evaluation requires:

```text
G_P(sigma,T) <= G_T(sigma,T)
```

on the same resolved-stress grid, positive raw Peierls and Taylor barriers,
exact detailed balance, and the uncapped v9.6/v9.7 kinetics.  The optimization
contains no constitutive density, stress, hit-order, mobile-density,
jump-length, or rate cap.

The response labels are output targets only.  No bottleneck or mechanism-family
rule is imposed.

## Fidelity-0 response calculation

For each temperature the objective evaluates:

- exact monotonic cleavage first passage on a K grid;
- emission exposure accumulated before cleavage;
- self-consistent emission/escape/recovery retention using the uncapped
  detailed-balance Peierls–Taylor model;
- physical line-kernel shielding;
- slip-based blunting;
- an exponential development length based on the source-refresh length;
- initial toughness, plateau toughness, R-curve increment, early rise, and
  plateau slope.

The previous `N_sat` and back-stress fits are not used as production caps or
hard optimizer parameters.

## Smoke run

```bash
OUTROOT=runs/mpz_v9_8_joint_response_smoke_v1 \
TARGET_CLASSES="ceramic weakT DBTT" \
SEED_COUNT=2 \
SEED_POOL_SIZE=30 \
DE_MAXITER=3 \
DE_POPSIZE=4 \
LOCAL_MAXITER=30 \
MAX_JOBS=2 \
DK=1.0 \
bash run_mpz_v9_8_joint_response_optimization.sh
```

## Production analytical search

```bash
nohup env \
OUTROOT=runs/mpz_v9_8_joint_response_optimization_v1 \
TARGET_CLASSES="ceramic weakT DBTT" \
TEMPERATURES="300 700 900 1200" \
SEED_COUNT=6 \
SEED_POOL_SIZE=120 \
DE_MAXITER=40 \
DE_POPSIZE=8 \
LOCAL_MAXITER=300 \
MAX_JOBS=2 \
DK=0.5 \
bash run_mpz_v9_8_joint_response_optimization.sh \
> runs/mpz_v9_8_joint_response_optimization_v1.nohup.log 2>&1 &
```

Completed basin checkpoints are reused.  Incomplete generation checkpoints are
removed automatically and that basin restarts from its deterministic seed.

## Outputs for each class

- `selected_atlas_seeds.csv`
- `joint_response_basin_results.csv`
- `joint_response_shortlist.csv`
- `joint_response_temperature_detail.csv`
- `joint_response_generation_history.csv`
- `joint_response_promotion_manifest.csv`
- `joint_response_optimization_config.json`
- `joint_response_optimization_summary.json`
- `checkpoints/basin_*.json`

The promotion manifest is not a final calibration.  The next gate is a spatial
moving-PZ run for the leading diverse candidates, followed by low-, transition-,
and high-temperature FEM/CZM validation.
