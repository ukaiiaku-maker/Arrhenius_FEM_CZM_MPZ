# MPZ v9.10 unified transport/retention search

## Why v9.10 is required

The v9.9 spatial promotion showed that the analytical retention proxy and the
spatial MPZ used different physics.  The analytical proxy treated slow
Peierls–Taylor escape as retained shielding, while the spatial state emitted
mobile dislocations and required a separate legacy Arrhenius trap barrier to
create retained lines.  Consequently, the analytical DBTT R-curve disappeared
in spatial crack growth and the weak-temperature/FCC-like search could not
combine fast Peierls motion with retained shielding.

Version 9.10 removes the independent trap barrier and uses the same mobile and
retained populations in the analytical screen and spatial promotion.

## Unified state equations

```text
dN_m/dt = R_emit - k_enc N_m + k_T N_r - k_esc N_m - k_mrec N_m
dN_r/dt =          k_enc N_m - k_T N_r             - k_rrec N_r
```

with

```text
v_P   = jump_length * lambda_P
k_enc = eta_enc * v_P * sqrt(rho_f)
k_T   = lambda_T_completion
k_esc = v_P / L_MPZ
```

Peierls motion therefore transports mobile dislocations and creates physical
forest encounters.  Taylor completion releases retained dislocations.  A
frozen Peierls branch cannot create retained shielding.

The reference barriers are ordered by construction:

```text
H_T = H_P + Delta H_PT,   Delta H_PT > 0
```

and the full free-energy surfaces are also checked on a common stress and
temperature grid.

## Broad search rather than prior down-selection

Every v9.10 optimization begins from a full Sobol population over the same
25-dimensional domain.  It does not use the v9.8.1 or v9.9 shortlist as an
initial population.  Multiple independent differential-evolution restarts are
run for each response class.

The free variables are:

- seven cleavage EXP-floor parameters;
- seven emission EXP-floor parameters;
- Peierls reference barrier and positive Taylor increment;
- independent Peierls and Taylor activation entropies;
- Taylor correlation density and scale;
- source inventory;
- geometric encounter efficiency;
- retained recovery rate;
- source-refresh length;
- blunting coefficient.

The microscopic attempt frequencies remain fixed at `1e12 s^-1` for Peierls
and `1e11 s^-1` for Taylor.  The independent entropies provide effective
prefactor variation without independently fitting an eleven-decade attempt
frequency range.

## Class handling

The same governing equations and broad bounds are used for all classes.

- `ceramic`: no Peierls-mobility requirement; a small spatial R-curve is
  acceptable.
- `weakT`: FCC-like screening requires process-zone-scale Peierls traversal,
  a weak temperature dependence, a moderate R-curve, and `H_T/H_P >= 2`.
- `DBTT`: acceptance is trend-based and requires the low-to-high temperature
  toughness increase and the high-temperature R-curve to persist during
  repeated crack advance.

## Verification

```bash
conda activate arrhenius-fem-czm
python -m pip install -e ".[dev]"
bash verify_mpz_v9_4.sh
```

The final message should be:

```text
MPZ v9.6 production through v9.10 unified broad-search verification passed.
```

## Global-search smoke test

The smoke search still uses the full 25-dimensional Sobol population; it only
reduces the number of generations and restarts.

```bash
OUTROOT=runs/mpz_v9_10_unified_global_search_smoke_v1 \
TARGET_CLASSES="ceramic weakT DBTT" \
RESTARTS=1 \
DE_MAXITER=2 \
DE_POPSIZE=3 \
LOCAL_MAXITER=20 \
MAX_JOBS=2 \
DK=1.0 \
KMAX=80 \
TARGET_EXTENSION_UM=50 \
DA_UM=5 \
bash run_mpz_v9_10_unified_global_search.sh
```

Each class writes:

- `unified_global_all_candidates.csv`
- `unified_global_accepted.csv`
- `unified_global_shortlist.csv`
- `unified_global_temperature_detail.csv`
- `unified_global_event_detail.csv`
- `unified_global_generation_history.csv`
- `spatial_promotion_manifest.csv`
- `unified_global_summary.json`
- `unified_global_config.json`
- `checkpoints/restart_*.json`

## Production global search

```bash
nohup env \
OUTROOT=runs/mpz_v9_10_unified_global_search_v1 \
TARGET_CLASSES="ceramic weakT DBTT" \
TEMPERATURES="300 700 900 1200" \
RESTARTS=3 \
DE_MAXITER=60 \
DE_POPSIZE=8 \
LOCAL_MAXITER=250 \
MAX_JOBS=2 \
DK=0.5 \
KDOT=0.005 \
KMAX=80 \
TARGET_EXTENSION_UM=500 \
DA_UM=5 \
bash run_mpz_v9_10_unified_global_search.sh \
> runs/mpz_v9_10_unified_global_search_v1.nohup.log 2>&1 &
```

## Spatial promotion smoke test

```bash
OUTROOT=runs/mpz_v9_10_unified_spatial_promotion_smoke_v1 \
MANIFEST_ROOT=runs/mpz_v9_10_unified_global_search_smoke_v1 \
CLASSES="ceramic weakT DBTT" \
MAX_PER_CLASS=1 \
TEMPERATURES="300 700 900 1200" \
TARGET_EXTENSION_UM=50 \
DA_UM=5 \
DK=0.5 \
KMAX=80 \
MPZ_LENGTH_UM=50 \
MPZ_N_BINS=80 \
bash run_mpz_v9_10_unified_spatial_promotion.sh
```

The spatial stage is still a reduced sharp-front MPZ validation.  Only
candidates whose class behavior persists there should proceed to the full 2-D
FEM/CZM solver.
