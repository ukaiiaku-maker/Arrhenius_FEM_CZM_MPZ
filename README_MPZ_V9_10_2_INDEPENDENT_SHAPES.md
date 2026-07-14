# MPZ v9.10.2 independent four-barrier EXP-floor shapes

## Motivation

The intended v9.10 hypothesis was

```text
alpha_e = alpha_P = alpha_T
n_e     = n_P     = n_T
```

with a separate cleavage shape.  The v9.10.1 audit instead imposed one common
shape on cleavage, emission, Peierls, and Taylor.  That audit produced stable
ceramic responses and spatially persistent DBTT candidates, but still produced
no accepted weak-temperature/FCC-like candidate.  Version 9.10.2 therefore
releases the complete shape space:

```text
(alpha_c, n_c)
(alpha_e, n_e)
(alpha_P, n_P)
(alpha_T, n_T)
```

The four barrier heights, thermal slopes/activation entropies, critical-stress
coordinates, and process-zone state parameters remain jointly optimized.

## Physics retained

The v9.10 unified mobile/retained closure is unchanged:

```text
v_P   = jump_length * lambda_P
k_enc = eta_enc * v_P * sqrt(rho_f)
k_T   = lambda_T_completion
k_esc = v_P / L_MPZ
```

There is no independent Arrhenius trapping barrier.  Peierls transport creates
forest encounters, and Taylor completion releases retained dislocations.

The search retains strict reference ordering

```text
H_T = H_P + Delta H_PT,  Delta H_PT > 0
```

and directly rejects candidates whose full same-stress free-energy surfaces
violate `G_P(sigma,T) <= G_T(sigma,T)` over the audit grid.

## Search size

The global vector contains 29 coordinates.  Every restart begins from a fresh
Sobol population over the complete common bounds and does not use a previous
shortlist.

With the production defaults

```text
RESTARTS=3
DE_POPSIZE=6
DE_MAXITER=45
```

SciPy rounds the 29-by-6 Sobol population to 256 members.  The approximate
global evaluation count is therefore

```text
256 * (45 + 1) * 3 = 35,328 evaluations per response class
```

before the local Powell refinements.  This deliberately matches the scale of
the earlier approximately 35,000-point searches while covering the larger
four-shape phase space.

## Smoke search

```bash
OUTROOT=runs/mpz_v9_10_2_independent_shape_global_search_smoke_v1
rm -rf "$OUTROOT"

env \
OUTROOT="$OUTROOT" \
TARGET_CLASSES="ceramic weakT DBTT" \
TEMPERATURES="300 700 900 1200" \
RESTARTS=1 \
DE_MAXITER=2 \
DE_POPSIZE=3 \
LOCAL_MAXITER=20 \
MAX_JOBS=2 \
DK=1.0 \
KDOT=0.005 \
KMAX=80 \
TARGET_EXTENSION_UM=50 \
DA_UM=5 \
bash run_mpz_v9_10_2_independent_shape_global_search.sh
```

## Production search

```bash
OUTROOT=runs/mpz_v9_10_2_independent_shape_global_search_v1
mkdir -p "$OUTROOT"

nohup env \
OUTROOT="$OUTROOT" \
TARGET_CLASSES="ceramic weakT DBTT" \
TEMPERATURES="300 700 900 1200" \
RESTARTS=3 \
DE_MAXITER=45 \
DE_POPSIZE=6 \
LOCAL_MAXITER=250 \
MAX_JOBS=2 \
DK=0.5 \
KDOT=0.005 \
KMAX=80 \
TARGET_EXTENSION_UM=500 \
DA_UM=5 \
bash run_mpz_v9_10_2_independent_shape_global_search.sh \
> "${OUTROOT}.nohup.log" 2>&1 &
```

## Spatial promotion

```bash
OUTROOT=runs/mpz_v9_10_2_independent_shape_spatial_promotion_v1
rm -rf "$OUTROOT"

env \
OUTROOT="$OUTROOT" \
MANIFEST_ROOT=runs/mpz_v9_10_2_independent_shape_global_search_v1 \
CLASSES="ceramic weakT DBTT" \
MAX_PER_CLASS=5 \
TEMPERATURES="300 700 900 1200" \
TARGET_EXTENSION_UM=500 \
DA_UM=5 \
DK=0.25 \
KDOT=0.005 \
KMAX=80 \
MPZ_LENGTH_UM=100 \
MPZ_N_BINS=200 \
bash run_mpz_v9_10_2_independent_shape_spatial_promotion.sh
```

The search and promotion workflows remain calibration gates.  The package-level
active state is still v9.5/v9.6 until a candidate passes spatial crack-growth
persistence and subsequent 2-D FEM/CZM validation.
