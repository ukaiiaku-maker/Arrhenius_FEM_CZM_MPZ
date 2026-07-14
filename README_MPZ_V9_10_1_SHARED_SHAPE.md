# MPZ v9.10.1 shared EXP-floor shape search

## Motivation

The initial v9.10 global search varied the EXP-floor shape parameters `alpha`
(`exp_a`) and `n` independently for cleavage and emission.  Peierls and Taylor
inherited the emission shape.  This was asymmetric and made the reported
Peierls/Taylor barrier heights appear more fundamental than they are.

Version 9.10.1 implements the simpler working hypothesis:

```text
alpha_c = alpha_e = alpha_P = alpha_T = alpha_shared
n_c     = n_e     = n_P     = n_T     = n_shared
```

The common shape is optimized globally.  The four mechanisms retain separate
barrier heights and thermal coordinates.  Peierls and Taylor still satisfy

```text
H_T = H_P + Delta H_PT
Delta H_PT > 0
```

so `H_P < H_T` is strict.

At fixed applied stress, the EXP-floor shape strongly affects the local free
energy and therefore the effective activation entropy and temperature
sensitivity.  The shared-shape search tests whether one stress-shape family,
with different barrier heights, is sufficient before introducing four
independent shape pairs.

## Search dimension

The shared-shape search contains 23 global variables.  It replaces the four
cleavage/emission shape variables

```text
cleave_exp_a, cleave_exp_n, emit_exp_a, emit_exp_n
```

with

```text
shared_exp_a, shared_exp_n.
```

The decoded candidate records all derived mechanism values explicitly:

```text
cleave_exp_a
emit_exp_a
peierls_exp_a
taylor_exp_a
cleave_exp_n
emit_exp_n
peierls_exp_n
taylor_exp_n
```

so the equality can be audited from each CSV row.

Every restart still starts from a fresh full Sobol population.  No v9.8, v9.9,
or v9.10 shortlist is used to initialize the search.

## Smoke search

```bash
OUTROOT=runs/mpz_v9_10_1_shared_shape_global_search_smoke_v1

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
bash run_mpz_v9_10_1_shared_shape_global_search.sh
```

## Spatial smoke

```bash
OUTROOT=runs/mpz_v9_10_1_shared_shape_spatial_promotion_smoke_v1

rm -rf "$OUTROOT"

env \
OUTROOT="$OUTROOT" \
MANIFEST_ROOT=runs/mpz_v9_10_1_shared_shape_global_search_smoke_v1 \
CLASSES="ceramic weakT DBTT" \
MAX_PER_CLASS=3 \
TEMPERATURES="300 700 900 1200" \
TARGET_EXTENSION_UM=50 \
DA_UM=5 \
DK=0.5 \
KDOT=0.005 \
KMAX=80 \
MPZ_LENGTH_UM=50 \
MPZ_N_BINS=80 \
bash run_mpz_v9_10_1_shared_shape_spatial_promotion.sh
```

## Escalation criterion

Do not immediately add four independent shape pairs.  First test whether the
shared-shape model produces a weak-temperature/FCC-like candidate and a DBTT
candidate whose R-curves survive spatial crack growth.

If a well-converged broad search still fails, the next model should release

```text
(alpha_c, n_c)
(alpha_e, n_e)
(alpha_P, n_P)
(alpha_T, n_T)
```

in one new full-space search.  That expanded search should again start from a
new Sobol population rather than from the shared-shape shortlist.
