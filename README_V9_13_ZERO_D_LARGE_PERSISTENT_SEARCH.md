# v9.13 persistent-site large zero-D DBTT search

## Purpose

This campaign replaces the legacy v9.10 zero-dimensional optimizer for DBTT
candidate generation. The legacy fidelity used a finite source inventory,
source depletion, crack-advance source refresh, and explicit retained recovery.
Those mechanisms are inactive in the calibrated v9.13 one-dimensional transfer
and therefore must not define the next candidate pool.

The new search follows the active v9.13 contract:

- persistent areal nucleation sites;
- no source depletion;
- no crack-advance source refresh;
- no explicit mobile or retained recovery;
- density-dependent source-front width;
- accumulated-slip crack-tip blunting;
- implicit backstress-limited emission;
- current independent cleavage, emission, Peierls, and Taylor EXP-floor shapes;
- the calibrated v10.2.22 common physics and stochastic loading map.

It remains a zero-dimensional screening model. Promotion means only that a
candidate should be evaluated by the full autonomous one-dimensional R-curve
model.

## Two screening fidelities

### Vectorized proxy

A scrambled Sobol population is evaluated with an analytical persistent-state
proxy. It uses the exact current barrier surfaces and common physical
normalizations to estimate the neutral cleavage load, persistent-site emission,
mechanical blocking density, front-width contraction, and blunting response.
The default population is 262,144 candidates.

### Event-clock zero-D replay

The best 4,096 proxy candidates are replayed with a scalar state model using the
actual loading-map geometry factors, stochastic threshold actions, event
lengths, displacement rate, multi-hit cleavage clock, implicit emission root,
and mean-field crack-tip translation. The default output is a diverse registry
of 512 candidates for full one-dimensional screening.

## Search objective

The temperature grid extends through the previously observed high-temperature
re-hardening branch:

```text
700 800 900 950 1000 1050 1100 1200 1300 1400 K
```

The search records separately:

- internal peak temperature;
- two-sided local prominence;
- post-peak drop;
- maximum high-temperature rebound;
- maximum backstress;
- maximum tip radius;
- minimum effective front width.

Default zero-D promotion gates are:

```text
two-sided prominence >= 5 MPa sqrt(m)
post-peak drop       >= 5 MPa sqrt(m)
high-T rebound       <= 3 MPa sqrt(m)
internal peak between 850 and 1100 K
```

The final selection is diversity-aware in all 26 variable active coordinates;
it is not simply the lowest-objective cluster around one anchor.

## Default production run

The required long stochastic map is:

```text
runs/v9_13_long_map_exponential_110um_v2/
  v10_2_22_long_rcurve_loading_map_exponential_110um.json
```

Run:

```bash
OUT=runs/v9_13_zero_d_large_persistent_search_v1 \
SAMPLES=262144 \
EXACT_COUNT=4096 \
PROMOTE_COUNT=512 \
MAX_JOBS=4 \
bash scripts/run_v913_zero_d_large_search.sh
```

The campaign is resumable. Vectorized proxy batches are written under
`proxy_batches/`, and each exact zero-D candidate is written independently under
`exact_cases/`.

## Important outputs

```text
run_contract.json
progress.json
proxy_summary.json
proxy_exact_input.csv
zero_d_ranked_candidates.csv
promoted_registry.csv
promoted_metrics.csv
summary.json
```

Use `promoted_registry.csv` as the source registry for the subsequent 300-500
candidate one-dimensional campaign. Do not interpret the zero-D ranking as a
validated 100 micrometre toughness prediction.
