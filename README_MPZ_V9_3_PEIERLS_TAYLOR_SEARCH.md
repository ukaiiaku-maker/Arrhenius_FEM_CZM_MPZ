# MPZ v9.3: emission-derived Peierls–Taylor plasticity and monotonicity search

## Purpose

Version 9.3 restores one common activated-plasticity architecture across the analytical atlas, moving process zone, cyclic controller, and bulk FEM update. Crack opening and crack-tip emission remain independent EXP-floor free-energy surfaces. Peierls-controlled motion and Taylor-type obstacle escape are scaled descendants of the selected emission surface.

```text
EXP-floor crack-tip emission
    -> scaled EXP-floor Peierls rate
    -> correlated multi-hit scaled EXP-floor Taylor completion
    -> sequential Peierls–Taylor rate
    -> plastic slip and Kocks–Mecking state evolution
```

The legacy additive `sigma_Peierls + sigma_Taylor` implementation remains selectable only through `legacy_additive_flow_stress`.

## High-density Poisson guardrail

The independent Arrhenius–Taylor prefactor can grow too rapidly with forest density and produce a flow-stress maximum followed by a nonphysical decrease. Version 9.3 does not use a fitted total-density cap to hide that turnover. Taylor completion is instead a correlated renewal:

```text
h1       = nu_T exp[-DeltaG_T/(kB T)]
m(rho_f) = cooperative hit order
h_T      = gammainc(m, h1 t_c) / t_c
h_PT     = (1/h_P + 1/h_T)^(-1)
```

A finite `m_cap` denotes the finite number of obstacles in one correlation domain; it is not a cap on total dislocation density. Every sampled closure is rejected when the flow-stress curve turns downward after hit-order saturation. The FEM density ceiling is a remote floating-point overflow guard and must remain inactive.

## Parent and mechanism scales

The effective emission candidate is the parent surface with scale 1.0. Prior centers are:

```text
Peierls energy / emission energy = 0.005
Taylor energy / emission energy  = 0.020
```

Temperature slopes are scaled independently through reported entropy ratios.

## Parameter search

The v9.2 analytical shortlist supplies intrinsic crack-opening/emission candidates. The v9.3 search maps:

- Peierls and Taylor relative energy scales;
- plastic temperature-slope multiplier;
- Taylor correlation density and renewal time;
- multi-hit exponent, scale, and finite correlation-domain obstacle count;
- mobile fraction and mobile-density saturation scale.

The code inverts the common kinetic law over temperature, strain-rate, and forest-density grids. `accepted=True` means that the curve is resolved and does not recover a high-density downturn. `strict_strength_window` additionally requires the 700 K, 1e-5 s^-1, 1e14 m^-2 reference strength to lie in the broad screening window.

## Production command

```bash
nohup env \
TRANSPORT_SAMPLES=1024 \
INTRINSIC_TOP_PER_REGION=5 \
TOP_PER_INTRINSIC=3 \
TEMPERATURES="300 700 900 1200" \
STRAIN_RATES="1e-5 1e-3" \
RHO_MIN=5e12 \
RHO_MAX=1e18 \
RHO_POINTS=65 \
bash run_mpz_v9_3_peierls_taylor_search.sh \
> runs/mpz_v9_3_peierls_taylor_search.nohup.log 2>&1 &
```

The input defaults to `runs/mpz_v9_2_analytic_first_passage_atlas/analytic_first_passage_atlas_shortlist_refined.csv` and falls back to the coarse shortlist when needed.

## Outputs

The default output directory is `runs/mpz_v9_3_peierls_taylor_search`.

- `peierls_taylor_search_all.csv.gz`: full intrinsic × transport map;
- `peierls_taylor_search_accepted.csv`: density-monotone resolved rows;
- `peierls_taylor_search_shortlist.csv`: ranked rows per intrinsic candidate;
- `mpz_v9_3_material_shortlist.csv`: complete rows for the next MPZ stage;
- `peierls_taylor_search_summary.csv`: search counts;
- `peierls_taylor_search_config.json`: exact grid and interpretation;
- `pt_strength_*.png`: representative fixed-rate strength curves.

## Next gate

Passing this search does not establish a DBTT or R-curve. The next stage calculates the developed moving-process-zone state, shielding/blunting response, and evolution length. Only transient virgin-tip runs that approach that developed branch proceed to 2-D FEM/CZM.
