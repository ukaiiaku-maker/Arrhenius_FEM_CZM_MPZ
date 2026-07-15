# MPZ v9.10.3 target-aware DBTT search

## Purpose

The v9.10.2 independent-shape search found many analytically accepted DBTT
candidates, but almost all exploited an essentially zero 300 K initiation
shelf.  The previous DBTT objective only checked whether the high-temperature
response rose sufficiently above the low-temperature response.  It did not rank
or constrain the absolute low-temperature toughness.

Version 9.10.3 changes only the DBTT response objective and acceptance gate.  It
retains the full 29-dimensional v9.10.2 parameter domain, the independent
cleavage/emission/Peierls/Taylor EXP-floor shapes, the strict Peierls-below-
Taylor barrier ordering, and the unified Peierls-transport/Taylor-retention MPZ
physics.

## Target-aware objective

The optimizer now uses the temperature-resolved DBTT rows in
`mpz_three_class_design_targets.csv` at 300, 700, 900, and 1200 K.  It scores:

- initiation toughness;
- plateau toughness;
- early R-curve development;
- late/plateau slope;
- the temperature-dependent `delta_KR` windows.

A smooth response-only guard removes the zero-shelf shortcut:

- 300 K initiation shelf: 8--25 MPa sqrt(m);
- 300 K plateau ceiling: 28 MPa sqrt(m);
- 300 K R-curve increment: at most 3 MPa sqrt(m);
- median 900/1200 K initiation: at least 25 MPa sqrt(m);
- median 900/1200 K plateau: at least 35 MPa sqrt(m);
- high-minus-low plateau rise: at least 15 MPa sqrt(m);
- median 900/1200 K R-curve increment: at least 5 MPa sqrt(m).

These are broad promotion gates around the existing design targets.  They do not
require a specific shielding fraction, retained population, blunting mechanism,
source-depletion mechanism, Peierls rate, or Taylor rate.

## Search initialization

Every restart begins from a fresh full Sobol population over all 29 parameters.
No v9.10.2 shortlist or prior DBTT candidate is used as an initial population.
The default production settings use approximately 35,000 global evaluations:

```text
29 parameters
DE_POPSIZE=6
Sobol population rounded to 256
DE_MAXITER=45
RESTARTS=3
```

## Workflow

Run the target-aware global search:

```bash
bash run_mpz_v9_10_3_dbtt_targeted_global_search.sh
```

Then screen the promoted candidates in the unchanged v9.10.2 independent-shape
spatial MPZ:

```bash
bash run_mpz_v9_10_3_dbtt_spatial_promotion.sh
```

The default spatial screen is 50 micrometres.  Candidates that preserve a
finite low-temperature shelf and a developed high-temperature branch should be
rerun to 500 micrometres before 2-D FEM/CZM validation.

## Version boundary

The ceramic and weak-temperature/FCC-like candidates already passed 500-
micrometre reduced spatial validation.  Version 9.10.3 is a DBTT-only calibration
gate.  It does not alter their selected parameters or the active installed
v0.9.6 production state.
