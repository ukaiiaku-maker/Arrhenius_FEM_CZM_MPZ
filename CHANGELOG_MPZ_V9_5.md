# MPZ v9.5 changelog

## Developed-state audit finding

The v9.4 common developed-state search completed all 133 closures, but the score was nearly flat and the best rows suppressed source activity. The reduced MPZ was not reaching the cooperative Taylor-density regime selected by the bulk Peierls–Taylor screen.

The cause was structural rather than a parameter-bound problem:

- forest density was calculated from total retained count divided by the square of the full MPZ length;
- one scalar tip stress was applied to every MPZ bin;
- for 50–200 micrometre zones, tens of near-tip retained lines changed the calculated density by far less than the 5e12 m^-2 floor;
- the Taylor hit-order transition near 1e14 m^-2 therefore never engaged, while the same retained lines still contributed directly to shielding.

## v9.5 correction

- Forest density is now calculated in every bin as retained line count divided by the physical bin area `dx * process_zone_width`.
- The existing blunting length supplies the effective process-zone width, avoiding an unconstrained new multiplier.
- Effective stress decays from the supplied tip stress over the one-dimensional MPZ using an LEFM-like square-root profile.
- Detailed-balance Peierls and correlated Taylor rates are evaluated using the local stress and local forest density arrays.
- Taylor release, trapping, and retained recovery are integrated locally by bin.
- Mobile advection remains conservative and uses the population-weighted/source-zone Peierls velocity.
- A physically dimensioned retained-density initializer was added only for developed-state branch-continuation audits; virgin production runs still begin from zero retained state.
- A continuation audit compares virgin and seeded forest-density branches to determine whether a developed state decays, persists, or attracts the virgin trajectory during crack growth.

The v9.4 developed-state shortlist is retained as a diagnostic of the global-density formulation and must not be interpreted as a calibrated MPZ parameter set.
