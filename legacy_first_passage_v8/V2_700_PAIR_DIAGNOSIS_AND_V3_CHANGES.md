# Diagnosis of the repeated v2 700 MPa pair and v3 corrections

The two uploaded `sn_stateful_pd_v2_700_pair` archives contain identical extracted files. The pasted stdout is consistent with the archived histories.

## What v2 did correctly

At 300 K, 700 MPa, R = 0.1, 1000 Hz, and 1e9 cycles:

| quantity | no shield | shielded |
|---|---:|---:|
| expected cumulative births | 0.256665 | 0.0088267 |
| expected stable defects | 0.256631 | 0.0088253 |
| realized births, seed 1 | 0 | 0 |
| realized stable defects, seed 1 | 0 | 0 |
| maximum bond damage | 0 | 0 |
| connected crack extent | 0 | 0 |

The discrete-site correction therefore worked: no realized stable defect means no cohesive softening.

## Remaining v2 problem

The spatial site measure was assigned through the full PD patch, including the displacement-coupling shell. Large nodal areas and a capped PD/FEM strain ratio concentrated much of the expected first-event weight near the artificial patch boundary rather than the scratch root.

The intact global FEM already resolves the scratch stress concentration. Applying the PD/FEM strain ratio before any cohesive damage therefore double-counted the localization and made low-denominator bonds reach the amplification cap.

## v3 changes

1. Candidate sites are excluded from the coupling shell.
2. The initiation population is localized to a root-centered process zone with a cosine taper.
3. Sites far behind the evolving root are excluded.
4. The initial physical candidate population is fixed and conserved even as the scratch geometry evolves.
5. PD/FEM amplification is exactly one while the bond network is intact and turns on continuously only after cohesive damage develops.
6. The bounded one-delivery-event-per-cycle process is represented as a continuous Poisson intensity, so the memory time remains dimensionally consistent in seconds.
7. Final NPZ output now includes initiation weight, candidate-site measure, hit rate, birth intensity, effective stress, and PD amplification.
8. A six-panel `pd_initiation_diagnostics_final.png` is written for every run.

With the production mesh and seed 1, the default v3 process-zone definition gives approximately 3.75e3 expected candidate sites and 3.68e3 realized sites. All coupling-shell sites are zero. Approximately 26% of the candidate population lies within 120 micrometers of the root, 61% within 180 micrometers, and 100% within 240 micrometers.

The v3 changes intentionally correct spatial bookkeeping and intact-field amplification before any barrier, site-density, or memory calibration is attempted.
