# DBTT / 1000 K PF-vs-FEM/CZM K-decomposition diagnostic

## Provenance correction — 2026-08-17

The first archived comparison on this branch was **not a valid PF-vs-FEM/CZM comparison**. It used two loose `steps_1000K.csv` histories with SHA-256 values:

```text
a11f28890af57d4ee33fddc0cf16298539efd61d13d455383c02cb36c2438eda
b1f7c74da4b7aa73d24fc946e8960fc8e0ae91d0d4fed07619f4ce4a3b6fbc20
```

Those two histories share the same initial mechanics row and the same first stochastic crack increment (`2.2973400956249022e-06 m`). The second was incorrectly treated as a PF reference candidate. It is **not** the authoritative PF DBTT/1000 K production trajectory. Results derived from that pairing must not be interpreted as PF-vs-FEM.

The user-supplied authoritative PF case is:

```text
T1000K_th0_seed1008666
option: v913_paper_dbtt01_0202500_persistent_sites
model entry: arrhenius_fracture.sharp_front_v10_2_30_hazard_energy_gated_audited
bulk_plasticity_mode: tip_only
theta: 0 deg
seed: 1008666
target projected extension: 1000 um
steps_1000K.csv sha256: 233824c0e7728ab05f9cf5aedb495e476c8c42696ecc14258bac4f065c326a5e
```

Its production metadata reports:

```text
Kc_first = 56.1881603166 MPa sqrt(m)
first accepted extension = 2.93110656568 um
final projected extension = 1000.82547666 um
final accepted-event KJ = 109.405272664 MPa sqrt(m)
```

The FEM/CZM history previously used as the FEM curve has:

```text
steps_1000K.csv sha256: a11f28890af57d4ee33fddc0cf16298539efd61d13d455383c02cb36c2438eda
KJ_first = 59.3006676307 MPa sqrt(m)
first accepted extension = 2.29734009562 um
final projected extension = 1000.09984489 um
final accepted-event KJ = 59.2513012872 MPa sqrt(m)
```

Therefore the corrected live-KJ comparison is qualitatively the expected one: the PF trajectory develops a strong rising R-curve from ~56.2 to ~109.4 MPa sqrt(m), whereas the FEM/CZM trajectory remains near ~59.3 MPa sqrt(m) through ~1 mm extension.

## Important consequence for apparent K

The true PF and FEM/CZM runs do **not** share the same fresh force/displacement/K row. Therefore the previous use of one shared fresh `(K/F)` calibration is invalid and must not be used to define experimental-equivalent apparent K.

The correct apparent-K calculation remains:

```text
K_app(a) = F_actual(a) * [K_reference/F](a_projected)
```

where `[K_reference/F](a_projected)` is obtained from an independent, common elastic reference-mechanics calculation using the same macroscopic specimen dimensions and boundary conditions but an ideal projected crack. This is the quantity to compare against each solver's live local `KJ`.

## Dynamic tip-site audit

The instrumentation objective is unchanged. The production emission engine already records accepted-event quantities including:

```text
multiplicity_at_event
rate_per_site_at_event_s
aggregate_hazard_at_event_s
tip_radius_before_event_m
tip_radius_after_event_m
front_width_at_event_m
```

These should be serialized without modifying the stochastic or constitutive evolution, then compared against projected crack extension together with `K_app`, live `KJ`, force, load/displacement response, backstress, and MPZ shielding.

## Interpretation

Do not use the prior loose-CSV PF label or its ~84.7 MPa sqrt(m) terminal curve. The corrected PF terminal value is ~109.4 MPa sqrt(m). The central scientific question is now sharper: why does the PF sharp-wake representation develop ~53 MPa sqrt(m) of live-KJ rise over 1 mm while the adaptive FEM/CZM representation remains essentially flat, and how much of that difference would be visible in a common macroscopic `K_app` measurement rather than only in the solver-local J mapping?
