# DBTT / 1000 K true PF-vs-FEM/CZM common-reference K audit

## Scope and provenance

This report supersedes the earlier loose-CSV PF candidate comparison.

PF input:

```text
T1000K_th0_seed1008666
v913_paper_dbtt01_0202500_persistent_sites
PF v10.2.30 hazard-energy-gated audited model
bulk_plasticity_mode = tip_only
theta = 0 deg
seed = 1008666
steps_1000K.csv sha256 = 233824c0e7728ab05f9cf5aedb495e476c8c42696ecc14258bac4f065c326a5e
```

FEM/CZM input:

```text
DBTT / 1000 K / theta=0 / seed1008666
tip-only moving persistent-site MPZ
adaptive CZM atomic-path corridor
steps_1000K.csv sha256 = a11f28890af57d4ee33fddc0cf16298539efd61d13d455383c02cb36c2438eda
```

No model parameters were fitted or changed for this audit.

## Common apparent-K definition

A common elastic projected-crack mechanics calculation was run on the diagnostic
branch. The reference uses the same 2 mm x 4 mm specimen dimensions, initial
0.5 mm notch, theta=0 cubic elastic constants and symmetric fixed-grip loading,
but no MPZ, plasticity, stochastic kinetics, backstress or evolved microscopic
tip morphology.

Rather than using the local J integral of the reference crack, apparent K is
constructed from global elastic compliance:

```text
G/F^2 = (1/(2B)) dC/da
K/F   = sqrt(E'/(2B) dC/da)
Kapp  = |F_actual| [K/F]_common_reference(a_projected)
```

The fixed-mesh compliance has small staircase discretization from discrete
straight-wake element removal. Polynomial degrees 1-7 were compared by BIC;
a quadratic is selected. The compliance-fit RMSE is about 0.93% of the total
compliance change and the fitted derivative is positive throughout the
calibration interval.

This `Kapp` is an experimental-style *common numerical reference* for the model
specimen, not an ASTM specimen-validity claim.

## Main result

At the first accepted fracture event:

```text
                         FEM/CZM          PF v10.2.30
extension (um)             2.297            2.931
local KJ                  59.301           56.188  MPa sqrt(m)
common-reference Kapp     58.085          154.679  MPa sqrt(m)
Kapp-KJ                   -1.215           98.491  MPa sqrt(m)
```

Near 1 mm projected extension:

```text
                         FEM/CZM          PF v10.2.30
extension (um)          1000.100         1000.825
local KJ                  59.251          109.405  MPa sqrt(m)
common-reference Kapp    751.256          160.066  MPa sqrt(m)
Kapp-KJ                  692.004           50.661  MPa sqrt(m)
```

Thus the two representations invert the apparent/local interpretation:

- PF shows a strong rise in **solver-local KJ** (~56 -> ~109 MPa sqrt(m)), while
  its common remote-load Kapp remains in a much narrower ~155-175 MPa sqrt(m)
  band through most of the trajectory.
- FEM/CZM shows approximately steady **solver-local KJ**, but the remote load
  required for propagation produces a very strong common-reference apparent
  Kapp rise, reaching ~751 MPa sqrt(m) near 1 mm.

The conclusion is not that one of these K values is automatically the
experimental truth. It is that the choice of crack representation and live J
mapping dramatically changes how the same macroscopic load history is projected
onto local KJ. A conventional experiment that infers K from remote load and
macroscopic projected crack length would be much closer in spirit to the common
`Kapp` calculation than to either solver's microscopic live-J definition.

## Tip radius

PF opening-field radius, recovered from the archived KJ and tip opening stress:

```text
first event:  7.57 um
near 1 mm:   10.58 um
```

FEM/CZM radius is exactly reconstructable from archived local slip and the
authoritative law `r = r0 + c_blunt*b*N_local_slip`:

```text
first event:  8.32 um
near 1 mm:    7.11 um
```

The two models therefore have qualitatively different developed tip states:
PF blunts further during long propagation, whereas the FEM/CZM moving-tip state
ends slightly sharper than at first fracture.

## Source area and site multiplicity

PF records the total effective tip-source multiplicity at every step in the
v10.2.30 hazard audit. With two active systems and unit source activity, the
hazard-equivalent source area per system follows directly from
`M_total/(2*rho_source)`.

```text
PF first event:
    source area/system = 2.86 um^2
    total multiplicity = 3000
PF near 1 mm:
    source area/system = 3.30 um^2
    total multiplicity = 3461
```

The archived FEM/CZM step table retains mean backstress but not the two
per-system densities required for an exact front width. Therefore the FEM/CZM
source area and multiplicity are reported as rigorous two-system bounds, not a
point estimate:

```text
FEM/CZM first event:
    source area/system = 2.11 - 2.98 um^2
    total multiplicity = 2209 - 3124
FEM/CZM near 1 mm:
    source area/system = 1.97 - 2.79 um^2
    total multiplicity = 2066 - 2922
```

The PF and FEM/CZM source-population scales overlap substantially. Site count
alone is therefore unlikely to explain the enormous long-growth K-mapping
difference.

## Backstress and signed MPZ K shielding

Backstress:

```text
FEM/CZM: 8.38 -> 7.66 GPa
PF:      7.94 -> 9.55 GPa
```

Signed MPZ shielding:

```text
FEM/CZM: -4.067 -> -0.012 MPa sqrt(m)
PF:      -3.658 -> -0.017 MPa sqrt(m)
```

The direct signed K-shielding term therefore becomes negligible in both models
at long extension. It cannot account for either the ~53 MPa sqrt(m) PF local-KJ
rise or the hundreds of MPa sqrt(m) FEM/CZM Kapp-KJ gap.

Backstress and source multiplicity are kinetic variables in the Arrhenius
emission law; they are not additive K terms. They must not be summed with
`Kapp-KJ` as if they were independent Delta-K contributions.

## Interpretation

The strongest conclusion from this audit is that the apparent disagreement
"PF has an R-curve, FEM/CZM does not" depends critically on which K is plotted.

Using live local J-derived K:

```text
PF       -> strong rising R-curve
FEM/CZM  -> rapidly developed fluctuating near-steady state
```

Using the same macroscopic remote-load/projected-crack calibration:

```text
PF       -> comparatively weak apparent-load evolution
FEM/CZM  -> very strong apparent R-curve
```

This means micron-scale crack representation is not a secondary bookkeeping
issue. It controls the remote-load-to-local-J mapping strongly enough to change
the qualitative R-curve interpretation.

The next exact step, if desired, is an instrumented rerun of the single
FEM/CZM DBTT/1000 K case using `scripts/run_with_tip_site_audit.py` so the
per-system source densities, exact source area, site multiplicity and eventwise
radius/front-width history are serialized directly rather than bounded from the
archived scalar backstress.
