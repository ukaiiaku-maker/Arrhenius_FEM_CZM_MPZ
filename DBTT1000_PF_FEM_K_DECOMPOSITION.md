# DBTT / 1000 K PF-vs-FEM/CZM K-decomposition diagnostic

## Purpose

This branch isolates whether the DBTT/1000 K R-curve discrepancy is primarily
(1) a change in the remote load needed to propagate the crack, (2) a change in
the mapping from remote load to the live local J-equivalent K, or (3) a change
in the tip kinetic state through tip radius/source area/site multiplicity.

The diagnostic is deliberately non-tuning. No cleavage barrier, emission
barrier, event length, RNG, mesh-quality gate, loading rate, tip MPZ law, or
bulk constitutive law is changed.

Branch:

```text
agent/dbtt1000-pf-fem-k-decomposition
```

Base:

```text
claude/v10.0.5.18.4.0-full-physics-advection-recovery
```

## Quantities

The postprocessor `scripts/analyze_dbtt1000_k_decomposition.py` reports:

### Live local K

```text
KJ_local = KJ recorded by the solver from its live J evaluation
```

This includes the solver's current crack representation and mechanical state.

### Remote-load diagnostic using the common fresh reference

```text
K_app_initial_geometry = F_actual * (KJ/F)_fresh
```

The two archived DBTT histories used for the current audit share the same fresh
elastic first row, so this calibration is common to the two calculations. This
quantity is intentionally a **load-only diagnostic**: it asks how the required
remote force changes while holding the microscopic fresh-tip load-to-K mapping
fixed. It is not the final experimental-equivalent K once the projected crack
length has changed substantially.

### Projected-crack experimental-style diagnostic

```text
K_app_SENT = K_app_initial_geometry * g(a_projected)/g(a0)
```

with a conventional single-edge-tension finite-width polynomial. This is only
a secondary diagnostic; it is not asserted to remain ASTM-valid at the largest
a/W. The preferred final apparent-K metric is an independently generated
elastic reference-FEM calibration for the actual specimen and boundary
conditions:

```text
K_app_reference = F_actual * [K_reference/F](a_projected)
```

The postprocessor already accepts such a calibration through
`--reference-calibration`.

### Mapping contribution

```text
DeltaK_map = K_app - KJ_local
M_map = KJ_local / K_app
```

The output also carries the recorded F/U and U/F ratios, tip stress, backstress,
MPZ shielding/counts where present, and a clearly labelled kinematic radius
proxy derived from KJ and sigma_tip. F/U is treated as a recorded global
load/displacement diagnostic until its exact interpretation is revalidated
against the solver boundary-condition bookkeeping; it is not automatically
called an experimental tangent stiffness.

## Dynamic tip site audit

The production emission engine already records, transactionally, for each
accepted emission event:

```text
multiplicity_at_event
rate_per_site_at_event_s
aggregate_hazard_at_event_s
tip_radius_before_event_m
tip_radius_after_event_m
front_width_at_event_m
```

`scripts/run_with_tip_site_audit.py` serializes this existing accepted history
after the normal driver returns. It does not alter site multiplicity, hazard,
RNG, radius evolution, source geometry, crack path, or constitutive state.

This is preferable to introducing new state into the production kinetics.

## Current archived-trajectory comparison

The decomposition was exercised on the two supplied DBTT/1000 K histories that
share the exact same fresh elastic first row:

- FEM/CZM history: final projected extension ~1000.10 um.
- strong-R-curve PF-reference candidate: final projected extension ~1000.01 um.

The strong-R-curve history must remain labelled a reference candidate until its
originating directory manifest is attached to the CSV itself; its behavior and
fresh first row are consistent with the PF comparison data supplied for this
audit.

The event-row selector uses either an accepted `n_fire > 0` event or a positive
committed crack-extension increment. With that event-resolved selection:

```text
first fracture:
    PF/FEM remote-force ratio          ~0.48833
    PF - FEM live KJ                  ~+0.2759 MPa sqrt(m)

last accepted event near 1 mm:
    PF/FEM remote-force ratio          ~0.9705
    PF - FEM live KJ                  ~+25.48 MPa sqrt(m)
```

The FEM/CZM event-state live KJ is ~59.30 MPa sqrt(m) at first fracture and
~59.25 MPa sqrt(m) at the last event near 1 mm. The recorded F/U ratio changes
by less than ~1% between those two events despite ~1 mm projected crack advance.

The strong-R-curve reference candidate begins near the same first-fracture KJ
but reaches ~84.73 MPa sqrt(m) at the last event near 1 mm. Its recorded F/U
history changes strongly.

Therefore the disagreement is not adequately described as a difference in
material resistance alone. The two crack representations assign very different
relationships among projected crack advance, global load/displacement response,
and local J-equivalent K.

A second important observation is that the mappings have already diverged at
first fracture: nearly equal local KJ is reached at very different remote
force. The final audit must therefore distinguish a pre-initiation
constitutive/J-map difference from the post-initiation sharp-wake versus
adaptive-CZM geometry mapping.

## Fresh-run fail-closed requirement

The authoritative theta=0 DBTT comparison requires the exact signed kernel
family whose production SHA-256 is:

```text
d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3
```

Local provenance identifies the source as the v10.2.28 kernel-cache family
with configuration fingerprint:

```text
1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5
```

The exact `family.json` bytes are not tracked in the current remote FEM repo or
the PF GitHub branch. A different Library `family.json` candidate was checked
and rejected because its SHA-256 did not match. The diagnostic therefore fails
closed rather than substituting another kernel.

Once the exact compact artifact is made available to the branch/runtime, the
fresh instrumented FEM/CZM rerun can use `run_with_tip_site_audit.py` and the
matched PF v10.4.1 theta=0 DBTT/1000 K case can be run without changing the
physics.

## CI

`.github/workflows/dbtt1000-k-decomposition.yml` runs the decomposition unit
tests and syntax/import checks and records the production-input preflight.

Workflow run `32057594901` completed successfully for commit
`17e02496ee88c8d839489ee1c4a1457851260f81`.

## Next exact calculation after the kernel is restored

For both PF and FEM/CZM, record at accepted fracture events:

```text
projected crack extension
Uapp
Ftop
live J and KJ
reference-FEM K_app
recorded F/U and U/F
sigma_tip
sigma_back
K_shield
r_tip
front width
source area
site multiplicity
per-site emission rate
aggregate emission hazard
```

Then compare separately:

1. reference-mechanics apparent R-curve, `K_app(a)`;
2. live local R-curve, `KJ(a)`;
3. `DeltaK_map(a) = K_app(a)-KJ(a)`;
4. tip-geometry/site evolution and its saturation length;
5. global load/displacement evolution;
6. changes that occur before first cleavage versus changes caused by crack/wake growth.
