# v9.12 full FEM/CZM and material-transfer audit

## Repository/branch scope

The criticism that “there is no FEM/CZM in the repository” does not describe this
branch.  `v9.12-full-field-material-rcurve` descends from the full
`v9.11-full-2d-three-class-integration` branch and retains:

- the elastic–plastic finite-element solve in `arrhenius_fracture/fem.py`;
- the domain/interaction-integral implementation in
  `arrhenius_fracture/j_integral.py`;
- adaptive cohesive/topological insertion in
  `arrhenius_fracture/crack_backend.py`;
- anisotropic mixed-mode, branching, coalescence, remeshing, fatigue/restart and
  field-history logic in `arrhenius_fracture/sharp_front.py`.

The separate `v1-v9102-parity-audit` / `ReusableSourceProcessZone` review quoted
in the project discussion concerns a different reduced constitutive repository or
branch.  Its findings must be mapped to an active call path before they are used to
modify this solver.

## Findings from the v9.11 700 K tip-only smoke

The three material classes used distinct manifests and produced distinct initiation
loads, but their post-initiation normalized topology-event histories were nearly
identical.  The common shape was not evidence that the barriers were identical.
It resulted from a combination of:

1. the same mesh/path and fixed physical advance increment;
2. the same stochastic threshold sequence being reused for all three classes;
3. an imposed class-independent `event_reload` increment after every topology
   event;
4. very small retained-line shielding compared with the applied K range;
5. discrete stochastic source emission changing the single-realization plastic
   response relative to the deterministic mean used during parameter transfer.

A serialized cohesive-edge history is therefore not accepted as a material
resistance curve solely because it contains many points.

## v9.12 protocol

The v9.12 material-transfer gate uses:

- full 2-D FEM/CZM mechanics;
- tip-source-only plasticity;
- stochastic cleavage first passage;
- deterministic expected finite-site emission by default, preserving the
  parameterized mean response;
- independent reproducible threshold streams across material classes by default;
- physical raw fixed-displacement propagation by default;
- cascade-aware classification of same-load topology events;
- mandatory full-field images at several accepted states.

`event_reload`, common random numbers, and discrete stochastic emission remain
available as explicit diagnostics, not production defaults.

## Required field output

Each completed case must write `field_snapshots_<T>K.png`.  The composite image
contains rows for:

1. damage with explicit crack-path/front overlay;
2. log10 dislocation density;
3. maximum principal FEM stress;
4. equivalent plastic strain.

The material-transfer audit fails its publication gate when these images are
missing or when pairwise normalized event shapes remain geometry dominated.

## Interpretation rule

If raw fixed-displacement propagation produces one large same-load cascade, the
result is reported as an unstable brittle jump, not plotted as a smooth R-curve.
A material R-curve requires multiple mechanically independent reload events whose
shape is not merely a scaled copy of the common finite-specimen compliance/path
response.
