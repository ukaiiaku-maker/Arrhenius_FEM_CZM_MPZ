# Multifront v4: clustered J decomposition for unresolved branches

This patch implements a two-scale crack-driving treatment for the sharp-front
multifront model.

## Motivation

The previous `localJ` version let newborn branch tips request an independent
small-radius J integral too early.  That can starve the daughter because a
sub-contour branch stub is still inside the parent process zone and is not yet
an independent crack-tip singularity.  It also accidentally made the primary
crack use the small local contour, changing the base driving response and in
some cases preventing initiation/forward growth.

## New decomposition

The code now distinguishes two J scales:

- `--rJ-cluster`: legacy/domain-integral length `ell` for the parent or
  unresolved branch cluster.  Actual outer radius is approximately
  `8*rJ_cluster`.  Default preserves the earlier multifront driving scale:
  `max(10*L_pz, 5e-6)`.
- `--rJ-outer`: actual outer radius for independent local daughter-tip J
  integrals after branch handoff.

With `--j-decomposition cluster` (default):

1. The original parent/root crack uses the outer cluster J.
2. A parent with unresolved daughter stubs also uses the outer cluster J.
3. Newborn/unresolved daughters co-grow under their parent cluster event budget.
4. Once a daughter exceeds `--branch-resolve-length`, it is promoted to an
   independent front and uses its own small local J contour.

The cluster contour excludes the parent and unresolved daughter segments from
line-of-sight blocking so that the contour encloses the whole unresolved kink
cluster rather than being split by its own branch stub.

## Energy allocation

The parent/daughter energy constraint remains local to a bifurcation.  Parent
and unresolved daughter split the parent event budget using the hazard weights,
but unrelated active fronts do not divide one global crack-increment budget.

## Diagnostics

`fronts_<T>K.csv` now includes:

- `J_source_code`: 0 = cluster/group J, 1 = independent local J,
  2 = unresolved daughter using parent-cluster drive.
- `cluster_id`: parent/root ID for the current cluster.
- `J_active_elems`: number of active elements in the J domain for resolved tips.

These diagnostics are intended to separate: (i) physical branch arrest,
(ii) under-resolved local J contours, and (iii) unresolved daughter co-growth.

## Recommended hierarchy

Use a scale hierarchy like:

    h_tip < da_phys < rJ_outer < branch_resolve_length < rJ_cluster_outer

Starting production values:

    --tip-h-fine 0.35e-6
    --da-phys 2e-6
    --rJ-outer 12e-6
    --branch-resolve-length 20e-6
    --rJ-cluster 10e-6       # outer radius ~= 80 um

No velocity cap or maximum-advance cap is introduced.
