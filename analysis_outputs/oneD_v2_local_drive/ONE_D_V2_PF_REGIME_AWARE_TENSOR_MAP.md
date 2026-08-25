# PF regime-aware tensor map

Status: **HYBRID PROVIDER CONTRACT QUALIFIED; EXTENSION-ONLY SURROGATE UNQUALIFIED**.

`PFProbeRegimeKey` excludes continuous tensor values and includes wake topology,
front basis, probe element support, probe weights, selected-system/sign
signatures, and reliability. The 17 exact source states span
0–1000 µm; all are reliable and retain killed-node IDs, realized event geometry,
probe IDs/weights, raw tensors, reaction, energy, native J/KJ, and profile bounds.
None of the 16 adjacent sparse intervals has an identical regime
key, so no interval is admitted for interpolation. Queries at a transition use
an exact deterministic mechanics query keyed by realized event geometry and are
cached; absence of that callable fails closed. Because the production event
length is a clipped continuous exponential draw, arbitrary trajectories cannot
be represented by a finite extension-only lookup.
