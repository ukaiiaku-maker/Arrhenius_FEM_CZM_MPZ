# One-dimensional V2 domain of usefulness

The unconditional screening domain is 300–1200 K and 0–100 µm projected extension. Both providers complete all 40 final class/temperature cases in that domain.

- PF loading-rate domain: 0.01×–100× the canonical opening rate.
- FEM/CZM loading-rate domain: 0.01×–1×; 100× is fail-closed.
- Mechanics/tensor maps: 0–1000 µm extension and 1–1000 µm tip radius.
- Long extension: Peak, weak-T, and ceramic-like complete 1000 µm in both providers. DBTT is limited to 825 µm PF and 448 µm FEM/CZM under the stated bounds.
- Parameter domain: the recorded ±25% local perturbations and the explicit finalist Sobol/interpolation hulls.

Outside the domain the code must invoke an exact deterministic oracle, enrich the map, or fail closed. It never clips extension, radius, front width, or tensor drive. Target completion is right-censoring, not demonstrated arrest. Full limits and observed state ranges are in `oneD_v2_domain_of_usefulness.json`.
