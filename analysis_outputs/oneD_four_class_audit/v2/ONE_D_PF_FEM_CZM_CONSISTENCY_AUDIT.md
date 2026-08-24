# One D Pf Fem Czm Consistency Audit

The autonomous v9.13 and generic PF 1-D implementations are not duplicate wrappers and are not intended frozen-state-equal: they differ in barrier/state definitions, external loading, threshold/RNG lifecycle, and post-event renewal. More importantly, the canonical four-class rows cannot be executed by the PF installation's 1-D interface at `0a340f6`. Frozen-state and matched-event numeric equality are therefore unavailable, not failed. Creating an adapter would invent a third model and was rejected. No FEM/CZM-specific 1-D engine exists.
