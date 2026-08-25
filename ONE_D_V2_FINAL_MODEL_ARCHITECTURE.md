# One-dimensional V2 final model architecture

## Outcome

The PF-consistent and FEM/CZM-consistent reduced models are fit for 100-µm screening inside the declared domain. They share barrier, hazard, and material-row physics, while mechanics, tensor drive, event length, renewal, translation, reload grouping, and veto behavior remain backend-owned.

The fast FEM/CZM lane is an explicitly versioned event surrogate trained against exact joint-K microfixtures. Its hazard-progress correction, DBTT precursor threshold, six-source-zone retained-state coordinate, and low-density plateau are backend reductions—not material parameters. The rejected dense-regime packet and the former common lifecycle are not used.

Native PF J/KJ and native FEM/CZM J/KJ drive kinetics. Qualified FEM/CZM structural G/K_G is reported separately and is never substituted into the kinetic law. Source-drive maps span 0–1000 µm extension and 1–1000 µm radius and fail closed outside that box.

Production PF/FEM formulas and canonical trajectories were not altered. Six new bounded PF validation trajectories were written only under `/private/tmp/oneD-v2-terminal-pf-transfer-runs`; no 2-D FEM/CZM run was launched.
