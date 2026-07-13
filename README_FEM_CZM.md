# Arrhenius FEM/CZM package

This package is a full-physics branch of the existing Arrhenius fracture/fatigue solver. It preserves the established hazard, fatigue, plasticity, process-zone, anisotropy, and branching machinery while introducing a replaceable FEM crack-geometry backend.

Start with [FEM_CZM_ARCHITECTURE.md](FEM_CZM_ARCHITECTURE.md).

## Verification commands

```bash
python -m compileall -q arrhenius_fracture
pytest -q tests/test_czm_backend.py
```

A small legacy backend run and a CZM backend run can be launched through `arrhenius_fracture.sharp_front`; the existing campaign drivers also accept `--crack-backend`.


## Oriented DBTT temperature sweep

The oriented DBTT driver now defaults to the angle-faithful `adaptive_czm` backend:

```bash
THETA=30 \
OUTROOT=runs/dbtt_czm_theta30 \
bash run_dbtt_czm_orientation_temperature_test.sh
```

Override `TEMPS`, `MAX_JOBS`, and `FORCE` as needed. The driver preserves the existing canonical DBTT barrier parameterization.
