# Stateful-PD v3.1 preflight overlay

This overlay contains the same root-localized v3 physics plus a strict source-path
and model-version preflight. It is intended to prevent an older v2 module or
runner from being used silently.

The runner:

- changes to the directory containing the script;
- prepends that directory to `PYTHONPATH`;
- removes stale `__pycache__` directories;
- runs Python with `-B`;
- prints the exact imported driver and PD-module paths;
- aborts unless `MODEL_ID` is the v3.1 root-localized identifier;
- aborts unless all v3 command-line options are present.

A valid v3.1 result must contain:

- model `SN_2D_intact_FEM_stateful_local_peridynamics_v3_1_root_localized`;
- `pd_initiation_diagnostics_final.png` in each case directory;
- initiation-radius/taper/back-extent fields in `summary.json`;
- `boundary`, `initiation_weight`, and `mean_candidate_sites` arrays in
  `pd_state_final.npz`;
- intact-state `pd_amplification_max = 1` until bond damage develops.
