# v10.0.5.11 — same-production-mesh J versus fixed-grip energy audit

v10.0.5.10 enlarged the physical refined disk from 100 to 330 micrometers. The
production J contour became stable over 180–300 micrometers, but the 240
micrometer `J/sigma^2` remained 10.54% below the v10.0.5.8 fixed-grip reference.
The 300 micrometer contour was already within 10%.

That remaining comparison mixes two mesh topologies: v10.0.5.8 used a symmetric
crack-perturbation audit mesh, while production uses a forward Mode-I corridor.
v10.0.5.11 therefore evaluates the thermodynamic identity on the exact production
mesh and boundary state.

## Audit definition

The probe retains the full v10.0.5.10 production path and, after the initial
maximum-load equilibrium, resolves three elastic starter-notch states on the same
mesh at the first three explicit corridor stations:

```text
a0, a0+h, a0+2h
```

It computes the second-order forward fixed-grip release rate

```text
G = [3 U(a0) - 4 U(a0+h) + U(a0+2h)] / (2 h)
```

and compares it with the median full-energy J over the accepted 180, 240, and 300
micrometer contour plateau.

The v10.0.5.8 cross-mesh reference remains in the output as a diagnostic only. No
tolerance is relaxed.

## Pass status

Only

```text
production_J_same_mesh_energy_parity_passed
```

passes. The audit also fails closed for inadequate physical support, contour
instability, unconverged forward energy derivative, nonquadratic G scaling, or a
same-mesh J/G mismatch.

## Validation

```bash
python -m py_compile \
  arrhenius_fracture/production_j_same_mesh_energy_v100511.py \
  arrhenius_fracture/mode_i_first_passage_v10_0_5_11_same_mesh_probe.py \
  run_v10_0_5_11_same_mesh_energy_audit.py

pytest -q \
  tests/test_v10056_kj_audit_bracket.py \
  tests/test_v10056_audited_wrappers.py \
  tests/test_v10058_fixed_grip_elastic_audit.py \
  tests/test_v10058_fixed_grip_audit_mesh.py \
  tests/test_v10059_production_j_parity.py \
  tests/test_v10059_v911_probe_contract.py \
  tests/test_v10059_straight_path_recorder.py \
  tests/test_v100510_refinement_support.py \
  tests/test_v100511_same_mesh_energy.py
```

## Run

Use the same command as v10.0.5.10, replacing the runner and output directory:

```bash
REFERENCE=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_8_fixed_grip_elastic_convergence/runs/v10_0_5_8_fixed_grip_elastic_convergence_v2/fixed_grip_elastic_convergence_v10_0_5_8.json

python run_v10_0_5_11_same_mesh_energy_audit.py \
  --reference-json "$REFERENCE" \
  --openings-um "1 2" \
  --contour-outer-um "100 140 180 240 300" \
  --accepted-contour-um "180 240 300" \
  --selected-outer-um 240 \
  --tip-refinement-radius-um 330 \
  --contour-stability-rel-tol 0.10 \
  --nx 36 --ny 72 \
  --tip-h-um 2.5 --tip-ratio 1.15 \
  --cluster-J-outer-um 240 --local-J-outer-um 100 \
  --mpz-length-um 100 --mpz-n-bins 80 \
  --da-phys-um 5 \
  --crack-backend adaptive_czm \
  --crystal-aniso --crystal-theta-deg 45 \
  --out runs/v10_0_5_11_same_mesh_energy_v1
```

Primary output:

```text
production_j_same_mesh_energy_v10_0_5_11.json
```
