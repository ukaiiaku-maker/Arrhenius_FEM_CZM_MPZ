# MPZ v9.12 reduced work and line-energy bookkeeping

## Purpose

The corrected v9.12 1-D state may show an intermediate-temperature peak in
`Delta_K_micro(T)` followed by a declining high-temperature shielding tail.  A
declining stationary-crack resistance increment does not imply that the system
absorbs less energy through plastic flow.  This diagnostic extension therefore
tracks reduced work and stored-line-energy quantities in parallel with
`Delta_K_micro`.

The diagnostics do **not** feed back into emission, Peierls transport, encounter
storage, Taylor release, recovery, annihilation, shielding, or cleavage.

## Reduced Orowan shear rate

For slip system `alpha`, the two signed dislocation populations move in opposite
spatial directions.  Their Burgers-sign/velocity products add in plastic shear:

```text
gamma_dot[alpha,i]
  = b * (rho_m[alpha,+,i] + rho_m[alpha,-,i]) * v_P[alpha,i].
```

The active-strip cell volume per unit crack-front thickness is

```text
dV/B = dx * active_strip_width.
```

All integrated work quantities therefore have units of `J/m` of crack-front
thickness.

## Stress and work decomposition

The bookkeeping separates four resolved stress contributions:

```text
tau_applied    from K_applied before shielding
tau_shield     = tau(K_eff) - tau(K_applied)
tau_GND        local signed GND stress
tau_eff        = tau_applied + tau_shield + tau_GND
```

where

```text
K_eff = max(K_applied - K_shield, 0).
```

The integrated signed work terms are

```text
W_external = integral tau_applied * gamma_dot dV dt
W_shield   = integral tau_shield  * gamma_dot dV dt
W_GND      = integral tau_GND     * gamma_dot dV dt
W_eff      = W_external + W_shield + W_GND.
```

`W_external` is the reduced applied resolved-work proxy.  `W_shield` and `W_GND`
are signed microstructural exchange terms.  They are not required to be
positive.  The nonnegative diagnostic dissipation is

```text
D_eff = integral max(tau_eff * gamma_dot, 0) dV dt.
```

Power is evaluated at the endpoints of each coupled half-step and integrated by
trapezoidal quadrature.  This bookkeeping does not change the state update.

## Per-crack-area normalization

At each accepted extension checkpoint,

```text
W_external / Delta_a
D_eff      / Delta_a
```

are reported in `J/m^2`.  They are reduced 1-D work-per-crack-area proxies under
the prescribed protocol.  They are useful for comparing temperature trends,
but they are not automatically a path-independent FEM `J`.

## Dislocation line-energy proxy

Mobile and retained line energies use the common isotropic logarithmic estimate

```text
E_line = G b^2 / [4 pi (1-nu)] * ln(R/r_c),
R       = 1 / [2 sqrt(rho_f)],
r_c     = core_regularization_b * b.
```

The line energy is integrated over the active strip per unit crack-front
thickness.  The reported fields are

```text
mobile_line_energy_J_per_m
retained_line_energy_J_per_m
total_line_energy_J_per_m.
```

This is a line-tension proxy.  It is not the complete elastic energy of the
signed GND field and does not replace a mechanically evaluated 2-D energy
balance.

## Output fields

Each `TxxxK.json` includes checkpoint histories for

```text
external_plastic_work_J_per_m
nonlocal_shielding_work_J_per_m
internal_stress_work_J_per_m
effective_plastic_work_J_per_m
effective_plastic_dissipation_J_per_m
external_plastic_work_per_crack_area_J_m2
effective_plastic_dissipation_per_crack_area_J_m2
mobile_line_energy_J_per_m
retained_line_energy_J_per_m
total_line_energy_J_per_m.
```

The integration metadata records

```text
energy_bookkeeping = reduced_1d_orowan_power_and_log_line_energy_v1
energy_bookkeeping_feedback_active = false
energy_units_scope = per_unit_crack_front_thickness.
```

## Interpretation limits

These diagnostics are intended to test whether high-temperature
`Delta_K_micro` softening coexists with increasing plastic dissipation.  They are
not, without the full mechanical solution:

- an ASTM-valid `J_IC`;
- a full `J-R` curve;
- a path-independent contour `J`;
- total specimen load-displacement work;
- Charpy impact energy.

The 2-D FEM/CZM stage should eventually compute external work, recoverable
elastic energy, cohesive/fracture work, plastic dissipation, and stored-state
changes from the complete mechanical fields.

## Corrected balanced-96 rerun

Run the four-candidate smoke first:

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
MODE=smoke \
bash scripts/run_mpz_v9_12_corrected_energy_balanced96.sh
```

Then run the complete balanced set:

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
MODE=full \
bash scripts/run_mpz_v9_12_corrected_energy_balanced96.sh
```

The production calculation uses

```text
MPZ_V912_COUPLED_OPERATOR_SUBSTEPS = 2
MPZ_V912_MAX_FEEDBACK_SUBSTEP_S    = 0.025 s
```

and retains the per-temperature JSON files required for work-versus-extension
analysis.
