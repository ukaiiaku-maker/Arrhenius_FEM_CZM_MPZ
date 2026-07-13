"""
Parameter sweep for Arrhenius fracture: map how barrier height H0 and
Taylor hardening alpha control the Kc(T) trend.

Three expected regimes:
  1. DBTT (W-like): moderate H0 ~ 1.7 eV, moderate alpha
     → Kc increases sharply with T as plastic zone grows
  2. Ceramic-like: high H0 >> kBT or low alpha
     → Kc decreases at high T due to global thermal softening
  3. Aluminum-like: low H0 << kBT, high alpha
     → Kc roughly constant, plastic fracture at all T

Usage:
    python -m arrhenius_fracture.sweep                   # default sweep
    python -m arrhenius_fracture.sweep --vary H0         # sweep H0 only
    python -m arrhenius_fracture.sweep --vary alpha      # sweep alpha only
    python -m arrhenius_fracture.sweep --vary both       # 2D sweep
    python -m arrhenius_fracture.sweep --quick            # fewer points
"""

import numpy as np
import os
import time
import argparse
import copy
from typing import Dict, List, Tuple

from .config import SimulationConfig, make_emergent_config, EV_TO_J, KB
from .main import run_simulation
from .diagnostics import history_summary


def selected_Kc_MPa(hist) -> float:
    """Selected apparent toughness used in sweep plots."""
    return float(history_summary(hist).get('KJ_selected_MPa_sqrt_m', 0.0))


def save_sweep_csv(results: Dict[float, Dict[float, float]], param_name: str,
                   outdir: str, filename: str):
    """Save compact sweep table with selected Kc values."""
    import csv
    os.makedirs(outdir, exist_ok=True)
    params = sorted(results.keys())
    temps = sorted(next(iter(results.values())).keys()) if results else []
    path = os.path.join(outdir, filename)
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([param_name] + [f'{T:.0f}K' for T in temps])
        for p in params:
            writer.writerow([p] + [results[p].get(T, 0.0) for T in temps])
    print(f"Saved sweep CSV to {path}")


def sweep_H0(
    H0_eV_list: List[float],
    T_list: List[float] = None,
    n_steps: int = None,
    base_cfg: SimulationConfig = None,
) -> Dict[float, Dict[float, float]]:
    """
    Sweep barrier height H0 at fixed v*.

    H0 controls the DBTT temperature and softening onset:
      - High H0: softening only at very high T → DBTT behavior
      - Low H0: softening at all T → ductile/constant Kc
      - Very high H0: no plasticity → brittle

    Returns dict: {H0_eV: {T: KJ_peak_MPa}}
    """
    if T_list is None:
        T_list = [300, 500, 700, 900, 1100]

    all_results = {}

    for H0_eV in H0_eV_list:
        cfg = copy.deepcopy(base_cfg) if base_cfg is not None else make_emergent_config()
        cfg.T_list = T_list
        if n_steps is not None:
            cfg.loading.n_steps = n_steps

        # Set barrier: H0 is value at sig0, H(0) = H0/(1-chi)
        cfg.plasticity_barrier.H0_J = H0_eV * EV_TO_J
        # Unique run directory: prevents parameter sweeps from overwriting
        # the history/snapshot files for previous parameter values.
        cfg.output_dir = os.path.join(cfg.output_dir, f'H0_{H0_eV:.3f}eV')

        print(f"\n{'#'*70}")
        print(f"  SWEEP: H0 = {H0_eV:.3f} eV  →  H(0) = {H0_eV/(1-cfg.plasticity_barrier.chiH):.3f} eV")
        print(f"{'#'*70}")

        results = run_simulation(cfg)

        kj = {}
        for T, hist in results.items():
            kj[T] = selected_Kc_MPa(hist)

        all_results[H0_eV] = kj

    return all_results


def sweep_v_star(
    v_star_b3_list: List[float],
    T_list: List[float] = None,
    H0_eV: float = 0.51,
    n_steps: int = None,
    base_cfg: SimulationConfig = None,
) -> Dict[float, Dict[float, float]]:
    """
    Sweep activation volume v* at fixed H0.

    v* controls the stress sensitivity and plastic zone size:
      - Small v*: high flow stress → small plastic zone
      - Large v*: low flow stress → large plastic zone

    Returns dict: {v_star_b3: {T: KJ_peak_MPa}}
    """
    if T_list is None:
        T_list = [300, 500, 700, 900, 1100]

    all_results = {}

    for v_b3 in v_star_b3_list:
        cfg = copy.deepcopy(base_cfg) if base_cfg is not None else make_emergent_config()
        cfg.T_list = T_list
        if n_steps is not None:
            cfg.loading.n_steps = n_steps

        cfg.plasticity_barrier.H0_J = H0_eV * EV_TO_J
        # v0_c is value at sig0, v(0) = v0/(1-psi)
        cfg.plasticity_barrier.v0_c = v_b3 * (cfg.material.b ** 3)
        cfg.output_dir = os.path.join(cfg.output_dir, f'vstar_{v_b3:.3f}b3')

        print(f"\n{'#'*70}")
        v0_actual = v_b3 / (1 - cfg.plasticity_barrier.psiV)
        print(f"  SWEEP: v* = {v_b3:.2f} b³ at σ₀  →  v*(0) = {v0_actual:.2f} b³")
        print(f"{'#'*70}")

        results = run_simulation(cfg)

        kj = {}
        for T, hist in results.items():
            kj[T] = selected_Kc_MPa(hist)

        all_results[v_b3] = kj

    return all_results


def plot_sweep(results: Dict[float, Dict[float, float]],
               param_name: str, param_unit: str,
               outdir: str = None):
    """Plot Kc(T) curves for each parameter value."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available")
        return

    if outdir:
        os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 7))
    cmap = plt.cm.viridis
    params = sorted(results.keys())
    n = len(params)

    for i, p in enumerate(params):
        kj = results[p]
        temps = sorted(kj.keys())
        vals = [kj[T] for T in temps]
        color = cmap(i / max(n-1, 1))
        ax.plot(temps, vals, 'o-', color=color, linewidth=2, markersize=7,
                label=f'{param_name} = {p:.3f} {param_unit}')

    ax.set_xlabel('Temperature (K)', fontsize=13)
    ax.set_ylabel('Apparent $K_c$ (MPa·√m)', fontsize=13)
    ax.set_title(f'Fracture toughness vs temperature\n(varying {param_name})', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)

    if outdir:
        path = os.path.join(outdir, f'sweep_{param_name.replace(" ", "_")}.png')
        fig.savefig(path, dpi=200, bbox_inches='tight')
        print(f"Saved sweep plot to {path}")
    plt.close(fig)


def print_sweep_table(results: Dict[float, Dict[float, float]],
                      param_name: str):
    """Print results as a table."""
    params = sorted(results.keys())
    temps = sorted(next(iter(results.values())).keys())

    header = f"{'':>12}" + "".join(f"  {T:7.0f}K" for T in temps)
    print(f"\n  Apparent Kc (MPa·√m) vs Temperature")
    print(f"  {param_name:>10}" + header)
    print("  " + "-" * (12 + 9 * len(temps)))

    for p in params:
        row = f"  {p:12.4f}"
        for T in temps:
            row += f"  {results[p].get(T, 0):8.3f}"
        print(row)


def main():
    parser = argparse.ArgumentParser(description='Parameter sweep')
    parser.add_argument('--vary', choices=['H0', 'v_star', 'both'],
                        default='H0', help='Which parameter to sweep')
    parser.add_argument('--quick', action='store_true',
                        help='Fewer sweep points')
    parser.add_argument('--steps', type=int, default=None,
                        help='Load steps per run')
    parser.add_argument('--temperatures', nargs='+', type=float, default=None)
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--memory-mode', choices=['off', 'weak_stage1', 'stage1'],
                        default='stage1', help='Reduced crack-tip memory ablation')
    parser.add_argument('--no-save-run-data', action='store_true',
                        help='Do not save individual run histories/snapshot images')
    parser.add_argument('--no-field-pngs', action='store_true',
                        help='Save .npz fields but skip field snapshot PNGs')
    parser.add_argument('--save-every', type=int, default=5,
                        help='Field snapshot cadence in load steps; final state is always saved')
    parser.add_argument('--enable-kinetic-damage-drive', action='store_true',
                        help='Enable extra Model-A kinetic drive in addition to variational AT2')
    parser.add_argument('--pf-damage-cap', type=float, default=None,
                        help='Maximum phase-field damage increment per stagger iteration')
    parser.add_argument('--no-auto-stop', action='store_true',
                        help='Disable force-drop auto-stop so full step history is saved')
    parser.add_argument('--rho-cap', type=float, default=None,
                        help='Maximum dislocation density [m^-2]')
    parser.add_argument('--dot-ep-max', type=float, default=None,
                        help='Maximum plastic strain rate [1/s]')
    parser.add_argument('--max-plastic-strain-increment', type=float, default=None,
                        help='Maximum equivalent plastic strain increment per stagger')
    parser.add_argument('--max-rho-relative-increment', type=float, default=None,
                        help='Maximum fractional rho change per stagger')
    parser.add_argument('--disable-plasticity', action='store_true',
                        help='Diagnostic ablation: mechanics/fracture only, no plastic update')
    parser.add_argument('--freeze-rho', action='store_true',
                        help='Diagnostic ablation: allow plastic strain but hold rho fixed')
    parser.add_argument('--disable-wp-gc-coupling', action='store_true',
                        help='Diagnostic ablation: set emergent Gc_local=Gc, no Wp toughening')
    parser.add_argument('--stop-on-invalid', action='store_true',
                        help='Fail fast when rho, Wp/Wext, Gc_local, KJ, or d_frac becomes invalid')
    parser.add_argument('--invalid-wp-wext-pct', type=float, default=None,
                        help='Stop-on-invalid threshold for 100*Wp/Wext')
    parser.add_argument('--invalid-min-step', type=int, default=None,
                        help='Do not apply stop-on-invalid before this step (default 3)')
    parser.add_argument('--invalid-wext-min', type=float, default=None,
                        help='Minimum Wext before Wp/Wext invalid check is active')

    args = parser.parse_args()

    outdir = os.path.join(os.getcwd(), 'results_sweep')
    os.makedirs(outdir, exist_ok=True)

    T_list = args.temperatures or [300, 500, 700, 900, 1100]

    base = make_emergent_config()
    base.tip_memory.mode = args.memory_mode
    base.tip_memory.enabled = args.memory_mode != 'off'
    if args.steps:
        base.loading.n_steps = args.steps
    if args.enable_kinetic_damage_drive:
        base.phase_field.use_kinetic_damage_drive = True
    if args.pf_damage_cap is not None:
        base.phase_field.max_damage_increment_per_stagger = args.pf_damage_cap
    if args.rho_cap is not None:
        base.dislocations.rho_cap = args.rho_cap
    if args.dot_ep_max is not None:
        base.dislocations.dot_ep_max = args.dot_ep_max
    if args.max_plastic_strain_increment is not None:
        base.dislocations.max_plastic_strain_increment = args.max_plastic_strain_increment
    if args.max_rho_relative_increment is not None:
        base.dislocations.max_rho_relative_increment = args.max_rho_relative_increment
    if args.disable_plasticity:
        base.dislocations.enable_plasticity = False
    if args.freeze_rho:
        base.dislocations.freeze_rho = True
    if args.disable_wp_gc_coupling:
        base.phase_field.plastic_work_to_Gc_efficiency = 0.0
        base.phase_field.Gc_local_cap_factor = 1.0
    if args.stop_on_invalid:
        base.stop_on_invalid = True
    if args.invalid_wp_wext_pct is not None:
        base.invalid_wp_wext_pct = args.invalid_wp_wext_pct
    if args.invalid_min_step is not None:
        base.invalid_min_step = args.invalid_min_step
    if args.invalid_wext_min is not None:
        base.invalid_Wext_min = args.invalid_wext_min
    if args.no_auto_stop:
        base.auto_stop.enabled = False

    # Save individual run data by default. This is intentionally enabled for
    # sweeps so every parameter/T point has raw fields, scalar CSV diagnostics,
    # JSON summary, and a field snapshot PNG.
    base.output_dir = os.path.join(outdir, 'runs')
    base.save_to_disk = not args.no_save_run_data
    base.diagnostics.save_fields = not args.no_save_run_data
    base.diagnostics.save_every = max(args.save_every, 1)
    base.diagnostics.save_field_pngs = (not args.no_save_run_data) and (not args.no_field_pngs)
    # Scalar time-series diagnostic plots are useful but numerous; keep them
    # off in large sweeps while still saving field snapshots and CSV/JSON data.
    base.diagnostics.make_plots = False

    if args.vary in ('H0', 'both'):
        if args.quick:
            H0_list = [0.10, 0.30, 0.51, 0.80]
        else:
            H0_list = [0.06, 0.15, 0.30, 0.51, 0.75, 1.00]

        results = sweep_H0(H0_list, T_list, base_cfg=base)
        print_sweep_table(results, 'H0 (eV)')
        save_sweep_csv(results, 'H0_eV', outdir, 'sweep_H0_table.csv')
        if not args.no_plots:
            plot_sweep(results, 'H0', 'eV', outdir)

    if args.vary in ('v_star', 'both'):
        if args.quick:
            v_list = [2.0, 5.0, 7.5, 15.0]
        else:
            v_list = [1.0, 3.0, 5.0, 7.5, 10.0, 20.0]

        results = sweep_v_star(v_list, T_list, base_cfg=base)
        print_sweep_table(results, 'v* (b³)')
        save_sweep_csv(results, 'v_star_b3', outdir, 'sweep_v_star_table.csv')
        if not args.no_plots:
            plot_sweep(results, 'v_star', 'b³', outdir)

    print(f"\nSweep complete. Results in {outdir}/")


if __name__ == '__main__':
    main()
