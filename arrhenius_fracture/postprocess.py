"""
Post-process and visualize results from Arrhenius fracture simulations.

Reads the saved .npz history files and produces summary plots:
  1. Kc(T) — apparent toughness vs temperature
  2. Load-displacement curves per temperature
  3. Plastic dissipation fraction vs T
  4. Peak Gc_local vs T (plastic toughening)
  5. Peak dislocation density vs T
  6. Energy balance breakdown

Usage:
    python -m arrhenius_fracture.postprocess                         # default dir
    python -m arrhenius_fracture.postprocess --dir results_arrhenius
    python -m arrhenius_fracture.postprocess --dir results_sweep
"""

import numpy as np
import os
import glob
import argparse

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def load_all_histories(results_dir: str) -> dict:
    """
    Load all history .npz files from a results directory.

    Returns dict: {T_K: npz_data}
    """
    histories = {}
    # Look for files like history_0700K.npz in subdirectories
    patterns = [
        os.path.join(results_dir, '*K', 'history_*K.npz'),
        os.path.join(results_dir, 'history_*K.npz'),
        os.path.join(results_dir, '*.npz'),
    ]

    files = []
    for pat in patterns:
        files.extend(glob.glob(pat))

    if not files:
        print(f"No .npz files found in {results_dir}")
        return histories

    for f in sorted(files):
        try:
            data = np.load(f, allow_pickle=True)
            T = float(data['T'][0])
            histories[T] = data
            print(f"  Loaded T={T:.0f}K from {os.path.basename(f)}")
        except Exception as e:
            print(f"  Warning: couldn't load {f}: {e}")

    return histories


def extract_summary(histories: dict) -> dict:
    """Extract key summary quantities from all temperature runs."""
    summary = {
        'T': [], 'Gc_eff': [], 'Kc_input': [],
        'KJ_peak': [], 'KJ_global_peak': [],
        'Fmax': [], 'Uapp_at_Fmax': [],
        'Wp_final': [], 'Wext_final': [], 'Wp_ratio': [],
        'rho_max': [], 'Gc_max': [],
        'n_steps': [],
    }

    for T in sorted(histories.keys()):
        d = histories[T]
        summary['T'].append(T)
        summary['Gc_eff'].append(float(d['Gc_eff'][0]))
        summary['Kc_input'].append(float(d['Kc_input'][0]))

        KJ_d = d['KJ_domain'] if 'KJ_domain' in d else np.zeros(1)
        KJ_g = d['KJ_global'] if 'KJ_global' in d else np.zeros(1)
        summary['KJ_peak'].append(float(np.max(KJ_d)) / 1e6)
        summary['KJ_global_peak'].append(float(np.max(KJ_g)) / 1e6)

        F = d['Freact'] if 'Freact' in d else np.zeros(1)
        U = d['Uapp'] if 'Uapp' in d else np.zeros(1)
        idx_max = np.argmax(np.abs(F))
        summary['Fmax'].append(float(np.max(np.abs(F))))
        summary['Uapp_at_Fmax'].append(float(U[idx_max]) * 1e6 if len(U) > 0 else 0)

        Wp = d['Wp'] if 'Wp' in d else np.zeros(1)
        Wext = d['Wext'] if 'Wext' in d else np.zeros(1)
        summary['Wp_final'].append(float(Wp[-1]))
        summary['Wext_final'].append(float(Wext[-1]))
        Wext_f = max(float(Wext[-1]), 1e-30)
        summary['Wp_ratio'].append(100 * float(Wp[-1]) / Wext_f)

        rho = d['rho_max'] if 'rho_max' in d else np.zeros(1)
        summary['rho_max'].append(float(np.max(rho)))

        Gc_max = d['Gc_local_max'] if 'Gc_local_max' in d else np.zeros(1)
        summary['Gc_max'].append(float(np.max(Gc_max)))

        summary['n_steps'].append(len(F))

    # Convert to arrays
    for k in summary:
        summary[k] = np.array(summary[k])

    return summary


def plot_summary(histories: dict, summary: dict, outdir: str):
    """Generate a 6-panel summary figure."""
    if not HAS_MPL:
        print("matplotlib not available")
        return

    T = summary['T']

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    # --- Panel 1: Apparent toughness vs T ---
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(T, summary['KJ_peak'], 'bo-', lw=2, ms=7, label='$K_J$ domain (peak)')
    ax1.plot(T, summary['KJ_global_peak'], 'rs--', lw=1.5, ms=6, label='$K_J$ global (peak)')
    ax1.axhline(summary['Kc_input'][0] / 1e6 if summary['Kc_input'][0] > 1 else summary['Kc_input'][0],
                color='k', ls=':', alpha=0.5, label='$K_c^{intrinsic}$')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('$K_c$ (MPa·√m)')
    ax1.set_title('Apparent Toughness')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # --- Panel 2: Load-displacement curves ---
    ax2 = fig.add_subplot(gs[0, 1])
    cmap = plt.cm.coolwarm
    for i, Tk in enumerate(sorted(histories.keys())):
        d = histories[Tk]
        U = d['Uapp'] * 1e6 if 'Uapp' in d else np.zeros(1)
        F = d['Freact'] / 1e3 if 'Freact' in d else np.zeros(1)
        color = cmap(i / max(len(histories) - 1, 1))
        ax2.plot(U, F, '-', color=color, lw=1.5, label=f'{Tk:.0f} K')
    ax2.set_xlabel('Opening (µm)')
    ax2.set_ylabel('Force (kN)')
    ax2.set_title('Load–Displacement')
    ax2.legend(fontsize=8, ncol=2)
    ax2.grid(True, alpha=0.3)

    # --- Panel 3: Peak force and displacement at peak ---
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(T, summary['Fmax'] / 1e3, 'go-', lw=2, ms=7)
    ax3.set_xlabel('Temperature (K)')
    ax3.set_ylabel('Peak force (kN)')
    ax3.set_title('Peak Force vs T')
    ax3.grid(True, alpha=0.3)
    ax3b = ax3.twinx()
    ax3b.plot(T, summary['Uapp_at_Fmax'], 'ms--', lw=1.5, ms=5, alpha=0.7)
    ax3b.set_ylabel('$U$ at $F_{max}$ (µm)', color='m')

    # --- Panel 4: Plastic dissipation ---
    ax4 = fig.add_subplot(gs[1, 0])
    # Cap ratio for display
    ratio = np.clip(summary['Wp_ratio'], 0, 1000)
    ax4.bar(T, ratio, width=60, color='coral', alpha=0.7, edgecolor='brown')
    ax4.set_xlabel('Temperature (K)')
    ax4.set_ylabel('$W_p / W_{ext}$ (%)')
    ax4.set_title('Plastic Dissipation Fraction')
    ax4.grid(True, alpha=0.3, axis='y')

    # --- Panel 5: Max Gc_local (plastic toughening) ---
    ax5 = fig.add_subplot(gs[1, 1])
    Gc_max = summary['Gc_max']
    if np.max(Gc_max) > 0:
        ax5.semilogy(T, np.maximum(Gc_max, 1), 'D-', color='darkgreen', lw=2, ms=7)
        ax5.axhline(summary['Gc_eff'][0], color='k', ls=':', alpha=0.5,
                     label=f'$G_c^{{intrinsic}}$ = {summary["Gc_eff"][0]:.1f}')
        ax5.legend(fontsize=8)
    else:
        ax5.plot(T, Gc_max, 'D-', color='darkgreen', lw=2, ms=7)
    ax5.set_xlabel('Temperature (K)')
    ax5.set_ylabel('Max local $G_c$ (J/m²)')
    ax5.set_title('Plastic Toughening ($G_c + w_p \\cdot \\ell$)')
    ax5.grid(True, alpha=0.3)

    # --- Panel 6: Peak dislocation density ---
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.semilogy(T, summary['rho_max'], '^-', color='purple', lw=2, ms=7)
    ax6.set_xlabel('Temperature (K)')
    ax6.set_ylabel('Peak $\\rho$ (m$^{-2}$)')
    ax6.set_title('Peak Dislocation Density')
    ax6.grid(True, alpha=0.3)

    fig.suptitle('Arrhenius Fracture — Simulation Summary', fontsize=14, y=0.98)

    path = os.path.join(outdir, 'summary_6panel.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f"\nSaved 6-panel summary to {path}")
    plt.close(fig)


def plot_field_snapshots(histories: dict, outdir: str):
    """Plot damage field snapshots at key steps for each temperature."""
    if not HAS_MPL:
        return

    for T in sorted(histories.keys()):
        d = histories[T]

        # Check for mesh and snapshots
        if 'mesh_nodes' not in d or 'mesh_elems' not in d:
            continue

        nodes = d['mesh_nodes']
        elems = d['mesh_elems']
        steps = d['snapshot_steps'] if 'snapshot_steps' in d else np.array([])

        if len(steps) == 0:
            continue

        # Pick up to 4 snapshots spread across the simulation
        if len(steps) <= 4:
            pick = steps
        else:
            idx = np.linspace(0, len(steps) - 1, 4, dtype=int)
            pick = steps[idx]

        import matplotlib.tri as mtri
        tri = mtri.Triangulation(nodes[:, 0] * 1e3, nodes[:, 1] * 1e3, elems)

        fig, axes = plt.subplots(1, len(pick), figsize=(4 * len(pick), 4))
        if len(pick) == 1:
            axes = [axes]

        for i, step in enumerate(pick):
            key = f'd_step_{step}'
            if key not in d:
                continue
            damage = d[key]
            ax = axes[i]
            tc = ax.tripcolor(tri, damage, cmap='hot', vmin=0, vmax=1,
                              shading='flat')
            ax.set_aspect('equal')
            ax.set_title(f'Step {step}', fontsize=10)
            ax.set_xlabel('x (mm)')
            if i == 0:
                ax.set_ylabel('y (mm)')

        fig.suptitle(f'Damage field — T = {T:.0f} K', fontsize=12)
        fig.colorbar(tc, ax=axes, shrink=0.6, label='damage $d$')

        path = os.path.join(outdir, f'damage_snapshots_{T:.0f}K.png')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved damage snapshots for T={T:.0f}K")


def print_summary_table(summary: dict):
    """Print a formatted summary table."""
    T = summary['T']
    print(f"\n{'='*90}")
    print(f"  {'T':>6}  {'Gc_eff':>8}  {'KJ_dom':>8}  {'KJ_glob':>8}  "
          f"{'Fmax':>8}  {'Wp/Wext':>8}  {'rho_max':>10}  {'Gc_max':>8}")
    print(f"  {'(K)':>6}  {'(J/m²)':>8}  {'(MPa√m)':>8}  {'(MPa√m)':>8}  "
          f"{'(kN)':>8}  {'(%)':>8}  {'(m⁻²)':>10}  {'(J/m²)':>8}")
    print(f"  {'-'*86}")

    for i in range(len(T)):
        Wp_r = min(summary['Wp_ratio'][i], 99999)
        print(f"  {T[i]:6.0f}  {summary['Gc_eff'][i]:8.2f}  "
              f"{summary['KJ_peak'][i]:8.3f}  {summary['KJ_global_peak'][i]:8.3f}  "
              f"{summary['Fmax'][i]/1e3:8.1f}  {Wp_r:8.1f}  "
              f"{summary['rho_max'][i]:10.2e}  {summary['Gc_max'][i]:8.1f}")
    print(f"{'='*90}")


def main():
    parser = argparse.ArgumentParser(description='Post-process fracture results')
    parser.add_argument('--dir', default='results_arrhenius',
                        help='Results directory')
    parser.add_argument('--snapshots', action='store_true',
                        help='Also plot damage field snapshots')
    args = parser.parse_args()

    results_dir = args.dir
    if not os.path.isdir(results_dir):
        print(f"Directory not found: {results_dir}")
        return

    print(f"Loading results from {results_dir}/")
    histories = load_all_histories(results_dir)
    if not histories:
        return

    summary = extract_summary(histories)
    print_summary_table(summary)

    outdir = os.path.join(results_dir, 'plots')
    os.makedirs(outdir, exist_ok=True)

    plot_summary(histories, summary, outdir)

    if args.snapshots:
        plot_field_snapshots(histories, outdir)

    print(f"\nDone. Plots saved to {outdir}/")


if __name__ == '__main__':
    main()
