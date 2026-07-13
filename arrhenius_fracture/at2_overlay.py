"""
AT2 OVERLAY MODEL  (parallel to the sharp-front engine)  -- PROTOTYPE
================================================================================
Branching, fragmentation, and spontaneous nucleation (AT2's strengths) WITHOUT
giving up the Arrhenius first-passage kinetics (DBTT/ceramic/peak/weak-T regime
physics) the sharp engine encodes.

DESIGN (agreed):
  * Two representations BOTH live, overlaid:
      - AT2 damage field phi  -> owns GEOMETRY: path, kinking, branching, merging,
        nucleus growth. Never carries a tip radius; never "sharpens".
      - per-tip sharp engines -> own KINETICS: each detected tip carries its own
        FrontEngine (ledger N_em, blunting r_eff, shielding, first-passage clock).
        Single source of truth for the RATE of energy release.
  * DRIVING FORCE per tip = isolated K (reduced/J-integral) + PAIRWISE elastic
    interaction (Kachanov, leading order) -> K_eff,i. Collinear shields, parallel
    amplifies: competing-hazard physics taken from elasticity, not invented.
  * LABELED LOCAL PURSES (budget handshake): each tip authorizes an advance da_auth
    from ITS OWN clock, ENERGY = Gc*da_auth, spendable ONLY in a ~L_pz neighbourhood
    of that tip. A branch SPLITS its parent's purse (no global pool). AT2 sets the
    SHAPE/direction; the ledger sets the MAGNITUDE.
  * PROPAGATION vs NUCLEATION is a financed, geometric decision: surface inside a
    tip purse and within budget = propagation; surface AT2 wants where no purse
    reaches (orphan / ahead-of-tip) is referred to a BULK nucleation hazard.
  * CONSERVATION GATE (type A, ledger-summed): created propagation surface energy
    == ledger authorization:  residual = Gc*dA_prop - sum_i Gc*da_auth,i  -> 0.
    Nucleation is a SEPARATE channel with its own (bulk-hazard) accounting.

REDUCED IN THIS PROTOTYPE (flagged again at call sites):
  [FEM]    reduced analytic K instead of FEM + configurational force. Hook:
           DrivingForce.k_iso  ->  replace with fem.py + j_integral.py.
  [ITER]   one pass per load step; the Delta2 fixed-point (re-solve K after AT2
           moves the geometry / compliance changed) is stubbed. See AT2Config.n_inner.
  [PZMERGE]pairwise elastic K used up to zone contact; AT2 label-merge handles
           coalescence; the smooth handoff is crude.
  [3BODY]  pairwise interaction only; degrades in dense knots.
  [KERNEL] isotropic purse kernel; anisotropic (forward-biased) deferred.
  [PLAST]  no spatial bulk plastic field; the per-tip engine still carries its
           internal shielding/embrittlement ledger, so tip KINETICS are intact.
First-cut choices; revisit if they bite.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from scipy import ndimage

from .sharp_front import build_engine, _build_parser
from .config import make_emergent_config


# ============================================================================ #
@dataclass
class Microstructure:
    """Smooth random Gc(x)/strength(x). Weak (low-Gc) patches site nucleation and
    bias branch selection. [no GBs / no cleavage anisotropy yet -- deferred]."""
    Gc0: float = 8.0
    Gc_amp: float = 0.35
    corr_len_cells: float = 6.0
    seed: int = 0

    def build(self, shape):
        rng = np.random.default_rng(self.seed)
        g = ndimage.gaussian_filter(rng.standard_normal(shape),
                                    self.corr_len_cells, mode='reflect')
        g /= (g.std() + 1e-12)
        Gc = np.clip(self.Gc0 * (1.0 + self.Gc_amp * g), 0.2 * self.Gc0, 3.0 * self.Gc0)
        weak = self.Gc0 / Gc
        return Gc, weak


# ============================================================================ #
@dataclass
class Tip:
    iy: int
    ix: int
    nx_dir: float
    ny_dir: float
    a_eff: float
    comp: int
    engine: object = None
    K_iso: float = 0.0
    K_eff: float = 0.0
    da_auth: float = 0.0


def extract_tips(phi, h, thresh=0.5):
    """Tips from the phi field, skimage-free and reproducible. One tip per crack
    component, located along the component's PCA PRINCIPAL AXIS (its own elongation
    direction), at the interior end furthest along that axis. Growth direction = the
    principal axis oriented toward that end. This is stable: the crack grows along
    its length, so the advance makes net forward progress and cannot wander
    sideways or stall. [multi-tip/true branching approximated -- see TODO]"""
    dmg = phi > thresh
    if not dmg.any():
        return [], np.zeros_like(phi, int)
    lab, n = ndimage.label(dmg, structure=np.ones((3, 3)))
    ny, nx = phi.shape
    tips = []
    min_cells = 6
    for comp in range(1, n + 1):
        ys, xs = np.nonzero(lab == comp)
        if len(xs) < min_cells:
            continue
        a_eff = h * np.hypot(np.ptp(xs) + 1, np.ptp(ys) + 1)
        cyc, cxc = ys.mean(), xs.mean()
        xm, ym = xs - cxc, ys - cyc
        cov = np.array([[(xm * xm).mean(), (xm * ym).mean()],
                        [(xm * ym).mean(), (ym * ym).mean()]])
        evals, evecs = np.linalg.eigh(cov)
        ax, ay = evecs[:, int(np.argmax(evals))]          # principal axis (x, y)
        proj = xm * ax + ym * ay
        interior = (xs > 1) & (xs < nx - 2) & (ys > 1) & (ys < ny - 2)
        if not interior.any():
            continue
        # the interior end with the largest |projection| is the tip; orient axis to it
        order = np.argsort(proj)
        cand = [order[-1], order[0]]                       # +axis end, -axis end
        tip_k = None
        for c in cand:
            if interior[c]:
                tip_k = c
                break
        if tip_k is None:
            tip_k = int(np.argmax(np.where(interior, np.abs(proj), -1.0)))
        iy, ix = int(ys[tip_k]), int(xs[tip_k])
        # orient growth direction along axis, toward the tip
        if (ix - cxc) * ax + (iy - cyc) * ay < 0:
            ax, ay = -ax, -ay
        tips.append(Tip(iy, ix, float(ax), float(ay), float(a_eff), int(comp)))
    return tips, lab


# ============================================================================ #
class DrivingForce:
    """Isolated K [FEM: replace k_iso with J-integral on a contour >> ell] plus the
    leading-order PAIRWISE elastic interaction (Kachanov): collinear shields,
    parallel amplifies. Supplies the competing-hazard physics AND the overlap
    weighting for the purse partition."""

    def __init__(self, Eprime, h):
        self.Eprime, self.h = Eprime, h

    def k_iso(self, tip, sigma_app):
        # [FEM HOOK] reduced edge-crack K = 1.12 sigma sqrt(pi a).
        return 1.12 * sigma_app * np.sqrt(np.pi * max(tip.a_eff, self.h))

    def _interaction_factor(self, ti, tj):
        rx = (tj.ix - ti.ix) * self.h
        ry = (tj.iy - ti.iy) * self.h
        d = np.hypot(rx, ry)
        if d < 1e-12:
            return 0.0
        cos = float(np.clip((rx * ti.nx_dir + ry * ti.ny_dir) / d, -1, 1))
        ang = (1.0 - 3.0 * cos * cos)                 # +1 parallel (amplify), -2 collinear (shield)
        strength = 0.5 * max(tj.a_eff, self.h) / d    # near-field 1/d  [3BODY: pairwise only]
        return float(np.clip(ang * strength, -0.9, 2.0))

    def k_eff(self, tips, sigma_app):
        for ti in tips:
            ti.K_iso = self.k_iso(ti, sigma_app)
        for i, ti in enumerate(tips):
            corr = sum(self._interaction_factor(ti, tj) for j, tj in enumerate(tips)
                       if j != i and tj.comp != ti.comp)   # only BETWEEN cracks
            corr = float(np.clip(corr, -0.7, 3.0))   # total shielding cannot annihilate K_eff
            ti.K_eff = max(ti.K_iso * (1.0 + corr), 0.0)
        return tips


# ============================================================================ #
class Nucleation:
    """Bulk Poisson hazard lambda_nuc(x)=nu0*weak*exp(-dG(sigma,T)/kT); cumulative
    Theta(x)+=lambda*dt; fires (seeds a phi nucleus) when Theta>=1. Sites the crack
    inventory on weak/stressed features. Also the channel that licenses ORPHAN AT2
    surface (ahead-of-tip / outside every purse)."""

    def __init__(self, shape, nu0=3e7, dG0_eV=1.55, sigma_star_GPa=6.0):
        self.Theta = np.zeros(shape)
        self.nu0, self.dG0_eV, self.sig_star = nu0, dG0_eV, sigma_star_GPa * 1e9

    def rate(self, sigma_field, weak, T):
        kB = 8.617e-5
        dG = self.dG0_eV * np.maximum(1.0 - sigma_field / self.sig_star, 0.05)
        r = self.nu0 * np.sqrt(weak) * np.exp(-dG / (kB * T))  # weak modulates, not dominates
        # HARD stress gate: nucleation is a high-stress process; below a real stress
        # fraction the hazard is zero (no slow accumulation at weak low-stress spots,
        # which produced the spurious early corner cracks).
        r[sigma_field < 0.45 * self.sig_star] = 0.0
        return r

    def advance(self, sigma_field, weak, T, dt, exclude):
        if self.nu0 <= 0:
            return []
        self.Theta += self.rate(sigma_field, weak, T) * dt
        self.Theta[exclude] = 0.0
        fired = self.Theta >= 1.0
        sites = []
        if fired.any():
            lab, n = ndimage.label(fired)
            for c in range(1, n + 1):
                ys, xs = np.nonzero(lab == c)
                k = int(np.argmax(self.Theta[ys, xs]))
                sites.append((int(ys[k]), int(xs[k])))
                self.Theta[ys, xs] = 0.0
        return sites


# ============================================================================ #
@dataclass
class AT2Config:
    Lx: float = 1.0e-3
    Ly: float = 1.0e-3
    nx: int = 160
    ny: int = 160
    ell: float = 2.0e-5
    L_pz: float = 1.0e-6
    purse_radius_mult: float = 8.0     # purse radius in units of ell (resolvable on
    #   the phi mesh). The PHYSICAL L_pz enters the tip hazard via the engine.
    T: float = 900.0
    n_inner: int = 1                   # [ITER] Delta2 fixed-point passes (stub=1)


class AT2Overlay:
    def __init__(self, cfg, micro, engine_args=None, nucleation=None):
        self.cfg = cfg
        self.h = cfg.Lx / cfg.nx
        self.shape = (cfg.ny, cfg.nx)
        mc = make_emergent_config()
        self.mat, self.Eprime = mc.material, mc.material.Eprime
        self.phi = np.zeros(self.shape)
        self.phi_max = np.zeros(self.shape)
        self.Gc, self.weak = micro.build(self.shape)
        self._Gc_gy, self._Gc_gx = np.gradient(self.Gc)   # for meander steering
        self.drive = DrivingForce(self.Eprime, self.h)
        self.nucl = nucleation or Nucleation(self.shape)
        self._engine_args = engine_args or self._default_engine_args()
        self._engines = []                 # spatial list: {eng,iy,ix,used,last}
        self.energy_log, self.inventory_log = [], []
        # initial edge notch
        cy, a0 = cfg.ny // 2, int(0.20 * cfg.nx)
        self.phi[cy - 1:cy + 1, :a0] = 1.0
        self.phi_max[:] = self.phi

    def _default_engine_args(self):
        a = _build_parser().parse_args(['--mode', '1d', '--temperatures', '900'])
        a.emit_S_T_c0_kB, a.emit_S_T_c1, a.emit_S_sigma_max_kB = -20.0, 0.02, 8.0
        a.multihit_m, a.multihit_tau = 3.0, 1e-6
        a.emb_sat_frac, a.cleave_H0_eV = 1.0, 3.0
        a.chi_shield, a.N_sat, a.v_rayleigh = 0.6, 2000.0, 2600.0
        return a

    def _surface_energy(self):
        gy, gx = np.gradient(self.phi, self.h)
        dens = self.phi ** 2 / (2 * self.cfg.ell) + 0.5 * self.cfg.ell * (gx ** 2 + gy ** 2)
        return float(np.sum(self.Gc * dens) * self.h ** 2)

    def _match_engine(self, tip, step_idx, match_cells):
        """Persist, per PHYSICAL crack, a FrontEngine ledger AND a FIXED growth axis,
        by nearest-neighbour matching across steps (component labels are unstable).
        The axis is set once at birth and never recomputed, so the forward direction
        cannot drift as the crack evolves (axis drift was the propagation-stall
        cause). Returns the record dict (carries 'eng' and 'axis')."""
        best, bestd = None, match_cells ** 2
        for rec in self._engines:
            if rec['used']:
                continue
            d2 = (rec['iy'] - tip.iy) ** 2 + (rec['ix'] - tip.ix) ** 2
            if d2 < bestd:
                best, bestd = rec, d2
        if best is None:
            eng = build_engine(self._engine_args, self.mat)
            eng.f.da = max(self.h, 2e-6)
            best = {'eng': eng, 'iy': tip.iy, 'ix': tip.ix, 'used': True,
                    'last': step_idx, 'axis': (tip.nx_dir, tip.ny_dir)}
            self._engines.append(best)
        else:
            best.update(used=True, iy=tip.iy, ix=tip.ix, last=step_idx)
        return best

    def step(self, sigma_app, dt, step_idx=0):
        cfg = self.cfg
        for rec in self._engines:
            rec['used'] = False

        tips, _ = extract_tips(self.phi, self.h)
        auth_total = dS_prop = resid = rel = auth_len = realized_len = 0.0
        if tips:
            match_cells = cfg.purse_radius_mult * cfg.ell / self.h
            for t in tips:
                rec = self._match_engine(t, step_idx, match_cells)
                t.engine = rec['eng']
                t.nx_dir, t.ny_dir = rec['axis']    # FIXED axis -> forward never drifts
            self.drive.k_eff(tips, sigma_app)          # K_iso + pairwise interaction
            # per-step advance cap: don't step over the AT2 instability. Once the
            # clock completes (n_fire can run to 1e7 = "ligament severs in one step"
            # in the sharp model), we cap the FIELD advance to a few ell and let an
            # unstable crack keep advancing over subsequent steps. Also bounded by
            # the Rayleigh speed. This is the Delta-t cap; total energy is unchanged,
            # it is just metered out over steps so the conservation gate can close.
            adv_cap = min(4.0 * cfg.ell, t.engine.f.v_rayleigh * dt
                          if np.isfinite(t.engine.f.v_rayleigh) else np.inf)
            for t in tips:                              # AUTHORIZE from first-passage clock
                info = t.engine.step(t.K_eff, cfg.T, dt)
                t.da_auth = min(t.engine.f.da * info['n_fire'], adv_cap)
            # PROPAGATION handshake (conservation gate): advance = ledger, shape = AT2
            S0 = self._surface_energy()
            realized_len = self._advance_fronts(tips)
            dS_prop = self._surface_energy() - S0       # AT2 FUNCTIONAL change (diagnostic;
            #   relates to Gc*length by the AT2 normalization constant c_w -- [TODO] calibrate)
            auth_len = sum(t.da_auth for t in tips)     # authorized advance LENGTH
            # CONSERVATION GATE (length form): realized advance == authorized advance.
            # Energy on both sides is Gc*length; they match iff lengths match. Branch
            # splits the parent purse, so branching cannot inflate realized_len.
            auth_total = sum(self.Gc[t.iy, t.ix] * t.da_auth for t in tips)  # energy, for reporting
            resid = realized_len - auth_len
            rel = abs(resid) / (abs(auth_len) + 1e-30)

        self._engines = [r for r in self._engines if step_idx - r['last'] <= 3]

        # NUCLEATION: separate channel, own accounting. Exclude existing cracks AND
        # every tip's purse (~ purse_radius): inside a purse, growth is propagation
        # (engine-authorized), so nucleation may only fire on stressed weak features
        # AWAY from the tips -- which is the spontaneous-nucleation observable.
        Sn0 = self._surface_energy()
        # exclude existing cracks AND a process-zone margin around them: nucleation
        # fires only in the genuine far field at weak features, never in the K-field
        # halo (which feeds tip propagation, not independent cracks).
        dmg = self.phi > 0.5
        margin_px = int(round(2.0 * cfg.purse_radius_mult * cfg.ell / self.h))
        exclude = dmg.copy() if not dmg.any() else (ndimage.distance_transform_edt(~dmg) <= margin_px)
        b = 16                                  # no nucleation in the boundary band
        exclude[:b, :] = exclude[-b:, :] = exclude[:, :b] = exclude[:, -b:] = True
        sites = self.nucl.advance(self._sigma_proxy(sigma_app, tips), self.weak,
                                  cfg.T, dt, exclude)
        for (iy, ix) in sites:
            self._seed_nucleus(iy, ix)
        dS_nucl = self._surface_energy() - Sn0

        self.energy_log.append({'sigma_app': sigma_app, 'dS_prop': dS_prop,
                                'auth_ledger': auth_total, 'prop_residual': resid,
                                'prop_rel_residual': rel, 'dS_nucl': dS_nucl,
                                'auth_len': auth_len, 'realized_len': realized_len,
                                'n_tips': len(tips), 'n_nucleated': len(sites)})
        lab, ncr = ndimage.label(self.phi > 0.5, structure=np.ones((3, 3)))
        sizes = ndimage.sum(np.ones_like(lab), lab, range(1, ncr + 1)) if ncr else []
        self.inventory_log.append({'n_cracks': int(ncr),
                                   'total_damaged_area': float((self.phi > 0.5).sum() * self.h ** 2),
                                   'sizes_cells': [int(s) for s in sizes]})
        return self.energy_log[-1], self.inventory_log[-1]

    def _advance_fronts(self, tips):
        """Advance each active tip forward by da_auth (ledger MAGNITUDE) along the
        AT2-preferred direction (SHAPE), then relax to width ell. Microstructure can
        split the advance into two lobes when both flanks are weaker -> emergent
        branching. [ITER] direction frozen within the step (n_inner=1)."""
        advanced = False
        realized = 0.0                  # total advance length actually painted [m]
        for t in tips:
            if t.da_auth <= 0.0:
                continue
            ncells = max(int(round(t.da_auth / self.h)), 1)
            dirx, diry = t.nx_dir, t.ny_dir
            perp = (-diry, dirx)
            fy, fx = t.iy + diry, t.ix + dirx
            g0 = self._sample_Gc(fy, fx)
            gA = self._sample_Gc(fy + perp[1], fx + perp[0])
            gB = self._sample_Gc(fy - perp[1], fx - perp[0])
            split = (gA < 0.85 * g0) and (gB < 0.85 * g0)     # weak flanks -> branch
            if split:
                headings = [(dirx + 0.6 * perp[0], diry + 0.6 * perp[1]),
                            (dirx - 0.6 * perp[0], diry - 0.6 * perp[1])]
                ncells_path = max(ncells // 2, 1)   # BRANCH SPLITS THE PARENT PURSE
            else:
                headings = [(dirx, diry)]
                ncells_path = ncells
            steer = 0.0                  # MEANDER OFF: the heuristic meander stalls on
            #   some microstructures (path must come from energy minimization, i.e.
            #   the FEM-driven AT2 evolution -- see TODO). Straight advance is robust.
            d0 = (t.nx_dir, t.ny_dir)    # FIXED forward axis (always present -> no stall)
            perp = (-d0[1], d0[0])
            for (hx0, hy0) in headings:
                cxf, cyf = float(t.ix), float(t.iy)
                for _k in range(ncells_path):
                    iyc = int(np.clip(round(cyf), 0, self.cfg.ny - 1))
                    ixc = int(np.clip(round(cxf), 0, self.cfg.nx - 1))
                    # MEANDER = bounded LATERAL nudge toward lower Gc, recomputed fresh
                    # from the fixed axis each cell (NOT accumulated), so the forward
                    # component always dominates and the crack cannot curl back / stall.
                    glat = self._Gc_gx[iyc, ixc] * perp[0] + self._Gc_gy[iyc, ixc] * perp[1]
                    gnorm = abs(self._Gc_gx[iyc, ixc]) + abs(self._Gc_gy[iyc, ixc]) + 1e-30
                    lat = -np.clip(glat / gnorm, -1.0, 1.0)          # toward lower Gc
                    hx = hx0 + steer * lat * perp[0]
                    hy = hy0 + steer * lat * perp[1]
                    hn = np.hypot(hx, hy) + 1e-12
                    hx, hy = hx / hn, hy / hn
                    cxf += hx; cyf += hy
                    py, px = int(round(cyf)), int(round(cxf))
                    if 0 <= py < self.cfg.ny and 0 <= px < self.cfg.nx:
                        if self.phi[py, px] < 1.0:
                            realized += self.h
                        self.phi[py, px] = 1.0
                        advanced = True
                        advanced = True
        if advanced:
            self.phi_max = np.maximum(self.phi_max, self.phi)   # lock in the advance FIRST
            for _ in range(3):
                self._at2_relax()                                # then regularize edges
            self.phi_max = np.maximum(self.phi_max, self.phi)
        return realized

    def _sample_Gc(self, iy, ix):
        return self.Gc[int(np.clip(iy, 0, self.cfg.ny - 1)), int(np.clip(ix, 0, self.cfg.nx - 1))]

    def _at2_relax(self):
        """Allen-Cahn relaxation toward the AT2 profile (regularizes the painted
        crack to width ell, irreversible). The SHAPE operator; no external drive."""
        lap = (np.roll(self.phi, 1, 0) + np.roll(self.phi, -1, 0)
               + np.roll(self.phi, 1, 1) + np.roll(self.phi, -1, 1) - 4 * self.phi) / self.h ** 2
        dtau = 0.2 * self.h ** 2 / self.cfg.ell ** 2
        self.phi = np.clip(self.phi - dtau * (self.phi - self.cfg.ell ** 2 * lap),
                           self.phi_max, 1.0)

    def _sigma_proxy(self, sigma_app, tips):
        """Reduced stress field: far-field sigma_app PLUS a near-tip K-field
        concentration sigma ~ K_eff/sqrt(2 pi r) ahead of each tip, so the tip is
        the dominant stress concentrator (otherwise nucleation fires uniformly on
        the far field and the tip never wins). [FEM] replace with the FEM principal
        stress field -- this proxy is the single biggest approximation and is what
        makes the demo's tip/nucleation balance only qualitative."""
        field = np.full(self.shape, sigma_app)
        yy, xx = np.mgrid[0:self.cfg.ny, 0:self.cfg.nx]
        for t in tips:
            r = np.hypot((xx - t.ix) * self.h, (yy - t.iy) * self.h)
            conc = t.K_eff / np.sqrt(2 * np.pi * np.maximum(r, 0.5 * self.h))
            taper = np.exp(-r / (0.5 * t.a_eff + 5 * self.cfg.ell))   # local to the tip
            field = np.maximum(field, sigma_app + conc * taper)
        # near-tip stress saturates at the cohesive/theoretical strength -- it is NOT
        # singular. Cap removes the spurious high-stress ring just outside the purse.
        field = np.minimum(field, 12.0e9)
        return field * (0.8 + 0.4 * self.weak)

    def _seed_nucleus(self, iy, ix, r=1):
        sl = (slice(max(iy - r, 0), iy + r + 1), slice(max(ix - r, 0), ix + r + 1))
        self.phi[sl] = 1.0
        self.phi_max[sl] = 1.0


# ============================================================================ #
def run_demo(out='at2_overlay_out', steps=34, sigma0=0.30e9, dsigma=0.04e9,
             nx=160, ny=160, T=650.0, dt=20.0, seed=2, nucleation_on=False, render=True,
             Lx=2.0e-3):
    import os
    os.makedirs(out, exist_ok=True)
    cfg = AT2Config(nx=nx, ny=ny, T=T, Lx=Lx, Ly=Lx, ell=2.5e-5)
    model = AT2Overlay(cfg, Microstructure(seed=seed),
                       nucleation=Nucleation((ny, nx), nu0=(3e7 if nucleation_on else 0.0)))
    print("=" * 74)
    print("  AT2 OVERLAY prototype  (sharp ledger authorizes, AT2 field spends)")
    print(f"  grid {nx}x{ny}  h={model.h*1e6:.2f}um  ell={cfg.ell*1e6:.2f}um "
          f"(={cfg.ell/model.h:.1f} cells)  T={T}K  nucleation {'ON' if nucleation_on else 'OFF'}")
    print("=" * 74)
    snaps = []
    for s in range(1, steps + 1):
        sigma_app = sigma0 + dsigma * s
        e, inv = model.step(sigma_app, dt=dt, step_idx=s)
        if s % max(1, steps // 8) == 0 or s == 1:
            snaps.append((s, sigma_app, model.phi.copy()))
            print(f"  step {s:3d} sig={sigma_app/1e9:5.2f}GPa cracks={inv['n_cracks']:2d} "
                  f"tips={e['n_tips']:2d} nucl={e['n_nucleated']} "
                  f"auth={e['auth_ledger']:8.2e} dS_prop={e['dS_prop']:8.2e} "
                  f"rel_resid={e['prop_rel_residual']:.3f}")
    rr = np.array([e['prop_rel_residual'] for e in model.energy_log if e['auth_ledger'] > 0])
    print("-" * 74)
    if rr.size:
        print(f"  PROPAGATION conservation residual (steps with advance): "
              f"median={np.median(rr):.3f}  p90={np.percentile(rr, 90):.3f}  n={rr.size}")
    else:
        print("  (no authorized advance occurred -- raise load or check K vs Kc(T))")
    print(f"  final crack inventory: {model.inventory_log[-1]['n_cracks']} cracks")
    if render:
        _render(out, model, snaps)
    return model


def _render(out, model, snaps):
    import os
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = max(len(snaps), 1)
    fig, axes = plt.subplots(2, n, figsize=(2.8 * n, 6.0))
    axes = np.atleast_2d(axes)
    if axes.shape[0] == 1:
        axes = axes.reshape(2, -1)
    for k, (s, sig, phi) in enumerate(snaps):
        axes[0, k].imshow(phi, origin='lower', cmap='inferno', vmin=0, vmax=1)
        axes[0, k].set_title(f'step {s}\n$\\sigma$={sig/1e9:.2f} GPa', fontsize=9)
        axes[0, k].set_xticks([]); axes[0, k].set_yticks([])
    axes[1, 0].imshow(model.Gc, origin='lower', cmap='viridis')
    axes[1, 0].set_title('Gc(x) microstructure', fontsize=9)
    axes[1, 0].set_xticks([]); axes[1, 0].set_yticks([])
    el = model.energy_log
    st = np.arange(1, len(el) + 1)
    if axes.shape[1] > 1:
        ax = axes[1, 1]
        ax.plot(st, [e.get('auth_len', 0)*1e6 for e in el], '--', label='ledger authorized', lw=1.8)
        ax.plot(st, [e.get('realized_len', 0)*1e6 for e in el], label='AT2 realized', lw=1.3)
        ax.set_title('advance handshake [um]', fontsize=9); ax.legend(fontsize=7)
        ax.set_xlabel('step', fontsize=8)
    if axes.shape[1] > 2:
        axes[1, 2].plot(st, [iv['n_cracks'] for iv in model.inventory_log], color='firebrick')
        axes[1, 2].set_title('crack inventory', fontsize=9); axes[1, 2].set_xlabel('step', fontsize=8)
    for k in range(3, axes.shape[1]):
        axes[1, k].axis('off')
    fig.tight_layout()
    p = os.path.join(out, 'at2_overlay.png')
    fig.savefig(p, dpi=120); plt.close(fig)
    print(f"  saved {p}")


# ============================================================================ #
#  Temperature sweep CLI                                                        #
# ============================================================================ #
def run_temperature_sweep(temperatures, out='at2_sweep', steps=34,
                          sigma0=0.25e9, dsigma=0.025e9, nx=160, ny=160, dt=20.0,
                          Lx=2.0e-3, ell=2.5e-5, seed=2, nucleation_on=False,
                          nu0=3e7, dG0_eV=1.55, sigma_star_GPa=6.0, render=True):
    """Run the AT2 overlay at a series of temperatures with everything else fixed,
    and summarise the response (final crack count, nucleation count, how far the
    main crack ran, conservation residual). Per-T phi snapshots are saved."""
    import os, json
    os.makedirs(out, exist_ok=True)
    rows = []
    per_T_final = []
    for T in temperatures:
        cfg = AT2Config(nx=nx, ny=ny, T=T, Lx=Lx, Ly=Lx, ell=ell)
        nucl = Nucleation((ny, nx), nu0=(nu0 if nucleation_on else 0.0),
                          dG0_eV=dG0_eV, sigma_star_GPa=sigma_star_GPa)
        m = AT2Overlay(cfg, Microstructure(seed=seed), nucleation=nucl)
        cy = ny // 2
        for s in range(1, steps + 1):
            m.step(sigma0 + dsigma * s, dt=dt, step_idx=s)
        xs = np.nonzero(m.phi[cy - 1:cy + 2, :].max(0) > 0.5)[0]
        tip_mm = (xs.max() * m.h * 1e3) if len(xs) else 0.0
        rr = [e['prop_rel_residual'] for e in m.energy_log if e.get('realized_len', 0) > 0]
        row = {'T_K': T,
               'main_tip_mm': round(float(tip_mm), 3),
               'final_cracks': int(m.inventory_log[-1]['n_cracks']),
               'total_nucleations': int(sum(e['n_nucleated'] for e in m.energy_log)),
               'damaged_area_mm2': round(m.inventory_log[-1]['total_damaged_area'] * 1e6, 4),
               'resid_median': round(float(np.median(rr)) if rr else 0.0, 3)}
        rows.append(row)
        per_T_final.append((T, m.phi.copy()))
        print(f"  T={T:6.1f}K  main_tip={row['main_tip_mm']:.2f}mm  "
              f"cracks={row['final_cracks']:2d}  nucl={row['total_nucleations']:2d}  "
              f"area={row['damaged_area_mm2']:.3f}mm^2  resid={row['resid_median']:.3f}")
    with open(os.path.join(out, 'sweep_summary.json'), 'w') as f:
        json.dump(rows, f, indent=2)
    if render:
        _render_sweep(out, per_T_final, rows)
    return rows


def _render_sweep(out, per_T_final, rows):
    import os
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = len(per_T_final)
    fig, axes = plt.subplots(2, max(n, 1), figsize=(2.8 * max(n, 1), 5.6))
    axes = np.atleast_2d(axes)
    if axes.shape[0] == 1:
        axes = axes.reshape(2, -1)
    for k, (T, phi) in enumerate(per_T_final):
        axes[0, k].imshow(phi, origin='lower', cmap='inferno', vmin=0, vmax=1)
        axes[0, k].set_title(f'{T:.0f} K', fontsize=9)
        axes[0, k].set_xticks([]); axes[0, k].set_yticks([])
    T = [r['T_K'] for r in rows]
    ax = axes[1, 0]
    ax.plot(T, [r['main_tip_mm'] for r in rows], 'o-', label='main tip [mm]')
    ax.plot(T, [r['final_cracks'] for r in rows], 's--', label='crack count')
    ax.plot(T, [r['total_nucleations'] for r in rows], '^:', label='nucleations')
    ax.set_xlabel('T [K]', fontsize=8); ax.legend(fontsize=7)
    ax.set_title('response vs temperature', fontsize=9)
    for k in range(1, axes.shape[1]):
        axes[1, k].axis('off')
    fig.tight_layout()
    p = os.path.join(out, 'at2_sweep.png')
    fig.savefig(p, dpi=120); plt.close(fig)
    print(f"  saved {p}")


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='AT2 overlay prototype')
    ap.add_argument('--mode', choices=['demo', 'sweep'], default='demo')
    ap.add_argument('--temperatures', type=float, nargs='+',
                    default=[500, 650, 800, 950, 1100])
    ap.add_argument('--out', default='at2_sweep')
    ap.add_argument('--steps', type=int, default=34)
    ap.add_argument('--nx', type=int, default=160)
    ap.add_argument('--ny', type=int, default=160)
    ap.add_argument('--Lx', type=float, default=2.0e-3)
    ap.add_argument('--ell', type=float, default=2.5e-5)
    ap.add_argument('--dt', type=float, default=20.0)
    ap.add_argument('--sigma0', type=float, default=0.25e9)
    ap.add_argument('--dsigma', type=float, default=0.025e9)
    ap.add_argument('--seed', type=int, default=2)
    ap.add_argument('--nucleation-on', action='store_true',
                    help='EXPERIMENTAL: secondary nucleation (placement qualitative, raises conservation residual -- needs the FEM driving force)')
    ap.add_argument('--nu0', type=float, default=1e8)
    ap.add_argument('--dG0-eV', type=float, default=1.7)
    ap.add_argument('--sigma-star-GPa', type=float, default=6.0)
    ap.add_argument('--no-render', action='store_true')
    a = ap.parse_args()
    if a.mode == 'demo':
        run_demo(out=a.out, steps=a.steps, nx=a.nx, ny=a.ny, dt=a.dt,
                 sigma0=a.sigma0, dsigma=a.dsigma, seed=a.seed, Lx=a.Lx,
                 nucleation_on=a.nucleation_on, render=not a.no_render)
    else:
        run_temperature_sweep(a.temperatures, out=a.out, steps=a.steps, nx=a.nx,
                              ny=a.ny, dt=a.dt, sigma0=a.sigma0, dsigma=a.dsigma,
                              Lx=a.Lx, ell=a.ell, seed=a.seed,
                              nucleation_on=a.nucleation_on, nu0=a.nu0,
                              dG0_eV=a.dG0_eV, sigma_star_GPa=a.sigma_star_GPa,
                              render=not a.no_render)
